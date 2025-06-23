#pragma once

#ifndef SECPROCESSOR_H
#define SECPROCESSOR_H

#include "mzUtils.h"
#include "mzSample.h"

class SECTraceSimilarityCosine;

class SECSearchParameters {

    public:

    /** =======================
     * PROGRAM
     * @param searchVersion: version of search protocol used to generate results.
     * ========================*/
    string searchVersion = "UNKNOWN";

    /** =======================
     * SEC TRACE filling
     * comments:
     * baseline smoothing is the same as peak smoothing.
     * @param traceNormalizeToSumIntensity is carried out after filling.
     * ========================*/
     float traceMissingIntensityFill = 0.0f;
     int traceMinFractionNumber = 1;
     int traceMaxFractionNumber = 72;
     bool traceNormalizeToSumIntensity = false;

     /** =======================
      * PEAK PICKING
      * comments:
      * baseline smoothing is the same as peak smoothing.
      * @param traceMinPeakIntensity and @param traceMinPeakSN are compared to the raw intensities, not smoothed values
      * ========================*/
     EIC::SmootherType traceSmoothingType = EIC::SmootherType::GAUSSIAN;
     int traceWindowSize = 5;
     float traceMinPeakIntensity = 1e6f;
     float traceMinSmoothedIntensity = 1e6f;
     float traceMinFracTopPeakIntensity = 0.0f; // fraction of max intensity
     float traceMinFracTopSmoothedIntensity = 0.0f; // fraction of max intensity
     float traceMinPeakSN = 0.0f;
     int traceMinPeakWidth = 0; // # of fractions in peak
     int traceBaselineDropTopX = 80; //EIC is full width
     float tracePeakBoundsMaxIntensityFraction = 0.0f;
     float traceRtBoundsSlopeThreshold = 0.01f; //minimum change in slope between peak points, as a fraction of the peak's max intensity (otherwise, reached edge of peak)
     bool traceIsPickEdgePeaks = false; // treat one-sided maxima at the edge of the sec trace as peaks


     /** =======================
      * PEAK Grouping
      * ========================*/
     float groupMaxFracDiff = 1.0f; // analogous to PeakPickingAndGroupingParameters.groupMaxRtDiff
     float groupMergeOverlap = 0.8f; // analogous to PeakPickingAndGroupingParameters.groupMergeOverlap
     bool groupIsMergeOverlappingPeakGroups = true; // analogous to PeakPickingAndGroupingParameters.groupIsMergeOverlappingPeakGroups
     int groupMinNumPeaks = 1;

     /** =======================
      * FRAGMENT
      * comments:
      * Parameters associated with generation of Fragment from SECTrace
      * Used primarily for similarity metrics.
      * @param fragmentIsSmoothedIntensity: use smoothed intensity for peak values.
      * ========================*/
     bool fragmentIsSmoothedIntensity = true;

     /** =======================
      * SIMILARITY SCORING
      * comments:
      * parameters specific to similarity scoring computations.
      * @param similarityMinNumPeaks: minimum number of peaks required in each SECTrace to perform comparison.
      * @param similarityFractionDiffTol: maximum number of fractions off where peaks can be matched.
      * ========================*/
     int similarityMinNumPeaks = 2;
     int similarityFractionDiffTol = 1;

     /** =======================
      * Peak similarity scoring
      * comments:
      * parameters used to compute/return peak similarity metric.
      * Once a comparison filter fails, the SECTracePeakSimilarity immediately stops being built.
      * SECTracePeakScorer checks a SECTracePeakSimilarity flag to decide to keep the SECTracePeakSimilarity or not.
      * ========================*/
     int peakSimMaxCenterDiff = 0;
     float peakSimMinSecFractionOverlap = 0.0f;
     float peakSimMinSecFractionJaccard = 0.0f;
     float peakSimMinSmoothedCorrelation = 0.0f;
     float peakSimMinRawCorrelation = 0.0f;
     int peakSimMinFractionNumber = -1;
     int peakSimMaxFractionNumber = -1;

     string encodeParams();
     shared_ptr<SECSearchParameters> static decode(string encodedParams);
     shared_ptr<PeakPickingAndGroupingParameters> toPeakPickingAndGroupingParams();
};

enum SECTraceType {Unset=0, Peptide=1, Protein=2};

/**
 * @brief The SECTrace class
 *
 * Contains raw data as well as processed data (interpolated points, smoothed points, peaks, etc).
 *
 * This is all kept together because in practice, every time a SECTrace is made,
 * data should be smoothed and peaks should be picked (no need to separate these responsibilities).
 */
class SECTrace {

public:

    //unique name-type string for this trace.
    //this ID must always be unique amongst all analytes and samples.
    string id;

    // name-type string for this trace.
    // This ID describes the analyte associated with this trace.
    // This ID may be shared between different samples / sample roll-ups.
    // For example, this could be a protein ID, gene name, or peptide sequence.
    string analyteId;

    //name-type string for this trace.
    //This ID describes a biological source for this trace such as a sample or group.
    //There may be multiple analyteIds for a given biologicalId
    //For example, the same sample might have multiple proteins.
    string biologicalId;

    SECTraceType type = SECTraceType::Unset;

    vector<int> fractionNums{}; //Includes missing fractions, e.g. from params

    vector<float> rawIntensities{}; //For missing intensities, use traceMissingIntensityFill from params
    vector<float> smoothedIntensities{}; //based on parameters

    shared_ptr<SECSearchParameters> params;

    vector<Peak> peaks{};

    Fragment *fragment = nullptr;
    EIC *eic = nullptr;

    //default constructor
    SECTrace(){}

    /**
     * @brief SECTrace
     * @param id
     * @param type
     * @param fractionNums
     * @param rawIntensities
     * @param params
     * @param debug
     */
    SECTrace(string id,
             SECTraceType type,
             vector<int> fractionNums,
             vector<float> rawIntensities,
             shared_ptr<SECSearchParameters> params,
             bool debug = false);

    vector<string> getPeakSummaryString(
            string empty = "",
            string leftPrefix = "L",
            string maxPrefix = "M",
            string rightPrefix = "R");

    Fragment* getFragment(
            shared_ptr<SECSearchParameters> params,
            bool debug = false);

    void computeTraceData(vector<int> fractionNums,
                          vector<float> rawIntensities,
                          shared_ptr<SECSearchParameters> params,
                          bool debug = false);

    void pickPeaks(bool debug);

    string getPeakPositionsString();
};

class SECTraceGroups {
public:
    string id;
    vector<SECTrace*> secTraces;
    shared_ptr<SECSearchParameters> params;

    vector<PeakGroup> groups{};
    vector<mzSample*> samples{};

    void computePeakGroups(bool debug);

    vector<int> getGroupIdsVector(SECTrace* trace, unsigned long groupIdOffset = 0);
    vector<float> getGroupsQuantVector(SECTrace* trace, string quantType);
};

//Introduced in 622
class SECTraceDiff : public SECTrace {

public:

    //Additional fields
    SECTrace *compare;
    SECTrace *reference;

    //compare - reference intensities
    vector<float> diffRawIntensities{};
    vector<float> diffSmoothedIntensities{};

    //additional comparisons
    SECTraceSimilarityCosine *similarityScore;

    //peaks, fragment, and absolute intensities are fields inherited by SECTrace

    SECTraceDiff(SECTrace *compare, SECTrace *reference, bool debug=false);
};

class SECTraceSimilarity {
public:

    float similarity = -1.0f; //unset

    virtual string getName() = 0;
    virtual float getSimilarity(bool debug = false);

    virtual ~SECTraceSimilarity(){}
};

class SECTraceSimilarityCosine : public SECTraceSimilarity {
public:

    SECTrace *first;
    SECTrace *second;
    shared_ptr<SECSearchParameters> params;
    string compareId;
    vector<int> ranks;

    //scores
    //float similarity = -1.0f; //inherited (value=matchedPeakCosineScore)
    int numPeakMatches = 0;
    float cosineScore = -1.0f;
    float matchedPeakCosineScore = -1.0f;
    float fractionPeaksMatched = 0.0f;
    float pearsonCorrelationRaw = -1.0f;
    float pearsonCorrelationSmoothed = -1.0f;

    SECTraceSimilarityCosine(SECTrace* first, SECTrace* second, shared_ptr<SECSearchParameters> params);

    string getName(){return "cosine";}
    float getSimilarity(bool debug = false);
    virtual ~SECTraceSimilarityCosine(){}
};

class SECTraceCosineSimilarityScorer {
public:
    static vector<SECTraceSimilarityCosine> scoreTraces(vector<SECTrace*> traces,
                                                        shared_ptr<SECSearchParameters> params,
                                                        bool debug = false);
};

class SECTracePeak {
public:
    SECTrace *trace=nullptr;
    int peakNum = -1;

    int getPeakFractionNum();
    int getMinFractionNum();
    int getMaxFractionNum();

    // position according to peak intensity vectors (getSmoothedIntensities(), getRawIntensities())
    // use peak.pos for position according to trace intensity vector
    int getPeakIndex();

    //comparableRange is <number left of peakIndex, number right of peakIndex>.
    //if unspecified, return the whole range (from peak min to max)
    vector<float> getSmoothedIntensities(pair<int, int> comparableRange = make_pair(-1, -1));
    vector<float> getRawIntensities(pair<int, int> comparableRange = make_pair(-1, -1));
    vector<int> getFractionNums(pair<int, int> comparableRange = make_pair(-1, -1));

    string getPeakId();

    SECTracePeak(SECTrace *trace, int peakNum);
    SECTracePeak(){}

private:
    bool isValid();
    vector<float> getComparableRangeSubset(const vector<float>&x, pair<int, int> comparableRange);
    vector<int> getComparableRangeSubset(const vector<int>&x, pair<int, int> comparableRange);
};

class SECTracePeakComparison {
public:
    SECTracePeak first;
    SECTracePeak second;

    // <number left of pivot, number right of pivot>
    pair<int, int> comparableRange = make_pair(-1, -1);

    float pearsonCorrelationSmoothed = -1.0f;
    float pearsonCorrelationRaw = -1.0f;

    float secFractionOverlap = -1.0f;
    float secFractionJaccard = -1.0f;

    int peakCenterDistance = -1;

    bool isPassesParameterFilters = false;

    int getMinFractionNum();
    int getMaxFractionNum();
    string getPeakComparisonId();

    void printSummary();

    SECTracePeakComparison(SECTrace *first, int firstPeakNum, SECTrace *second, int secondPeakNum, shared_ptr<SECSearchParameters> params);

private:
    void computeComparableRange();
    void computeSecFractionJaccard();

};

class SECTracePeakScorer {
public:
    static vector<SECTracePeakComparison> scorePeaks(
            vector<SECTrace*> traces,
            shared_ptr<SECSearchParameters> params,
            bool debug);
};

class SECTraceDiffGenerator {
public:
    static vector<SECTraceDiff*> generateSECTraceDiffs(
            vector<SECTrace*> referenceTraces,
            vector<SECTrace*> compareTraces,
            bool debug);
};

#endif // SECPROCESSOR_H
