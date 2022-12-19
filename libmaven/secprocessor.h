#pragma once

#ifndef SECPROCESSOR_H
#define SECPROCESSOR_H

#include "mzUtils.h"
#include "mzSample.h"

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
      * @param csimilarityFractionDiffTol: maximum number of fractions off where peaks can be matched.
      * ========================*/
     int similarityMinNumPeaks = 2;
     int similarityFractionDiffTol = 1;

     /** =======================
      * Peak similarity scoring
      * comments:
      * ========================*/
     int peakSimMaxCenterDiff = 0;

     string encodeParams();
     shared_ptr<SECSearchParameters> static decode(string encodedParams);
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

    string id; //unique name-type string

    SECTraceType type = SECTraceType::Unset;

    vector<int> fractionNums{}; //Includes missing fractions, e.g. from params

    vector<float> rawIntensities{}; //For missing intensities, use traceMissingIntensityFill from params
    vector<float> smoothedIntensities{}; //based on parameters

    shared_ptr<SECSearchParameters> params;

    vector<Peak> peaks{};

    Fragment *fragment = nullptr;

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

#endif // SECPROCESSOR_H

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

    int getMinFractionNum();
    int getMaxFractionNum();
    string getPeakComparisonId();

    void printSummary();

    SECTracePeakComparison(SECTrace *first, int firstPeakNum, SECTrace *second, int secondPeakNum);

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
