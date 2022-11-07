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
     float traceMinPeakSN = 0.0f;
     int traceBaselineDropTopX = 80; //EIC is full width
     float tracePeakBoundsMaxIntensityFraction = 0.1f;

     /** =======================
      * FRAGMENT
      * comments:
      * Parameters associated with =eneration of Fragment from SECTrace
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
     * @param type
     * @param fractionNums
     * @param rawIntensities
     * @param params
     */
    SECTrace(SECTraceType type,
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

    virtual string getName();
    virtual float getSimilarity(SECTrace* first, SECTrace* second, shared_ptr<SECSearchParameters> params, bool debug = false) = 0;

    virtual ~SECTraceSimilarity(){}
};

class SECTraceSimilarityCosine : public SECTraceSimilarity {
public:
    string getName(){return "cosine";}
    float getSimilarity(SECTrace* first, SECTrace* second, shared_ptr<SECSearchParameters> params, bool debug = false);
    virtual ~SECTraceSimilarityCosine(){}
};

#endif // SECPROCESSOR_H
