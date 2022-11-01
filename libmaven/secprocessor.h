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
     * ========================*/
     float traceMissingIntensityFill = 0.0f;
     int traceMinFractionNumber = 1;
     int traceMaxFractionNumber = 72;

     /** =======================
      * PEAK PICKING
      * comments:
      * baseline smoothing is the same as peak smoothing.
      * @param traceMinPeakIntensity and @param traceMinPeakSN are compared to the raw intensities, not smoothed values
      * ========================*/
     EIC::SmootherType traceSmoothingType = EIC::SmootherType::GAUSSIAN;
     int traceWindowSize = 5;
     float traceMinPeakIntensity = 1e6f;
     float traceMinPeakSN = 0.0f;
     int traceBaselineDropTopX = 80; //EIC is full width
     float tracePeakBoundsMaxIntensityFraction = 0.1f;

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

    SECTraceType type = SECTraceType::Unset;

    vector<int> fractionNums{}; //Includes missing fractions, e.g. from params

    vector<float> rawIntensities{}; //For missing intensities, use traceMissingIntensityFill from params
    vector<float> smoothedIntensities{}; //based on parameters

    shared_ptr<SECSearchParameters> params;

    vector<Peak> peaks{};

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
};

#endif // SECPROCESSOR_H
