#pragma once

#include "mzSample.h"
#include "mzUtils.h"

class mzSample;

/**
 * @brief The DirectInfusionAnnotation class
 *
 * MS/MS scans from an @param sample are agglomerated,
 * and compared to a compound database to identify matches,
 * and relative abundance of various compounds in those scans.
 *
 * Essential to the annotation are the @param precMzMin and @param precMzMax
 * fields, which describe the m/z range scanned for this annotation.
 */
class DirectInfusionAnnotation {

public:

        /**
         * @brief sample
         * source sample
         */
        mzSample* sample = nullptr;

        /**
         * @brief precMzMin, precMzMax
         * refers to the m/z of precursors.
         */
        float precMzMin;
        float precMzMax;

        /**
         * @brief scan
         * a representative scan.
         * Potentially, an average of many scans collected over different infusion times.
         * Or, may just be a scan taken right out of the mzSample.
         */
        Scan *scan = nullptr;

        /**
         * @brief fragmentationPattern and fragMatchScore
         * fragmentation data. used ultimately for peak group display.
         *
         * TODO: refactor to pointers for speed / storage,
         * but need to worry about memory leaks
         */
        Fragment* fragmentationPattern;
        FragmentationMatchScore fragMatchScore;

        /**
         * each tuple refers to the compound, adduct, and estimated proportion of the spectrum
         * associated with the match.
         *
         * FragmentationMatchScores are also provided.
         */
        vector<tuple<Compound*, Adduct*, double, FragmentationMatchScore>> compounds;
};

class DirectInfusionProcessor {

    public:
         static std::vector<DirectInfusionAnnotation*> processSingleSample(mzSample *sample, const vector<Compound*>& compounds, const vector<Adduct*>& adducts, bool debug);
};


