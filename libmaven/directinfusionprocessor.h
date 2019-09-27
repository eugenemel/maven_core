#pragma once

#include "mzSample.h"
#include "mzUtils.h"
#include <memory>

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

/**
 * @brief The DirectInfusionSearchSet class
 * Data container class
 */
class DirectInfusionSearchSet {

public:

        /**
         *  set of all DIA ranges represented in experiment.
         */
        set<int> mapKeys = {};

        /**
         * mapping of each DIA range to actual m/z values
         * <first, second> = <minMz, maxMz>
         */
        map<int, pair<float,float>> mzRangesByMapKey = {};

        /**
         * <key, value> = <map_key, valid pair<Compound,Adduct>>
         */
        multimap<int, pair<Compound*, Adduct*>> compoundsByMapKey = {};

};

class DirectInfusionProcessor {

    public:
    /**
          * @brief getSearchSet
          * @param sample
          * a representative sample - may be anything
          * @param compounds
          * @param adducts
          * @param isRequireAdductPrecursorMatch
          * @param debug
          * @return DirectInfusionSearchSet
          * --> all compound, adducts, organized into m/z bins.
          *
          * This structure can be reused if all samples in an experiment have the same organization.
          */
         static shared_ptr<DirectInfusionSearchSet> getSearchSet(
                 mzSample *sample,
                 const vector<Compound*>& compounds,
                 const vector<Adduct*>& adducts,
                 bool isRequireAdductPrecursorMatch,
                 bool debug);

         static vector<DirectInfusionAnnotation*> processSingleSample(
                 mzSample *sample,
                 shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
                 bool debug);
};


