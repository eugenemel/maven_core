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

/**
 * @brief The DirectInfusionSearchParameters class
 *
 * single class to contain all parameters used in direct infusion search analysis.
 */
class DirectInfusionSearchParameters {

public:

    /**
     * @brief minNumMatches
     * mininum number of matches for a single <Compound*, Adduct*>
     * to match to a spectrum in order to retain this <Compound*, Adduct*>
     * as a component of the observed spectrum
     */
    int minNumMatches = 3;

    /**
     * @brief minNumUniqueMatches
     * minimum number of matches for a single <Compound*, Adduct*>
     * with unique fragment m/zs, given the universe of all <Compound*, Adduct*>
     * matches searched.
     *
     * TODO: should this be considered before or after @param minNumMatches?
     * Could use @param minNumMatches as a first filter, then this as a subsequent filter.
     */
    int minNumUniqueMatches = 0;

    /**
     * @brief isRequireAdductPrecursorMatch
     * The compound's associated adduct must be the adduct in the supplied list of adducts,
     * otherwise the match will be ignored.
     */
    bool isRequireAdductPrecursorMatch = true;

    /**
     * @brief productPpmTolr
     * tolerance value used for matching library fragment m/z s to Scan m/z s
     */
    float productPpmTolr = 20;

};

/**
 * @brief The DirectInfusionProcessor class
 * All methods should be static - functional programming paradigm
 */
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
                 shared_ptr<DirectInfusionSearchParameters> params,
                 bool debug);

         /**
          * @brief processSingleSample
          * @param sample
          * @param directInfusionSearchSet
          * @param debug
          * @return
          *
          * Returns DirectInfusionAnnotation assessments for a single sample.
          * TODO: think about how to agglomerate these across samples?
          * What to do when there are different compositions in different sample?
          * eg, sample 1 has 70% A, 20% B, 10% C, and sample 2 and 50% A, 0% B, 0% C, and 50% D?
          *
          * Definitely some choices to be made here
          */
         static vector<DirectInfusionAnnotation*> processSingleSample(
                 mzSample *sample,
                 shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
                 shared_ptr<DirectInfusionSearchParameters> params,
                 bool debug);
};


