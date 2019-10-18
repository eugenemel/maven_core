#pragma once

#include "mzSample.h"
#include "mzUtils.h"
#include <memory>

class mzSample;
class DirectInfusionAnnotation;
enum class SpectralCompositionAlgorithm;

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
 * @brief The SpectralDeconvolutionAlgorithm enum
 *
 * Short description of different approaches for spectral deconvolution
 * (goal is to determine the relative proportions of different compounds)
 */
enum class SpectralCompositionAlgorithm {
    ALL_CANDIDATES,
    MEDIAN_UNIQUE
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
    int minNumMatches = 5;

    /**
     * @brief minNumUniqueMatches
     * minimum number of matches for a single <Compound*, Adduct*>
     * with unique fragment m/zs, given the universe of all <Compound*, Adduct*>
     * matches searched.
     *
     * Considered after @param minNumMatches? - the idea is that @param minNumMatches is used
     * to find likely IDs, and @param minNumUniqueMatches might be mroe useful for determining
     * relative composition
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

    /**
     * @brief spectralCompositionAlgorithm
     * By default, do nothing, just return all matches, without doing any elimination or quantitation
     * of spectral components.
     */
    SpectralCompositionAlgorithm spectralCompositionAlgorithm = SpectralCompositionAlgorithm::ALL_CANDIDATES;

};

/**
 * @brief The DirectInfusionMatchData struct
 *
 * A container class for organizing association data
 */
struct DirectInfusionMatchData {
    Compound* compound;
    Adduct* adduct;
    FragmentationMatchScore fragmentationMatchScore;
    double proportion = 0;
};

/**
 * @brief The DirectInfusionMatchInformation structure
 *
 * A structure to organize all fragment matches from all compound, adduct pairs that match to a single
 * direct infusion spectrum.
 *
 */
struct DirectInfusionMatchInformation {

public:
    map<int, vector<shared_ptr<DirectInfusionMatchData>>> fragToMatchData = {};
    map<shared_ptr<DirectInfusionMatchData>, vector<int>> matchDataToFrags = {};
    map<pair<int, shared_ptr<DirectInfusionMatchData>>,float> fragToTheoreticalIntensity = {};
    map<pair<int, shared_ptr<DirectInfusionMatchData>>,float> fragToObservedIntensity = {};

    float getNormalizedTheoreticalIntensity(int fragId, shared_ptr<DirectInfusionMatchData> matchData){return fragToTheoreticalIntensity.at(make_pair(fragId, matchData));}
    float getObservedIntensity(int fragId, shared_ptr<DirectInfusionMatchData> matchData){return fragToObservedIntensity.at(make_pair(fragId, matchData));}

    float getIntensityRatio(int fragId, shared_ptr<DirectInfusionMatchData> matchData){
        return (getObservedIntensity(fragId, matchData) / getNormalizedTheoreticalIntensity(fragId, matchData));
    }
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
     static map<int, DirectInfusionAnnotation*> processSingleSample(
             mzSample *sample,
             shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
             shared_ptr<DirectInfusionSearchParameters> params,
             bool debug);

     /**
      * @brief deconvolveAllShared
      * @param allCandidates
      * @return
      *
      * Triggered by SpectralDeconvolutionAlgorithm::ALL_SHARED
      *
      * Removes compounds with all shared fragments, computes relative abundance based on
      * unshared fragments.
      *
      * Input is the list of all candidates, plus the observed spectrum they all matched to
      */
     static vector<shared_ptr<DirectInfusionMatchData>> determineComposition(
             vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
             Fragment *observedSpectrum,
             SpectralCompositionAlgorithm algorithm,
             bool debug
             );

     /**
      * @brief getMatches
      * @param allCandidates
      *
      * @return
      * function to return all compound matches organized into maps, either with key as compound
      * or fragment m/z.
      *
      * fragments converted m/z <--> int keys using mzToIntKey(mz, 1000) and intKeyToMz(intKey, 1000).
      *
      * Note that this function does no processing, filtering, or analysis - it simply reorganizes
      * the compound match data into maps.
      */
     static unique_ptr<DirectInfusionMatchInformation> getMatchInformation(
             vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
             Fragment *observedSpectrum,
             bool debug);
};

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
     * @brief precMzMin, precMzMax
     * refers to the m/z of precursors.
     */
    float precMzMin;
    float precMzMax;

    /**
     * @brief sample
     * source sample
     */
    mzSample* sample = nullptr;

    /**
     * @brief scan
     * a single representative scan from the sample.
     */
    Scan *scan = nullptr;

    /**
     * @brief fragmentationPattern
     * Agglomeration of multiple DI scans (if they exist),
     * otherwise same data as 'scan'.
     */
    Fragment* fragmentationPattern;

    /**
     * compound, adduct, and estimated proportion of the spectrum
     * associated with the match.
     *
     * FragmentationMatchScores are also provided.
     */
    vector<shared_ptr<DirectInfusionMatchData>> compounds;
};

/**
 * @brief The DirectInfusionGroupAnnotation class
 * group of DirectInfusionAnnotation results across many samples.
 */
class DirectInfusionGroupAnnotation : DirectInfusionAnnotation {

public:

    /**
     * @brief annotationBySample
     * retain original samples for reference.
     */
    map<mzSample*, DirectInfusionAnnotation*> annotationBySample = {};

    void clean();
    static unique_ptr<DirectInfusionGroupAnnotation> createByAverageProportions(vector<DirectInfusionAnnotation*> singleSampleAnnotations);

};

typedef map<int, vector<shared_ptr<DirectInfusionMatchData>>>::iterator fragToMatchDataIterator;
typedef map<shared_ptr<DirectInfusionMatchData>, vector<int>>::iterator matchDataToFragIterator;
typedef map<shared_ptr<DirectInfusionMatchData>, vector<float>>::iterator matchDataToFragIntensityIterator;
typedef map<shared_ptr<DirectInfusionMatchData>, float>::iterator matchDataToFloatIterator;


