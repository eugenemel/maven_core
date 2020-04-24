#pragma once

#include "mzSample.h"
#include "mzUtils.h"
#include <memory>
#include <algorithm>

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
     map<int, vector<pair<Compound*, Adduct*>>> compoundsByMapKey = {};

};

/**
 * @brief The SpectralDeconvolutionAlgorithm enum
 *
 * Short description of different approaches for spectral deconvolution
 * (goal is to determine the relative proportions of different compounds)
 */
enum class SpectralCompositionAlgorithm {
    ALL_CANDIDATES,
    AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE
};

/**
 * @brief The DirectInfusionSearchParameters class
 *
 * single class to contain all parameters used in direct infusion search analysis.
 */
class DirectInfusionSearchParameters {

public:

    /** =======================
     * SCAN FILTER ASSOCIATED
     * All parameters are arguments of Fragment::Fragment() constructor.
     * ========================*/

    float scanFilterMinFracIntensity = 0;
    float scanFilterMinSNRatio = 0;
    unsigned int scanFilterMaxNumberOfFragments = UINT_MAX;
    int scanFilterBaseLinePercentile = 0;
    bool scanFilterIsRetainFragmentsAbovePrecursorMz = true;
    float scanFilterPrecursorPurityPpm = 0;
    float scanFilterMinIntensity = 0;

    /** =======================
     * CONSENSUS SPECTRUM ASSOCIATED
     * All parameters are arguments Fragment::buildConsensus() method.
     * ========================*/

    float consensusPpmTolr = 10;
    bool consensusIsIntensityAvgByObserved = true;
    int consensusMinNumMs2Scans = 0;
    float consensusMinFractionMs2Scans = 0;
    bool consensusIsNormalizeTo10K = false;

    /** ===================
     * MS2 SEARCH RELATED
     * @param ms2MinNumMatches: min number of reference peaks found in observed spectrum
     * @param ms2MinNumDiagnosticMatches: min number of reference peaks labeled as diagnostic found in observed spectrum
     * @param ms2MinNumUniqueMatches: CURRENTLY UNUSED PARAMETER
     * @param ms2PpmTolr: m/z tolerance value used for matching reference <--> observed spectra
     * @param ms2MinIntensity: minimum intensity value for an MS2 spectral peak to be considered real
     * ==================== */

    int ms2MinNumMatches = 5;
    int ms2MinNumDiagnosticMatches = 0;
    int ms2MinNumUniqueMatches = 0;
    float ms2PpmTolr = 20;
    float ms2MinIntensity = 0;

    /** ===================
     * MS1 SEARCH RELATED
     * @param ms1IsRequireAdductPrecursorMatch:reference compound associated adduct must == query adduct
     * @param ms1IsFindPrecursorIon: only retain matches where precursor peak is found in MS1 scan(s).
     * @param ms1PpmTolr: tolerance value used for matching theoretical ion m/z to an m/z peak in an MS1 scan
     * @param ms1MinIntensity: min intensity value for a MS1 spectral peak to be considered real
     * @param ms1ScanFilter: consider only MS1 scans that substring match in their filterString field to this
     * ==================== */

    bool ms1IsRequireAdductPrecursorMatch = true;
    bool ms1IsFindPrecursorIon = false;
    float ms1PpmTolr = 5;
    float ms1MinIntensity = 0;
    string ms1ScanFilter = "";

    /** ===================
     * AGGLOMERATION
     * @param isAgglomerateAcrossSamples: If true, align results across all samples. Otherwise, return each sample results individually.
     * @param spectralCompositionAlgorithm:
     *      SpectralCompositionAlgorithm::ALL_CANDIDATES: Return all matches without any elimination or quantitation
     *      SpectralCompositionAlgorithm::AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE:
     *          automatically summarize results to higher level if possible.
     *          Apply parsimony to spectral matches.
     * ==================== */

    bool isAgglomerateAcrossSamples = false;
    SpectralCompositionAlgorithm spectralCompositionAlgorithm = SpectralCompositionAlgorithm::ALL_CANDIDATES;

    void printParams(){
        string encodedParams = encodeParams();
        replace(encodedParams.begin(), encodedParams.end(), ';', '\n');
        replace(encodedParams.begin(), encodedParams.end(), '=', ' ');
        cout << encodedParams << endl;
    }

    string encodeParams() {

        string encodedParams;

        //scan filter params
        encodedParams = encodedParams + "scanFilterMinFracIntensity" + "=" + to_string(scanFilterMinFracIntensity) + ";";
        encodedParams = encodedParams + "scanFilterMinSNRatio" + "=" + to_string(scanFilterMinSNRatio) + ";";
        encodedParams = encodedParams + "scanFilterMaxNumberOfFragments" + "=" + to_string(scanFilterMaxNumberOfFragments) + ";";
        encodedParams = encodedParams + "scanFilterBaseLinePercentile" + "=" + to_string(scanFilterBaseLinePercentile) + ";";
        encodedParams = encodedParams + "scanFilterIsRetainFragmentsAbovePrecursorMz" + "=" + to_string(scanFilterIsRetainFragmentsAbovePrecursorMz) + ";";
        encodedParams = encodedParams + "scanFilterPrecursorPurityPpm" + "=" + to_string(scanFilterPrecursorPurityPpm) + ";";
        encodedParams = encodedParams + "scanFilterMinIntensity" + "=" + to_string(scanFilterMinIntensity) + ";";

        //consensus spectrum params
        encodedParams = encodedParams + "consensusPpmTolr" + "=" + to_string(consensusPpmTolr) + ";";
        encodedParams = encodedParams + "consensusIsIntensityAvgByObserved" + "=" + to_string(consensusIsIntensityAvgByObserved) + ";";
        encodedParams = encodedParams + "consensusMinNumMs2Scans" + "=" + to_string(consensusMinNumMs2Scans) + ";";
        encodedParams = encodedParams + "consensusMinFractionMs2Scans" + "=" + to_string(consensusMinFractionMs2Scans) + ";";
        encodedParams = encodedParams + "consensusIsNormalizeTo10K" + "=" + to_string(consensusIsNormalizeTo10K) + ";";

        //ms2 search params
        encodedParams = encodedParams + "ms2MinNumMatches" + "=" + to_string(ms2MinNumMatches) + ";";
        encodedParams = encodedParams + "ms2MinNumDiagnosticMatches" + "=" + to_string(ms2MinNumDiagnosticMatches) + ";";
        encodedParams = encodedParams + "ms2MinNumUniqueMatches" + "=" + to_string(ms2MinNumUniqueMatches) + ";";
        encodedParams = encodedParams + "ms2PpmTolr" + "=" + to_string(ms2PpmTolr) + ";";
        encodedParams = encodedParams + "ms2MinIntensity" + "=" + to_string(ms2MinIntensity) + ";";

        //ms1 search params
        encodedParams = encodedParams + "ms1IsRequireAdductPrecursorMatch" + "=" + to_string(ms1IsRequireAdductPrecursorMatch) + ";";
        encodedParams = encodedParams + "ms1IsFindPrecursorIon" + "=" + to_string(ms1IsFindPrecursorIon) + ";";
        encodedParams = encodedParams + "ms1PpmTolr" + "=" + to_string(ms1PpmTolr) + ";";
        encodedParams = encodedParams + "ms1MinIntensity" + "=" + to_string(ms1MinIntensity) + ";";
        encodedParams = encodedParams + "ms1ScanFilter" + "=" + ms1ScanFilter + ";";

        //agglomeration params
        encodedParams = encodedParams + "isAgglomerateAcrossSamples" + "=" + to_string(isAgglomerateAcrossSamples) + ";";

        string spectralCompositionAlgorithmStr = "UNSPECIFIED";
        if (spectralCompositionAlgorithm == SpectralCompositionAlgorithm::ALL_CANDIDATES) {
            spectralCompositionAlgorithmStr = "ALL_CANDIDATES";
        } else if (spectralCompositionAlgorithm == SpectralCompositionAlgorithm::AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE) {
            spectralCompositionAlgorithmStr = "AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE";
        }
        encodedParams = encodedParams + "spectralCompositionAlgorithm" + "=" + spectralCompositionAlgorithmStr + ";";

        return encodedParams;
    }

    static shared_ptr<DirectInfusionSearchParameters> decode(string encodedParams){
        shared_ptr<DirectInfusionSearchParameters> directInfusionSearchParameters = shared_ptr<DirectInfusionSearchParameters>(new DirectInfusionSearchParameters());

        unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams);

        //scan filter params
        if (decodedMap.find("scanFilterMinFracIntensity") != decodedMap.end()){
            directInfusionSearchParameters->scanFilterMinFracIntensity = stof(decodedMap["scanFilterMinFracIntensity"]);
        }
        if (decodedMap.find("scanFilterMinSNRatio") != decodedMap.end()){
            directInfusionSearchParameters->scanFilterMinSNRatio = stof(decodedMap["scanFilterMinSNRatio"]);
        }
        if (decodedMap.find("scanFilterMaxNumberOfFragments") != decodedMap.end()) {
            string scanFilterMaxNumberOfFragmentsStr = decodedMap["scanFilterMaxNumberOfFragments"];
            try {
                directInfusionSearchParameters->scanFilterMaxNumberOfFragments = stoul(decodedMap["scanFilterMaxNumberOfFragments"]);
            } catch (std::out_of_range exception ){}
        }
        if (decodedMap.find("scanFilterBaseLinePercentile") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterBaseLinePercentile = stoi(decodedMap["scanFilterBaseLinePercentile"]);
        }
        if (decodedMap.find("scanFilterIsRetainFragmentsAbovePrecursorMz") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterIsRetainFragmentsAbovePrecursorMz = decodedMap["scanFilterIsRetainFragmentsAbovePrecursorMz"] == "1";
        }
        if (decodedMap.find("scanFilterPrecursorPurityPpm") != decodedMap.end()){
            directInfusionSearchParameters->scanFilterPrecursorPurityPpm = stof(decodedMap["scanFilterPrecursorPurityPpm"]);
        }
        if (decodedMap.find("scanFilterMinIntensity") != decodedMap.end()){
            directInfusionSearchParameters->scanFilterMinIntensity = stof(decodedMap["scanFilterMinIntensity"]);
        }

        //consensus spectrum params
        if (decodedMap.find("consensusPpmTolr") != decodedMap.end()){
            directInfusionSearchParameters->consensusPpmTolr = stof(decodedMap["consensusPpmTolr"]);
        }
        if (decodedMap.find("consensusIsIntensityAvgByObserved") != decodedMap.end()){
            directInfusionSearchParameters->consensusIsIntensityAvgByObserved = decodedMap["consensusIsIntensityAvgByObserved"] == "1";
        }
        if (decodedMap.find("consensusMinNumMs2Scans") != decodedMap.end()){
            directInfusionSearchParameters->consensusMinNumMs2Scans = stoi(decodedMap["consensusMinNumMs2Scans"]);
        }
        if (decodedMap.find("consensusMinFractionMs2Scans") != decodedMap.end()){
            directInfusionSearchParameters->consensusMinFractionMs2Scans = stof(decodedMap["consensusMinFractionMs2Scans"]);
        }
        if (decodedMap.find("consensusIsNormalizeTo10K") != decodedMap.end()){
            directInfusionSearchParameters->consensusIsNormalizeTo10K = decodedMap["consensusIsNormalizeTo10K"] == "1";
        }

        //ms2 search params
        if (decodedMap.find("ms2MinNumMatches") != decodedMap.end()){
            directInfusionSearchParameters->ms2MinNumMatches = stoi(decodedMap["ms2MinNumMatches"]);
        }
        if (decodedMap.find("ms2MinNumDiagnosticMatches") != decodedMap.end()){
            directInfusionSearchParameters->ms2MinNumDiagnosticMatches = stoi(decodedMap["ms2MinNumDiagnosticMatches"]);
        }
        if (decodedMap.find("ms2MinNumUniqueMatches") != decodedMap.end()){
            directInfusionSearchParameters->ms2MinNumUniqueMatches = stoi(decodedMap["ms2MinNumUniqueMatches"]);
        }
        if (decodedMap.find("ms2PpmTolr") != decodedMap.end()){
            directInfusionSearchParameters->ms2PpmTolr = stof(decodedMap["ms2PpmTolr"]);
        }
        if (decodedMap.find("ms2MinIntensity") != decodedMap.end()){
            directInfusionSearchParameters->ms2MinIntensity = stof(decodedMap["ms2MinIntensity"]);
        }

        //ms1 search params
        if (decodedMap.find("ms1IsRequireAdductPrecursorMatch") != decodedMap.end()){
            directInfusionSearchParameters->ms1IsRequireAdductPrecursorMatch = decodedMap["ms1IsRequireAdductPrecursorMatch"] == "1";
        }
        if (decodedMap.find("ms1IsFindPrecursorIon") != decodedMap.end()){
            directInfusionSearchParameters->ms1IsFindPrecursorIon = decodedMap["ms1IsFindPrecursorIon"] == "1";
        }
        if (decodedMap.find("ms1PpmTolr") != decodedMap.end()){
            directInfusionSearchParameters->ms1PpmTolr = stof(decodedMap["ms1PpmTolr"]);
        }
        if (decodedMap.find("ms1MinIntensity") != decodedMap.end()){
            directInfusionSearchParameters->ms1MinIntensity = stof(decodedMap["ms1MinIntensity"]);
        }
        if (decodedMap.find("ms1ScanFilter") != decodedMap.end()){
            directInfusionSearchParameters->ms1ScanFilter = decodedMap["ms1ScanFilter"];
        }

        //agglomeration params
        if (decodedMap.find("isAgglomerateAcrossSamples") != decodedMap.end()){
            directInfusionSearchParameters->isAgglomerateAcrossSamples = decodedMap["isAgglomerateAcrossSamples"] == "1";
        }
        if (decodedMap.find("spectralCompositionAlgorithm") != decodedMap.end()){
            string spectralCompositionAlgorithmStr = decodedMap[""];
            if (spectralCompositionAlgorithmStr == "ALL_CANDIDATES") {
                directInfusionSearchParameters->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::ALL_CANDIDATES;
            } else if (spectralCompositionAlgorithmStr == "ALL_CANDIDATES") {
                directInfusionSearchParameters->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE;
            }
        }

        return directInfusionSearchParameters;
    }
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
    float fragmentMaxObservedIntensity = 0;
    double proportion = 0;
};

/**
 * @brief The DirectInfusionMatchDataCompare struct
 *
 * lhs->compound == rhs->compound iff lhs and rhs point to the same data, even if memory addresses are different
 *
 * Special class for comparisons
 */
struct DirectInfusionMatchDataCompare {
    bool operator() (const shared_ptr<DirectInfusionMatchData>& lhs, const shared_ptr<DirectInfusionMatchData>& rhs) const {
        if (lhs->compound == rhs->compound) {
            return lhs->adduct < rhs->adduct;
        } else {
            return lhs->compound < rhs->compound;
        }
    }
};

struct DirectInfusionMatchDataCompareByNames {
    bool operator() (const shared_ptr<DirectInfusionMatchData>& lhs, const shared_ptr<DirectInfusionMatchData>& rhs) const {
        if (lhs->compound && rhs->compound) {
            if (lhs->compound->name == rhs->compound->name) {
                if (lhs->adduct && rhs->adduct) {
                    return lhs->adduct->name < rhs->adduct->name;
                } else {
                    return false;
                }
            } else {
                return lhs->compound->name < rhs->compound->name;
            }
        } else {
            return false;
        }
    }
};

/**
 * @brief The DirectInfusionSinglePeakMatchData struct
 *
 * container for peak match data, used by @link DirectInfusionMatchInformation
 */
struct DirectInfusionSinglePeakMatchData {
    float normalizedTheoreticalIntensity;
    float observedIntensity;
    float getIntensityRatio() {return (observedIntensity / normalizedTheoreticalIntensity); }
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

    //observed
    map<pair<int, shared_ptr<DirectInfusionMatchData>>,float> fragToObservedIntensity = {};

    //unsummarized
    map<int, vector<shared_ptr<DirectInfusionMatchData>>> fragToMatchData = {};
    map<shared_ptr<DirectInfusionMatchData>, vector<int>> matchDataToFrags = {};
    map<pair<int, shared_ptr<DirectInfusionMatchData>>,float> fragToTheoreticalIntensity = {};

    /**
      * If
      *
      * 1. multiple DirectInfusionMatchData match to exactly the same fragments,
      * 2. the DirectInfusionMatchData's associated compounds can be readily summarized to a common form,
      * and
      * 3. the DirectInfusionMatchData's matched adducts are of the same type,
      *
      * these DirectInfusionMatchData can be summarized.
      *
      * When these conditions are met, the summary mappings are saved in these maps
      */
    map<shared_ptr<DirectInfusionMatchData>, string> originalMatchToSummaryString = {};
    map<string, set<shared_ptr<DirectInfusionMatchData>>> chainLengthSummaries = {}; //LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()
    map<string, set<shared_ptr<DirectInfusionMatchData>>> compositionSummaries = {}; //LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()

    //summarized
    map<int, vector<shared_ptr<DirectInfusionMatchData>>> fragToMatchDataSummarized = {};
    map<shared_ptr<DirectInfusionMatchData>, vector<int>> matchDataToFragsSummarized  = {};
    map<pair<int, shared_ptr<DirectInfusionMatchData>>,float> fragToTheoreticalIntensitySummarized = {};


    float getNormalizedTheoreticalIntensity(int fragId, shared_ptr<DirectInfusionMatchData> matchData){
        pair<int, shared_ptr<DirectInfusionMatchData>> pair = make_pair(fragId, matchData);
        if (fragToTheoreticalIntensitySummarized.find(pair) != fragToTheoreticalIntensitySummarized.end()) {
            return fragToTheoreticalIntensitySummarized.at(pair);
        } else if (fragToTheoreticalIntensity.find(pair) != fragToTheoreticalIntensity.end()) {
            return fragToTheoreticalIntensity.at(pair);
        } else {
            return 0.0f; //TODO: should this throw an error?
        }
    }

    float getObservedIntensity(int fragId, shared_ptr<DirectInfusionMatchData> matchData){
        return fragToObservedIntensity.at(make_pair(fragId, matchData));
    }

    float getIntensityRatio(int fragId, shared_ptr<DirectInfusionMatchData> matchData){
        return (getObservedIntensity(fragId, matchData) / getNormalizedTheoreticalIntensity(fragId, matchData));
    }

    shared_ptr<DirectInfusionSinglePeakMatchData> getSinglePeakMatchData(int fragId, shared_ptr<DirectInfusionMatchData> matchData){

        shared_ptr<DirectInfusionSinglePeakMatchData> directInfusionSinglePeakMatchData = shared_ptr<DirectInfusionSinglePeakMatchData>(new DirectInfusionSinglePeakMatchData());

        directInfusionSinglePeakMatchData->normalizedTheoreticalIntensity = getNormalizedTheoreticalIntensity(fragId, matchData);
        directInfusionSinglePeakMatchData->observedIntensity = getObservedIntensity(fragId, matchData);

        return directInfusionSinglePeakMatchData;
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
             shared_ptr<DirectInfusionSearchParameters> params,
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
      * fragments converted m/z <--> int keys using mzToIntKey(mz, 1000000) and intKeyToMz(intKey, 1000000).
      *
      * Note that this function does no processing, filtering, or analysis - it simply reorganizes
      * the compound match data into maps.
      */
     static unique_ptr<DirectInfusionMatchInformation> getMatchInformation(
             vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
             Fragment *observedSpectrum,
             shared_ptr<DirectInfusionSearchParameters> params,
             bool debug);


     /**
      * Divide up search into smaller, more manageable parts
      * mapKey --> compounds in the database
      *        --> relevant parts of each individual sample to search
      *
      * --> a function to divide the raw file into relevant pieces
      *     relevant pieces are MS2 scans of interest
      *     ie, scansByPrecursorMap
      */
     //static DirectInfusionAnnotation*
     static DirectInfusionAnnotation* processBlock(int blockNum,
                              const pair<float,float>& mzRange,
                              mzSample* sample,
                              const vector<Scan*>& ms2Scans,
                              const vector<Scan*>& ms1Scans,
                              const vector<pair<Compound*, Adduct*>> library,
                              const shared_ptr<DirectInfusionSearchParameters> params,
                              const bool debug);

     /**
      * @brief assessMatch
      * @param ms2Scans
      * @param ms1Scans
      * @param libraryMatch
      * @param params
      * @param debug
      *
      * Designed to be multithreaded, work of comparing / evaluating individual matches
      */
     static pair<FragmentationMatchScore, float> assessMatch(const Fragment *f,
                             const vector<Scan*>& ms1Scans,
                             const pair<Compound*, Adduct*>& libraryMatch,
                             const shared_ptr<DirectInfusionSearchParameters> params,
                             const bool debug);
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
class DirectInfusionGroupAnnotation : public DirectInfusionAnnotation {

public:

    /**
     * @brief annotationBySample
     * retain original samples for reference.
     */
    map<mzSample*, DirectInfusionAnnotation*> annotationBySample = {};

    void clean();

    /**
     * @brief createByAverageProportions
     * @param singleSampleAnnotations
     *
     * @return pointer to DirectInfusionGroupAnnotation object.
     * This pointer must be deleted explicitly! Cannot use smart pointers b/c of QMetaType rules.
     */
    static DirectInfusionGroupAnnotation* createByAverageProportions(
            vector<DirectInfusionAnnotation*> singleSampleAnnotations,
            shared_ptr<DirectInfusionSearchParameters> params,
            bool debug);

};

typedef map<int, vector<shared_ptr<DirectInfusionMatchData>>>::iterator fragToMatchDataIterator;
typedef map<shared_ptr<DirectInfusionMatchData>, vector<int>>::iterator matchDataToFragIterator;
typedef map<shared_ptr<DirectInfusionMatchData>, vector<shared_ptr<DirectInfusionSinglePeakMatchData>>>::iterator matchDataToFragIntensityIterator;
typedef map<shared_ptr<DirectInfusionMatchData>, float>::iterator matchDataToFloatIterator;
typedef map<string, std::set<shared_ptr<DirectInfusionMatchData>>>::iterator stringToMatchDataIterator;
typedef map<shared_ptr<DirectInfusionMatchData>, string>::iterator matchToStringIterator;

