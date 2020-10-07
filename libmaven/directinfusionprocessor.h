#pragma once

#include "mzSample.h"
#include "mzUtils.h"
#include "Fragment.h"
#include <memory>
#include <algorithm>
#include <sstream>
#include <numeric>

class mzSample;
class DirectInfusionAnnotation;
class Ms3SingleSampleMatch;
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
    ALL_CANDIDATES,                                     // no summarization, no quant
    AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE,   // uses lipid summarizations + quant approach
    AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION,        // lipid summarization without quant
    AUTO_SUMMARIZED_IDENTICAL_FRAGMENTS                 // summarize compounds with identical matches together, regardless of structural relationships
};

enum class Ms3IntensityType{
    ALL_MATCHES,                //sum of all intensity of all Ms3 matches in tolerance
    MAX_INTENSITY,              // highest intensity peak in tolerance
    CLOSEST_MZ                  // closest m/z to theoretical match
};

/**
 * @brief The DirectInfusionSearchParameters class
 *
 * single class to contain all parameters used in direct infusion search analysis.
 */
class DirectInfusionSearchParameters : public SearchParameters {

public:

    /** ===================
     * MS1 SEARCH RELATED
     * @param ms1IsRequireAdductPrecursorMatch:reference compound associated adduct must == query adduct
     * @param ms1IsFindPrecursorIon: only retain matches where precursor peak is found in MS1 scan(s).
     * @param ms1PpmTolr: tolerance value used for matching theoretical ion m/z to an m/z peak in an MS1 scan
     * @param ms1MinIntensity: min intensity value for a MS1 spectral peak to be considered real
     * @param ms1ScanFilter: consider only MS1 scans that substring match in their filterString field to this
     * @param isRequireMonoisotopic: disqualify candidate MS1 peak if there is a peak exactly one 13C-12C width behind
     * @param mMinusOnePeakMaxIntensityFraction: Avoid disqualifying candidate MS1 peak
     *          if candidate monoisotopic peak intensity / candidate MS1 peak intensity <= this fraction
     * ==================== */

    bool ms1IsRequireAdductPrecursorMatch = true;
    bool ms1IsFindPrecursorIon = false;
    float ms1PpmTolr = 5;
    float ms1MinIntensity = 0;
    string ms1ScanFilter = "";
    bool ms1IsRequireMonoisotopic = true;
    float ms1MMinusOnePeakMaxIntensityFraction = 1.0f;

    /** ===================
     * MS3 SEARCH RELATED
     * @param ms3IsMs3Search: if experimental data contains MS3 scans
     * @param ms3MinNumMatches: Minimum number of reference MS3 peaks found in observed MS3 scans
     * @param ms3MinNumMs3MzMatches: Minimum number of reference MS3 m/zs found in observed MS3 scans, regardless of precursor m/zs
     * @param ms3AnalysisMs1PrecursorPpmTolr: m/z tolerance used for matching reference MS1 <--> MS3 scan precursor m/z
     * @param ms3PrecursorPpmTolr: m/z tolerance value used for matching reference MS2 m/z <--> MS3 scan precursor m/z
     * @param ms3PpmTolr: m/z tolerance value used for matching reference <--> observed spectral peaks in MS3 spectrum
     * ==================== */
    bool ms3IsMs3Search = false;
    int ms3MinNumMatches = 1;
    int ms3MinNumMs3MzMatches = 1;
    float ms3AnalysisMs1PrecursorPpmTolr = 20;
    float ms3PrecursorPpmTolr = 20;
    float ms3MatchTolrInDa = 0.5f;
    float ms3MinIntensity = 0;
    Ms3IntensityType ms3IntensityType = Ms3IntensityType::MAX_INTENSITY;

    /** =======================
     * MS3 SCAN ASSOCIATED
     * @param scanFilterMs3MinRt: min RT for valid MS3 scan (otherwise excluded). -1 to ignore.
     * @param scanFilterMs3MaxRt: max RT for valid MS3 scan (otherwise excluded). -1 to ignore.
     * @param isPreferSmallestMassWindow: used in DirectInfusionUtils::findNormalizedIntensity().
     *              Agglomerate results from all individual scans based on mass window lengths.
     * ========================*/
    float scanFilterMs3MinRt = -1.0f;
    float scanFilterMs3MaxRt = -1.0f;

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

    /** ===================
     * INTENSITY COMPUTATION
     * @param ms1PartitionIntensityByFragments: list of fragments to sum to determine a fractional value
     *      for splitting ms1 intensity, when multiple compounds map to the same ms1 intensity value.
     * @param isPreferSmallestMassWindow: used in DirectInfusionUtils::findNormalizedIntensity().
     *       Agglomerate results from all individual scans based on mass window lengths.
     * ==================== */
    vector<string> ms1PartitionIntensityByFragments{"sn1","sn2"};
    bool isPreferSmallestScanMassWindow = true;

    //Issue 270
    bool isReduceBySimpleParsimony = false;

    bool isDiagnosticFragmentMapAgreement(map<string, int> observedNumDiagnosticMatchesMap){

        for (auto it = ms2MinNumDiagnosticMatchesMap.begin(); it != ms2MinNumDiagnosticMatchesMap.end(); ++it) {
            string key = it->first;

            if (observedNumDiagnosticMatchesMap.find(key) != observedNumDiagnosticMatchesMap.end()){
                int numObservedDiagnosticMatches = observedNumDiagnosticMatchesMap[key];
                if (numObservedDiagnosticMatches < it->second){
                    return false;   //insufficient count
                }
            } else {
                if (it->second > 0) {
                    return false; //not finding a key implies count of 0
                }
            }
        }

        return true;
    }

    void printParams(){
        string encodedParams = encodeParams();
        replace(encodedParams.begin(), encodedParams.end(), ';', '\n');
        replace(encodedParams.begin(), encodedParams.end(), '=', ' ');
        cout << encodedParams << endl;
    }

    //RESERVED DELIMITERS - DO NOT CHANGE!
    static constexpr const char* const INTERNAL_MAP_DELIMITER = "|,|";

    string encodeParams() {

        string encodedParams;

        //scan filter params (all ms levels)
        encodedParams = encodedParams + "scanFilterMinFracIntensity" + "=" + to_string(scanFilterMinFracIntensity) + ";";
        encodedParams = encodedParams + "scanFilterMinSNRatio" + "=" + to_string(scanFilterMinSNRatio) + ";";
        encodedParams = encodedParams + "scanFilterMaxNumberOfFragments" + "=" + to_string(scanFilterMaxNumberOfFragments) + ";";
        encodedParams = encodedParams + "scanFilterBaseLinePercentile" + "=" + to_string(scanFilterBaseLinePercentile) + ";";
        encodedParams = encodedParams + "scanFilterIsRetainFragmentsAbovePrecursorMz" + "=" + to_string(scanFilterIsRetainFragmentsAbovePrecursorMz) + ";";
        encodedParams = encodedParams + "scanFilterPrecursorPurityPpm" + "=" + to_string(scanFilterPrecursorPurityPpm) + ";";
        encodedParams = encodedParams + "scanFilterMinIntensity" + "=" + to_string(scanFilterMinIntensity) + ";";

        //scan filter for MS1 scans
        encodedParams = encodedParams + "scanFilterMs1MinRt" + "=" + to_string(scanFilterMs1MinRt) + ";";
        encodedParams = encodedParams + "scanFilterMs1MaxRt" + "=" + to_string(scanFilterMs1MaxRt) + ";";

        //scan filter for MS2 scans
        encodedParams = encodedParams + "scanFilterMs2MinRt" + "=" + to_string(scanFilterMs2MinRt) + ";";
        encodedParams = encodedParams + "scanFilterMs2MaxRt" + "=" + to_string(scanFilterMs2MaxRt) + ";";

        //scan filter for MS3 scans
        encodedParams = encodedParams + "scanFilterMs3MinRt" + "=" + to_string(scanFilterMs3MinRt) + ";";
        encodedParams = encodedParams + "scanFilterMs3MaxRt" + "=" + to_string(scanFilterMs3MaxRt) + ";";

        //consensus spectrum params (all ms levels)
        encodedParams = encodedParams + "consensusIsIntensityAvgByObserved" + "=" + to_string(consensusIsIntensityAvgByObserved) + ";";
        encodedParams = encodedParams + "consensusIsNormalizeTo10K" + "=" + to_string(consensusIsNormalizeTo10K) + ";";
        string consensusIntensityAgglomerationTypeStr = "UNSPECIFIED";
        if (consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
            consensusIntensityAgglomerationTypeStr = "MEAN";
        } else if (consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
            consensusIntensityAgglomerationTypeStr = "MEDIAN";
        }
        encodedParams = encodedParams + "consensusIntensityAgglomerationType" + "=" + consensusIntensityAgglomerationTypeStr + ";";

        //ms1 consensus spectrum params
        encodedParams = encodedParams + "consensusMs1PpmTolr" + "=" + to_string(consensusMs1PpmTolr) + ";";
        encodedParams = encodedParams + "consensusMinNumMs1Scans" + "=" + to_string(consensusMinNumMs1Scans) + ";";
        encodedParams = encodedParams + "consensusMinFractionMs1Scans" + "=" + to_string(consensusMinFractionMs1Scans) + ";";

        //ms2 consensus spectrum params
        encodedParams = encodedParams + "consensusPpmTolr" + "=" + to_string(consensusPpmTolr) + ";";
        encodedParams = encodedParams + "consensusMinNumMs2Scans" + "=" + to_string(consensusMinNumMs2Scans) + ";";
        encodedParams = encodedParams + "consensusMinFractionMs2Scans" + "=" + to_string(consensusMinFractionMs2Scans) + ";";

        //ms3 search params
        encodedParams = encodedParams + "ms3IsMs3Search" + "=" + to_string(ms3IsMs3Search) + ";";
        encodedParams = encodedParams + "ms3MinNumMatches" + "=" + to_string(ms3MinNumMatches) + ";";
        encodedParams = encodedParams + "ms3MinNumMs3MzMatches" + "=" + to_string(ms3MinNumMs3MzMatches) + ";";
        encodedParams = encodedParams + "ms3AnalysisMs1PrecursorPpmTolr" + "=" + to_string(ms3AnalysisMs1PrecursorPpmTolr) + ";";
        encodedParams = encodedParams + "ms3PrecursorPpmTolr" + "=" + to_string(ms3PrecursorPpmTolr) + ";";
        encodedParams = encodedParams + "ms3MatchTolrInDa" + "=" + to_string(ms3MatchTolrInDa) + ";";
        encodedParams = encodedParams + "ms3MinIntensity" + "=" + to_string(ms3MinIntensity) + ";";

        string ms3IntensityTypeStr = "UNSPECIFIED";
        if (ms3IntensityType == Ms3IntensityType::CLOSEST_MZ) {
            ms3IntensityTypeStr = "CLOSEST_MZ";
        } else if (ms3IntensityType == Ms3IntensityType::MAX_INTENSITY) {
            ms3IntensityTypeStr = "MAX_INTENSITY";
        } else if (ms3IntensityType == Ms3IntensityType::ALL_MATCHES) {
            ms3IntensityTypeStr = "ALL_MATCHES";
        }
        encodedParams = encodedParams + "ms3IntensityType" + "=" + ms3IntensityTypeStr + ";";

        //ms2 search params
        encodedParams = encodedParams + "ms2MinNumMatches" + "=" + to_string(ms2MinNumMatches) + ";";
        encodedParams = encodedParams + "ms2MinNumDiagnosticMatches" + "=" + to_string(ms2MinNumDiagnosticMatches) + ";";
        encodedParams = encodedParams + "ms2MinNumUniqueMatches" + "=" + to_string(ms2MinNumUniqueMatches) + ";";
        encodedParams = encodedParams + "ms2PpmTolr" + "=" + to_string(ms2PpmTolr) + ";";
        encodedParams = encodedParams + "ms2MinIntensity" + "=" + to_string(ms2MinIntensity) + ";";

        encodedParams = encodedParams + "ms2MinNumDiagnosticMatchesMap" + "=" + "{";

        for (auto it = ms2MinNumDiagnosticMatchesMap.begin(); it != ms2MinNumDiagnosticMatchesMap.end(); ++it) {
            string key = it->first;
            string value = to_string(it->second);
            encodedParams = encodedParams + key + "=" + value + INTERNAL_MAP_DELIMITER;
        }

        encodedParams = encodedParams + "};";

        //ms1 search params
        encodedParams = encodedParams + "ms1IsRequireAdductPrecursorMatch" + "=" + to_string(ms1IsRequireAdductPrecursorMatch) + ";";
        encodedParams = encodedParams + "ms1IsFindPrecursorIon" + "=" + to_string(ms1IsFindPrecursorIon) + ";";
        encodedParams = encodedParams + "ms1PpmTolr" + "=" + to_string(ms1PpmTolr) + ";";
        encodedParams = encodedParams + "ms1MinIntensity" + "=" + to_string(ms1MinIntensity) + ";";
        encodedParams = encodedParams + "ms1ScanFilter" + "=" + ms1ScanFilter + ";";
        encodedParams = encodedParams + "ms1IsRequireMonoisotopic" + "=" + to_string(ms1IsRequireMonoisotopic) + ";";
        encodedParams = encodedParams + "ms1MMinusOnePeakMaxIntensityFraction" + "=" + to_string(ms1MMinusOnePeakMaxIntensityFraction) + ";";

        //DIMS intensity options
        encodedParams = encodedParams + "ms1PartitionIntensityByFragments" + "=" + "{";
        for (auto fragmentLabel : ms1PartitionIntensityByFragments) {
            encodedParams = encodedParams + fragmentLabel + INTERNAL_MAP_DELIMITER;
        }
        encodedParams = encodedParams + "};";
        encodedParams = encodedParams + "isPreferSmallestScanMassWindow" + "=" + to_string(isPreferSmallestScanMassWindow) + ";";

        //agglomeration params
        encodedParams = encodedParams + "isAgglomerateAcrossSamples" + "=" + to_string(isAgglomerateAcrossSamples) + ";";

        string spectralCompositionAlgorithmStr = "UNSPECIFIED";
        if (spectralCompositionAlgorithm == SpectralCompositionAlgorithm::ALL_CANDIDATES) {
            spectralCompositionAlgorithmStr = "ALL_CANDIDATES";
        } else if (spectralCompositionAlgorithm == SpectralCompositionAlgorithm::AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE) {
            spectralCompositionAlgorithmStr = "AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE";
        } else if (spectralCompositionAlgorithm == SpectralCompositionAlgorithm::AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION) {
            spectralCompositionAlgorithmStr = "AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION";
        } else if (spectralCompositionAlgorithm == SpectralCompositionAlgorithm::AUTO_SUMMARIZED_IDENTICAL_FRAGMENTS) {
            spectralCompositionAlgorithmStr = "AUTO_SUMMARIZED_IDENTICAL_FRAGMENTS";
        }
        encodedParams = encodedParams + "spectralCompositionAlgorithm" + "=" + spectralCompositionAlgorithmStr + ";";

        encodedParams = encodedParams + "isReduceBySimpleParsimony" + "=" + to_string(isReduceBySimpleParsimony) + ";";

        return encodedParams;
    }

    static void addMs2MinNumDiagnosticMatchesMap(shared_ptr<DirectInfusionSearchParameters> directInfusionSearchParameters, string encodedMs2MinNumDiagnosticMatchesMap){
        unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedMs2MinNumDiagnosticMatchesMap, INTERNAL_MAP_DELIMITER);
        for (auto it = decodedMap.begin(); it != decodedMap.end(); ++it){
            string key = it->first;
            int value = stoi(it->second);
            directInfusionSearchParameters->ms2MinNumDiagnosticMatchesMap.insert(make_pair(key, value));
        }
    }

    static shared_ptr<DirectInfusionSearchParameters> decode(string encodedParams){
        shared_ptr<DirectInfusionSearchParameters> directInfusionSearchParameters = shared_ptr<DirectInfusionSearchParameters>(new DirectInfusionSearchParameters());

        unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

        //scan filter params
        if (decodedMap.find("scanFilterMinFracIntensity") != decodedMap.end()){
            directInfusionSearchParameters->scanFilterMinFracIntensity = stof(decodedMap["scanFilterMinFracIntensity"]);
        }
        if (decodedMap.find("scanFilterMinSNRatio") != decodedMap.end()){
            directInfusionSearchParameters->scanFilterMinSNRatio = stof(decodedMap["scanFilterMinSNRatio"]);
        }
        if (decodedMap.find("scanFilterMaxNumberOfFragments") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterMaxNumberOfFragments = stoi(decodedMap["scanFilterMaxNumberOfFragments"]);
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

        //scan filter for MS1 scans
        if (decodedMap.find("scanFilterMs1MinRt") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterMs1MinRt = stof(decodedMap["scanFilterMs1MinRt"]);
        }
        if (decodedMap.find("scanFilterMs1MaxRt") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterMs1MaxRt = stof(decodedMap["scanFilterMs1MaxRt"]);
        }

        //scan filter for MS2 scans
        if (decodedMap.find("scanFilterMs2MinRt") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterMs2MinRt = stof(decodedMap["scanFilterMs2MinRt"]);
        }
        if (decodedMap.find("scanFilterMs2MaxRt") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterMs2MaxRt = stof(decodedMap["scanFilterMs2MaxRt"]);
        }

        //scan filter for MS3 scan
        if (decodedMap.find("scanFilterMs3MinRt") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterMs3MinRt = stof(decodedMap["scanFilterMs3MinRt"]);
        }
        if (decodedMap.find("scanFilterMs3MaxRt") != decodedMap.end()) {
            directInfusionSearchParameters->scanFilterMs3MaxRt = stof(decodedMap["scanFilterMs3MaxRt"]);
        }

        //consensus spectrum params (all ms levels)

        if (decodedMap.find("consensusIsIntensityAvgByObserved") != decodedMap.end()){
            directInfusionSearchParameters->consensusIsIntensityAvgByObserved = decodedMap["consensusIsIntensityAvgByObserved"] == "1";
        }
        if (decodedMap.find("consensusIsNormalizeTo10K") != decodedMap.end()){
            directInfusionSearchParameters->consensusIsNormalizeTo10K = decodedMap["consensusIsNormalizeTo10K"] == "1";
        }
        if (decodedMap.find("consensusIntensityAgglomerationType") != decodedMap.end()) {
            string consensusIntensityAgglomerationTypeStr = decodedMap["consensusIntensityAgglomerationType"];
            if (consensusIntensityAgglomerationTypeStr == "MEAN") {
                directInfusionSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
            } else if (consensusIntensityAgglomerationTypeStr == "MEDIAN") {
                directInfusionSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
            }
        }

        //ms1 consensus spectrum params
        if (decodedMap.find("consensusMs1PpmTolr") != decodedMap.end()){
            directInfusionSearchParameters->consensusMs1PpmTolr = stof(decodedMap["consensusMs1PpmTolr"]);
        }
        if (decodedMap.find("consensusMinNumMs1Scans") != decodedMap.end()){
            directInfusionSearchParameters->consensusMinNumMs1Scans = stoi(decodedMap["consensusMinNumMs1Scans"]);
        }
        if (decodedMap.find("consensusMinFractionMs1Scans") != decodedMap.end()){
            directInfusionSearchParameters->consensusMinFractionMs1Scans = stof(decodedMap["consensusMinFractionMs1Scans"]);
        }

        //ms2 consensus spectrum params
        if (decodedMap.find("consensusPpmTolr") != decodedMap.end()){
            directInfusionSearchParameters->consensusPpmTolr = stof(decodedMap["consensusPpmTolr"]);
        }
        if (decodedMap.find("consensusMinNumMs2Scans") != decodedMap.end()){
            directInfusionSearchParameters->consensusMinNumMs2Scans = stoi(decodedMap["consensusMinNumMs2Scans"]);
        }
        if (decodedMap.find("consensusMinFractionMs2Scans") != decodedMap.end()){
            directInfusionSearchParameters->consensusMinFractionMs2Scans = stof(decodedMap["consensusMinFractionMs2Scans"]);
        }

        //ms3 search params
        if (decodedMap.find("ms3IsMs3Search") != decodedMap.end()) {
            directInfusionSearchParameters->ms3IsMs3Search = decodedMap["ms3IsMs3Search"] == "1";
        }
        if (decodedMap.find("ms3MinNumMatches") != decodedMap.end()) {
            directInfusionSearchParameters->ms3MinNumMatches = stoi(decodedMap["ms3MinNumMatches"]);
        }
        if (decodedMap.find("ms3MinNumMs3MzMatches") != decodedMap.end()) {
            directInfusionSearchParameters->ms3MinNumMs3MzMatches = stoi(decodedMap["ms3MinNumMs3MzMatches"]);
        }
        if (decodedMap.find("ms3AnalysisMs1PrecursorPpmTolr") != decodedMap.end()) {
            directInfusionSearchParameters->ms3AnalysisMs1PrecursorPpmTolr = stof(decodedMap["ms3AnalysisMs1PrecursorPpmTolr"]);
        }
        if (decodedMap.find("ms3PrecursorPpmTolr") != decodedMap.end()) {
            directInfusionSearchParameters->ms3PrecursorPpmTolr = stof(decodedMap["ms3PrecursorPpmTolr"]);
        }
        if (decodedMap.find("ms3MatchTolrInDa") != decodedMap.end()) {
            directInfusionSearchParameters->ms3MatchTolrInDa = stof(decodedMap["ms3MatchTolrInDa"]);
        }
        if (decodedMap.find("ms3MinIntensity") != decodedMap.end()) {
            directInfusionSearchParameters->ms3MinIntensity = stof(decodedMap["ms3MinIntensity"]);
        }
        if (decodedMap.find("ms3IntensityType") != decodedMap.end()) {
            string ms3IntensityTypeStr = decodedMap["ms3IntensityType"];
            if (ms3IntensityTypeStr == "CLOSEST_MZ") {
                directInfusionSearchParameters->ms3IntensityType = Ms3IntensityType::CLOSEST_MZ;
            } else if (ms3IntensityTypeStr == "MAX_INTENSITY") {
                directInfusionSearchParameters->ms3IntensityType = Ms3IntensityType::MAX_INTENSITY;
            } else if (ms3IntensityTypeStr == "ALL_MATCHES") {
                directInfusionSearchParameters->ms3IntensityType = Ms3IntensityType::ALL_MATCHES;
            }
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
        if (decodedMap.find("ms2MinNumDiagnosticMatchesMap") != decodedMap.end()) {
            string encodedDiagnosticFragmentsMap = decodedMap["ms2MinNumDiagnosticMatchesMap"];
            addMs2MinNumDiagnosticMatchesMap(directInfusionSearchParameters, encodedDiagnosticFragmentsMap);
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
        if (decodedMap.find("ms1IsRequireMonoisotopic") != decodedMap.end()) {
            directInfusionSearchParameters->ms1IsRequireMonoisotopic = decodedMap["ms1IsRequireMonoisotopic"] == "1";
        }
        if (decodedMap.find("ms1MMinusOnePeakMaxIntensityFraction") != decodedMap.end()) {
            directInfusionSearchParameters->ms1MMinusOnePeakMaxIntensityFraction = stof(decodedMap["ms1MMinusOnePeakMaxIntensityFraction"]);
        }

        //DIMS intensity options
        if (decodedMap.find("ms1PartitionIntensityByFragments") != decodedMap.end()){
            string encodedMs1PartitionIntensityByFragments = decodedMap["ms1PartitionIntensityByFragments"];
            directInfusionSearchParameters->ms1PartitionIntensityByFragments = mzUtils::decodeParameterVector(encodedMs1PartitionIntensityByFragments, INTERNAL_MAP_DELIMITER);
        }
        if (decodedMap.find("isPreferSmallestScanMassWindow") != decodedMap.end()) {
            directInfusionSearchParameters->isPreferSmallestScanMassWindow = decodedMap["isPreferSmallestScanMassWindow"] == "1";
        }

        //agglomeration params
        if (decodedMap.find("isAgglomerateAcrossSamples") != decodedMap.end()){
            directInfusionSearchParameters->isAgglomerateAcrossSamples = decodedMap["isAgglomerateAcrossSamples"] == "1";
        }
        if (decodedMap.find("spectralCompositionAlgorithm") != decodedMap.end()){
            string spectralCompositionAlgorithmStr = decodedMap["spectralCompositionAlgorithm"];
            if (spectralCompositionAlgorithmStr == "ALL_CANDIDATES") {
                directInfusionSearchParameters->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::ALL_CANDIDATES;
            } else if (spectralCompositionAlgorithmStr == "AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE") {
                directInfusionSearchParameters->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE;
            } else if (spectralCompositionAlgorithmStr == "AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION") {
                directInfusionSearchParameters->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION;
            } else if (spectralCompositionAlgorithmStr == "AUTO_SUMMARIZED_IDENTICAL_FRAGMENTS") {
                directInfusionSearchParameters->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::AUTO_SUMMARIZED_IDENTICAL_FRAGMENTS;
            }
        }
        if (decodedMap.find("isReduceBySimpleParsimony") != decodedMap.end()) {
            directInfusionSearchParameters->isReduceBySimpleParsimony = decodedMap["isReduceBySimpleParsimony"] == "1";
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

    //Issue 232
    float observedMs1Intensity = 0;

    //Issue 210
    int numUniqueFragments = 0;
    vector<bool> isFragmentUnique; //follows m/z-sorted Compound* fragment vectors

    //Issue 288
    int ms1IntensityCoord = -1;
    float ms1PartitionFraction = 1.0f;

    //Issue 292
    float ms1PartitionFractionByScan = 1.0f;
};

/**
 * @brief The DirectInfusionMatchAssessment struct
 *
 * as returned by DirectInfusionProcessor::assessMatch()
 */
struct DirectInfusionMatchAssessment {
    FragmentationMatchScore fragmentationMatchScore;
    map<string, int> diagnosticFragmentMatchMap = {};
    float fragmentMaxObservedIntensity = 0;
    float observedMs1Intensity = 0;
    int ms1IntensityCoord = -1;
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
 * @brief The DirectInfusionMatchInformation structure
 *
 * A structure to organize all fragment matches from all (compound, adduct) pairs that match to a single
 * direct infusion ms2 spectrum.
 *
 * Provides maps of individual fragments to compound matches, and compound to fragment matches.
 *
 */
struct DirectInfusionMatchInformation {

public:

    //group of identified fragment m/zs ==> all identified compounds
    map<vector<int>, vector<shared_ptr<DirectInfusionMatchData>>> fragListToCompounds = {};

    //single compound ==> all identified fragment mz/s
    map<shared_ptr<DirectInfusionMatchData>, vector<int>> matchDataToFrags = {};

    //single fragment m/z ==> all identified compounds containing fragment
    map<int, unordered_set<shared_ptr<DirectInfusionMatchData>>> fragToMatchData = {};

    /**
     * @brief getCompounds
     *
     * method to compute current set of compounds, based on the current status of the
     * data within this struct
     *
     * This is intentionally left as a method instead of a field to avoid duplicated state
     * dependencies.
     *
     * @return
     */
    vector<shared_ptr<DirectInfusionMatchData>> getCompounds();

    //Issue 275
    string getFragmentGroupId(shared_ptr<DirectInfusionMatchData> compound, int precision=2);

    //Issue 288, 292
    void computeMs1PartitionFractions(const vector<Scan*>& ms2Scans,
                                      const Fragment *ms1Fragment,
                                      const shared_ptr<DirectInfusionSearchParameters> params,
                                      const bool debug);

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
      * @brief getMs3CompoundSet
      * @param compounds
      * @param debug
      * @return
      *
      * vector of Ms3Compound* objects, where each Ms3Compound* derives from an input Compound*.
      */
     static vector<Ms3Compound*> getMs3CompoundSet(const vector<Compound*>& compounds,
                                                   bool debug);

    /**
     * @brief processSingleSample
     * @param sample
     * @param directInfusionSearchSet
     * @param debug
     * @return
     *
     * Returns DirectInfusionAnnotation assessments for a single sample.
     */
     static map<int, DirectInfusionAnnotation*> processSingleSample(
             mzSample *sample,
             shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
             shared_ptr<DirectInfusionSearchParameters> params,
             bool debug);

     /**
      * @brief processSingleMs3Sample
      * @param sample
      * @param ms3Compounds
      * @param params
      * @param debug
      * @return
      *
      * DirectInfusionAnnotation assessments for a single sample.
      *
      * TODO: MS3 spectra are organized into groups to build consensus spectra,
      * based on MS2 precursor m/z.
      */
     static vector<Ms3SingleSampleMatch*> processSingleMs3Sample(
             mzSample* sample,
             const vector<Ms3Compound*>& ms3Compounds,
             shared_ptr<DirectInfusionSearchParameters> params,
             bool debug);

     /**
      * @brief organizeMs3ScansByPrecursor
      * @param allMs3Scans
      * @param debug
      * @return
      *
      * Group MS3 Scans together by MS2 precursors, according to ms3 precursor ppm tolerance parameter.
      *
      * This parameter is often stored in DirectInfusionSearchParameters as
      *
      * params->ms3PrecursorPpmTolr
      *
      */
     static vector<vector<tuple<double, double, Scan*>>> organizeMs3ScansByPrecursor(
             vector<tuple<double, double, Scan*>> allMs3Scans,
             double ms2PrecursorPpmTolr,
             double ms3PrecursorPpmTolr,
             bool debug);

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
      * @brief DirectInfusionProcessor::getFragmentMatchMaps
      * @param allCandidates
      * @param observedSpectrum
      * @param params
      * @param debug
      * @return
      *
      * Determine maps of fragments to compounds based on current match information.
      */
     static unique_ptr<DirectInfusionMatchInformation> getFragmentMatchMaps(
             vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
             Fragment *observedSpectrum,
             shared_ptr<DirectInfusionSearchParameters> params,
             bool debug);

     /**
      * @brief reduceBySimpleParsimony
      * @param matchInfo
      * @param params
      * @param debug
      * @return
      * Issue 270
      *
      * Adjusts data structure and reduces based on simple parsimony
      *
      * Remove a fragment group if the fragments in that group are all
      * also found in fragment groups that contain more fragments.
      *
      * Also need to make corresponding adjustments to other data structures
      * in DirectInfusionMatchInformation.
      */
     static unique_ptr<DirectInfusionMatchInformation> reduceBySimpleParsimony(
             unique_ptr<DirectInfusionMatchInformation> matchInfo,
             shared_ptr<DirectInfusionSearchParameters> params,
             bool debug);

     /**
      * @brief DirectInfusionProcessor::summarizeFragmentGroups
      * @param matchInfo
      * @param observedSpectrum
      * @param params
      * @param debug
      * @return
      */
     static unique_ptr<DirectInfusionMatchInformation> summarizeFragmentGroups(
             unique_ptr<DirectInfusionMatchInformation> matchInfo,
             Fragment *observedSpectrum,
             shared_ptr<DirectInfusionSearchParameters> params,
             bool debug);


    static unique_ptr<DirectInfusionMatchInformation> reduceByUniqueMatches(
            unique_ptr<DirectInfusionMatchInformation> matchInfo,
            Fragment* observedSpectrum,
            shared_ptr<DirectInfusionSearchParameters> params,
            bool debug);


     static void addBlockSpecificMatchInfo(
             DirectInfusionMatchInformation *matchInfo,
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
     static DirectInfusionAnnotation* processBlock(int blockNum,
                              const pair<float,float>& mzRange,
                              mzSample* sample,
                              const vector<Scan*>& ms2Scans,
                              const Fragment *ms1Fragment, //only one per sample, computed at the same time that ms1 scans are retrieved.
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
     static unique_ptr<DirectInfusionMatchAssessment> assessMatch(const Fragment *f, //ms2 fragment
                             const Fragment *ms1Fragment,
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
     * @brief matchInformation
     * series mappings between compounds and identified ms2 fragments
     */
    unique_ptr<DirectInfusionMatchInformation> matchInformation;

    /**
     * compound, adduct, and estimated proportion of the spectrum
     * associated with the match.
     *
     * FragmentationMatchScores are also provided.
     */
    vector<shared_ptr<DirectInfusionMatchData>> compounds;
};

class DirectInfusionUtils {
public:
    static float findNormalizedIntensity(const vector<Scan*>& scans,
                                         float queryMz,
                                         float standardMz,
                                         shared_ptr<DirectInfusionSearchParameters> params,
                                         bool debug = false);

    static constexpr double C_13_MASS = 1.00335483521;

};

//Issue 226
class Ms3Annotation {
public:
    map<mzSample*,Ms3SingleSampleMatch*> matchesBySample{};
};

/**
 * @brief The Ms3SingleSampleMatch class
 *
 * Unlike @class{DirectInfusionAnnotation}, which is data-centric, this class is compound-centric.
 *
 * This class contains an Ms3 compound, a sample, and a path to all matching information
 * associated with the compound.
 */
class Ms3SingleSampleMatch {

public:

    Ms3Compound *ms3Compound = nullptr;

    mzSample *sample = nullptr;

    int numMs3Matches = 0;
    int numMs3MzMatches = 0;
    float observedMs1Intensity = 0;

    //MS3 VALUES ARE MAPPED BASED ON POSITIONS IN Ms3Compound

    //Issue 295: Keep track of al fragment intensities from all scans, preserving precursors
    //
    //<ms2 m/z, vector position in Ms3Compound maps>
    map<pair<int, int>, vector<float>> scanIntensitiesByMs1Ms2Ms3Mzs{};

    //Issue 295: summarize each ms3 m/z values into a single intensity
    //
    // <ms2 m/z, vector position in Ms3Compound maps>
    map<pair<int, int>, float> intensityByMs1Ms2Ms3Mzs{};

    // ---------------------------------------------------- //

    // MS3 VALUES ARE MAPPED BASED ON MS3 FRAGMENT m/z KEY

    //Issue 295: Keep track of all fragment intensities from all scans,
    //Even if these ms3 fragments are detected from different precursors.
    //
    //ms3MzKey
    map<int, vector<float>> scanIntensitiesByMs3Mz{};

    //Issue 295: summarize each ms3 m/z values into a single intensity
    //
    //ms3MzKey
    map<int, float> intensityByMs3Mz{};

    float sumMs3MzIntensity = 0.0f;

    // ---------------------------------------------------- //

    // SUMMARIES BASED ON MS2 m/z KEY

    // <ms2 m/z>
    map<int, unordered_set<int>> matchingCoordsByMs2{};
    map<int, int> ms3MatchesByMs2Mz{};
    map<int, float> sumMs3IntensityByMs2Mz{};
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

/**
 * @brief The DIPipelineSampleData struct
 * container for a single DI sample,
 * which encompasses a search vs. spiked-in standards
 * and a search vs. a regular library.
 *
 * Results from these two searches are combined to perform normalizations.
 */
struct DIPipelineSampleData {

    //Data

    mzSample* sample = nullptr;
    vector<Scan*> validMs1Scans = {};
    Fragment* ms1Fragment = nullptr;
    map<int, vector<Scan*>> ms2ScansByPrecursorRangeId = {};

    //Search Results

    //is results
    long searchNumOutputRows = 0;
    map<int, DirectInfusionAnnotation*> searchAnnotationsByPrecursorRangeId = {};

    //search results
    long isNumOutputRows = 0;
    map<int, DirectInfusionAnnotation*> isAnnotationsByPrecursorRangeId = {};

    //Internal Standards Normalization

    //precursor
    map<pair<string, string>, float> precursorQuantNormalizationIntensityMap = {};
    map<pair<string, string>, float> precursorQuantNormalizationMzMap = {};

    //fragments
    map<tuple<string, string, string>, float> fragmentQuantNormalizationMap = {};
};
