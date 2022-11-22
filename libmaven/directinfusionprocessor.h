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

    //Issue 313: reserved map key for compounds that cannot have any MS2 information
    static inline int getNoMs2ScansMapKey(){return -1;}

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

     //Issue 319:
     //Every adduct will become a column in a matrix of observed intensities.
     set<Adduct*> allAdducts{};

     //Issue 365:
     //use class-specific adducts for adducts table
     map<string, set<Adduct*>> adductsByClass{};

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
class DirectInfusionSearchParameters : public SearchParameters, public LipidSearchParameters {

public:

    /** ===================
     * MS1 SEARCH RELATED
     * @param ms1IsRequireAdductPrecursorMatch:reference compound associated adduct must == query adduct
     * @param ms1IsFindPrecursorIon: only retain matches where precursor peak is found in MS1 scan(s).
     * @param ms1IsMPlusOneValidPrecursor: if the [M+0] of a precursor peak cannot be found, fall back to look for the [M+1] (13C).
     * @param ms1MinIntensity: min intensity value for an MS1 spectral peak in a consensus MS1 spectrum to be considered real.
     * @param ms1ScanFilter: consider only MS1 scans that substring match in their filterString field to this
     * @param isRequireMonoisotopic: disqualify candidate MS1 peak if there is a peak exactly one 13C-12C width behind
     * @param mMinusOnePeakMaxIntensityFraction: Avoid disqualifying candidate MS1 peak
     *          if candidate monoisotopic peak intensity / candidate MS1 peak intensity <= this fraction
     * @param ms1MinScanIntensity: min intensity value for an MS1 spectral peak in an MS1 scan to be considered real.
     *          Note the difference between this parameter and 'scanFilterMinIntensity',
     *          which filters spectral peaks based on intensity as a part of consensus spectrum formation.
     * ==================== */

    bool ms1IsRequireAdductPrecursorMatch = true;
    bool ms1IsFindPrecursorIon = false;
    bool ms1IsMPlusOneValidPrecursor = false;
    float ms1MinIntensity = 0;
    string ms1ScanFilter = "";
    bool ms1IsRequireMonoisotopic = true;
    float ms1MMinusOnePeakMaxIntensityFraction = 1.0f;
    float ms1MinScanIntensity = 0;

    /** ===================
     * MS2 SEARCH RELATED
     * labels
     * If a fragment label starts with any of these substrings, it is flagged with the appropriate type.
     * Once a fragment label encounters a character that is not covered by any of the fragment labels,
     * all labels have been assigned and the check for labeling fragments stops.
     * @param ms2DiagnosticFragmentLabelTag: label indicates this ms2 fragment is diagnostic. [DEPRECATED]
     * @param ms2sn1FragmentLabelTag: label indicates this ms2 fragment is associated with an sn1 acyl chain. [DEPRECATED]
     * @param ms2sn2FragmentLabelTag: label indicates this ms2 fragment is associated with an sn2 acyl chain. [DEPRECATED]
     * @param ms2sn1MinNumMatches: min num sn1-associated reference peaks found in observed spectrum.
     * @param ms2sn2MinNumMatches: min num sn2-associated reference peaks found in observed spectrum.
     * @param ms2IsRequirePrecursorMatch: Require that a fragment m/z matched in the MS2 spectrum matches within MS2 tolerance to the compound MS1 precursor m/z.
     * ==================== */
    string ms2DiagnosticFragmentLabelTag = "*"; //here instead of LipidSearchParameters only because deprecated
    string ms2sn1FragmentLabelTag = "@"; //here instead of LipidSearchParameters only because deprecated
    string ms2sn2FragmentLabelTag = "$"; //here instead of LipidSearchParameters only because deprecated
    bool ms2IsRequirePrecursorMatch = false; //Issue 390

    /** ===================
     * MS3 SEARCH RELATED
     * @param ms3IsMs3Search: if experimental data contains MS3 scans
     * @param ms3MinNumMatches: Minimum number of reference MS3 peaks found in observed MS3 scans
     * @param ms3MinNumMs3MzMatches: Minimum number of reference MS3 m/zs found in observed MS3 scans, regardless of precursor m/zs
     * @param ms3AnalysisMs1PrecursorPpmTolr: m/z tolerance used for matching reference MS1 <--> MS3 scan precursor m/z
     * @param ms3PrecursorPpmTolr: m/z tolerance value used for matching reference MS2 m/z <--> MS3 scan precursor m/z
     * @param ms3PpmTolr: m/z tolerance value used for matching reference <--> observed spectral peaks in MS3 spectrum
     * @param ms3MinNumScans: ms3 fragment peak must be found in this many ms3 scans to count as a match
     * @param ms3MinFractionScans: ms3 fragment peak must be found in this proportion of all appropriate ms3 scans to count as a match
     * ==================== */
    bool ms3IsMs3Search = false;
    int ms3MinNumMatches = 1;
    int ms3MinNumMs3MzMatches = 1;
    float ms3AnalysisMs1PrecursorPpmTolr = 20;
    float ms3PrecursorPpmTolr = 20;
    float ms3MatchTolrInDa = 0.5f;
    float ms3MinIntensity = 0;
    int ms3MinNumScans = 0;
    float ms3MinFractionScans = 0.0f;
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
     * @param minNumScansNearestScanNormalizedIntensity: Minimum number of scans required to use this
     *       class of scans for nearest scan normalized intensity quant. Classes of scans are defined by
     *       the combination of the m/z width of the window, and the distance in scan num between scans.
     * @param minNumScansMs1ScanIntensity: Minimum number of scans require dto use this class of scans
     *       for scans ms1 observed intensity.  Classes of scans are defined by the m/z width of the window
     *       and the proximity of the query scan to the center of the window.
     * @param normClassMap: lipid class to normalize analyte to when no IS of a lipid class is available.
     *       <key, value>: <class_with_no_IS, class_with_IS>
     *       Operates as a backup, so will default to looking for an IS of the correct lipid class first.
     * ==================== */
    vector<string> ms1PartitionIntensityByFragments{"sn1","sn2"};
    bool isPreferSmallestScanMassWindow = true;
    int minNumScansNearestScanNormalizedIntensity = 2;
    int minNumScansMs1ScanIntensity = 2;
    map<string, string> normClassMap{{"Alkyl_PC","Alkenyl_PC"},{"Alkyl_PE","Alkenyl_PE"}};

    /** ===================
      * PARTITION INTENSITY COMPUTATION
      * @param partFragMinNumScans: minimum number of scans to trust fragment measurement.
      * @param partFragMaxCV: maximum CV (from scan measurements) to trust fragment measurement.
      * @param partFragMinIntensity: minimum absolute intensity value of consensus fragment intensity measurement
      *       to trust for use in partitioning.
      * =================== */
    int partFragMinNumScans = 0;
    float partFragMaxCV = 0;
    float partFragMinIntensity = 0;

    //Issue 270
    bool isReduceBySimpleParsimony = false;

    bool isDiagnosticFragmentMapAgreement(map<string, int> observedNumDiagnosticMatchesMap);

    void printParams();

    string encodeParams();

    /**
     * @brief addMs2MinNumDiagnosticMatchesMap
     * @param directInfusionSearchParameters
     * @param encodedMs2MinNumDiagnosticMatchesMap
     *
     * @deprecated
     * TODO: remove this if it starts to get in the way.
     */
    static void addMs2MinNumDiagnosticMatchesMap(shared_ptr<DirectInfusionSearchParameters> directInfusionSearchParameters, string encodedMs2MinNumDiagnosticMatchesMap);

    static shared_ptr<DirectInfusionSearchParameters> decode(string encodedParams);
};

struct ScanQuantOutput {

    bool isValid = false;
    bool isSummarized = false; // to handle combination of multiple ScanQuantOutputs during summarization

    float intensity = 0.0f;
    int numMeasurements = 0.0f;
    float medianAbsoluteDeviation = 0.0f;
    int scanDiff = 0;
    int scanWidth = 0;

    // <minMz, maxMz>
    pair<int, int> scanMzRange = make_pair(0,0);
    inline string getScanMzRangeString() {return "(" + to_string(scanMzRange.first) + ", " + to_string(scanMzRange.second) + ")";}

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
    float observedMs1Intensity = 0.0f;

    //Issue 309, 393
    ScanQuantOutput observedMs1ScanIntensityQuant;

    //Issue 522
    //Field is not filled out unless the [M+0] intensity is missing, e.g. ion coalescence
    ScanQuantOutput observedMs1ScanIntensityQuantMPlusOne;

    //Issue 210
    int numUniqueFragments = 0;
    vector<bool> isFragmentUnique; //follows m/z-sorted Compound* fragment vectors

    //Issue 288
    int ms1IntensityCoord = -1;

    //Issue 416
    string partitionGroupId;

    //Issue 416
    //lazily default to single compound ID, otherwise called by DirectInfusionMatchInformation::computeMs1PartitionFractions2()
    string getPartitionGroupId() {
        if (partitionGroupId.empty()) {
            stringstream s;
            s << std::fixed << setprecision(4);
            s << adduct->computeAdductMass(compound->getExactMass()) << "_" << static_cast<long>(round(observedMs1ScanIntensityQuant.intensity));
            partitionGroupId = s.str();
        }
        return partitionGroupId;
    }

    //Issue 486
    float acylChainPartitionFraction = 1.0f;
    float acylChainPartitionFractionSAF = 1.0f;
    float diagnosticPartitionFraction = 1.0f;
    float diagnosticPartitionFractionSAF = 1.0f;

    //Issue 492: new quant types
    float acylFragmentSum = -1.0f;
    float acylFragmentSumSAF = -1.0f;
    float diagnosticFragmentSum = -1.0f;
    float diagnosticFragmentSumSAF = -1.0f;
};

/**
 * @brief The DirectInfusionMatchAssessment struct
 *
 * as returned by DirectInfusionProcessor::assessMatch()
 */
struct DirectInfusionMatchAssessment {

    //ms1-associated
    //consensus ms1 spectrum
    float observedMs1Intensity = 0.0f;
    int ms1IntensityCoord = -1;

    //ms1 scans
    ScanQuantOutput observedMs1ScanIntensityQuant;

    //Issue 522
    //Field is not filled out unless the [M+0] intensity is missing, e.g. ion coalescence
    ScanQuantOutput observedMs1ScanIntensityQuantMPlusOne;

    //ms2-associated
    FragmentationMatchScore fragmentationMatchScore;
    map<string, int> diagnosticFragmentMatchMap = {};
    float fragmentMaxObservedIntensity = 0.0f;


    //Issue 303:
    //If this is ever true, disqualify the match.
    //Note that matches can be disqualified for other reasons, even if this flag remains false.
    bool isDisqualifyThisMatch = false;

    //Issue 313: handles various special cases
    //match information writes to ms2-associated fields
    void computeMs2MatchAssessment(Fragment *observedSpectrum,
                                   const Compound* compound,
                                   const shared_ptr<DirectInfusionSearchParameters> params,
                                   const bool debug);

    //Issue 313: for a given fragmentLabel, this method returns a vector of all labels associated
    //with the fragmentLabel.
    //the strings match the names of the parameter in the DirectInfusionSearchParameters:
    //
    //ms2DiagnosticFragmentLabel
    //ms2sn1FragmentLabel
    //ms2sn2FragmentLabel
    //
    vector<string> static getFragmentLabelTags(string fragmentLabel,
                                                 const shared_ptr<DirectInfusionSearchParameters> params,
                                                 const bool debug);

    string static getFragmentLabelWithoutTags(string fragmentLabel,
                                         const shared_ptr<DirectInfusionSearchParameters> params,
                                         const bool debug);
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

struct PartitionInformation {
    map<shared_ptr<DirectInfusionMatchData>, float> partitionFractions{};

    //These are compounds that have an MS1 m/z shared by multiple compounds,
    //and the partition fragments associated with those compounds are not all
    //unambiguous. This means that a compound that itself does not have any
    //ambiguous fragments may still have an adjusted SAF, if it shares an
    //MS1 m/z with a compound that does have ambiguous fragments.
    //If a compound has a unique MS1 m/z, it does not need SAF adjustment, and so
    //would not be included in this set.
    set<shared_ptr<DirectInfusionMatchData>> compoundsWithAdjustedSAFs{};

    map<shared_ptr<DirectInfusionMatchData>, vector<int>> matchDataToPartitionFrags{};

    //Issue 491: Add more ambiguity categories

    //Compounds with an MS1 m/z where at least one identified compound has no
    //known way to partition intensity.
    //If a compound is ever in this set, partitioning cannot be performed.
    //A compound should only ever end up in this set if it would need to be partitioned,
    //and it was found that it could not be.
    set<shared_ptr<DirectInfusionMatchData>> unpartitionableCompounds{};
    //TODO: currently unused - need to think carefully about this one,
    //If we move away from partitioning, may not be necessary
    //TODO: Issue 493

    //Compounds where at least one associated partition fragment is associated
    //with some other compound.
    //
    //This is distinct from compoundsWithAdjustedSAFs, as that is tied to the MS1 m/z.
    //In this case, the MS1 m/z is not considered at all.
    set<shared_ptr<DirectInfusionMatchData>> compoundsWithAmbiguousFragments{};

    //Map of the sum of all partition fraction-associated MS2 intensities,
    //with intensities adjusted according to SAF (when appropriate).
    //To check which compounds do not have ambiguous fragments, see compoundsWithAmbiguousFragments.
    map<shared_ptr<DirectInfusionMatchData>, float> partitionFragmentIntensitySum{};

    //m/z of fragments designated as partition fragments.  If a fragment is not a partition fragment,
    //ambiguity rules do not apply to it.
    map<int, float> partitionFragToIntensity{};
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
    map<vector<int>, vector<shared_ptr<DirectInfusionMatchData>>> fragListToCompounds{};

    //single compound ==> all identified fragment mz/s
    map<shared_ptr<DirectInfusionMatchData>, vector<int>> matchDataToFrags{};

    //single fragment m/z ==> all identified compounds containing fragment
    map<int, unordered_set<shared_ptr<DirectInfusionMatchData>>> fragToMatchData{};

    //Issue 311: get all fragments that had no fragment matches
    vector<shared_ptr<DirectInfusionMatchData>> compoundsNoFragMatches{};

    //Issue 470: keep track of compound m/zs used for acyl chain partitioning
    map<shared_ptr<DirectInfusionMatchData>, vector<int>> matchDataToAcylPartitionFrags{};

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

    //Issue 470
    string getPartitionFragmentGroupId(shared_ptr<DirectInfusionMatchData> compound, int precision=2);

    // <precursor m/z, fragment m/z> = sum of all observed ms1 scan intensity
    map<int, float> getFragToSumObservedMs1ScanIntensity(const bool debug);

    PartitionInformation getPartitionFractions(const Fragment *ms2Fragment,
                                             const shared_ptr<DirectInfusionSearchParameters> params,
                                             vector<string> fragmentLabelTags,
                                             const bool debug);
    //Issue 491
    void computeMs1PartitionFractions3(const Fragment *ms2Fragment,
                                       const shared_ptr<DirectInfusionSearchParameters> params,
                                       const bool debug);
private:

    // <precursor m/z, compound match data> = for use with partitioning intensity
    map<int, vector<shared_ptr<DirectInfusionMatchData>>> getPrecMzPartitionMap(const bool debug);

    void computeMs1PartitionFractionFromMap(
            const map<shared_ptr<DirectInfusionMatchData>, float>& totalFragIntensityByCompound,
            const float allFragIntensity,
            bool isSAF,
            const bool debug);

    void computeScanMs1PartitionFractionFromMap(
            const map<Scan*, float>& totalFragIntensityByScan,
            const map<Scan*, map<shared_ptr<DirectInfusionMatchData>, float>>& compoundFragIntensityByScan,
            bool isSAF,
            const shared_ptr<DirectInfusionSearchParameters> params,
            const bool debug);

    string getGroupId(map<shared_ptr<DirectInfusionMatchData>, vector<int>>& map, shared_ptr<DirectInfusionMatchData> compound, int precision);
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
                              const map<pair<int, int>, vector<Scan*>>& ms1Scans,
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
     static unique_ptr<DirectInfusionMatchAssessment> assessMatch(
                             const map<pair<int, int>, vector<Scan*>>& ms1Scans,
                             const Fragment *ms1Fragment,
                             const Fragment *ms2Fragment,
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

    static ScanQuantOutput findNearestScanNormalizedIntensity(const vector<Scan*>& scans,
                                                    float queryMz,
                                                    float standardMz,
                                                    shared_ptr<DirectInfusionSearchParameters> params,
                                                    int scanWidthInDa = -1, // -1 indicates all scan widths
                                                    bool debug = false);


    static map<pair<int, int>, vector<Scan*>> computeValidMs1ScansByMzRange(vector<Scan*>& validMs1Scans);

    //Issue 319: easier to factor this out
    //Issue 393: Update to respect scan types
    static ScanQuantOutput getObservedMs1ScanIntensity(
            const map<pair<int, int>, vector<Scan*>>& validMs1Scans,
            float queryMz,
            shared_ptr<DirectInfusionSearchParameters> params,
            bool debug);

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

//Issue 319: functors to compare objects instead of memory addresses in maps
//See https://stackoverflow.com/questions/25122932/pointers-as-keys-in-map-c-stl
struct adduct_less {
    bool operator()(Adduct* lhs, Adduct *rhs) const {
        return lhs->name < rhs->name;
    }
};

struct compound_less {
    bool operator()(Compound* lhs, Compound* rhs) const {
        return lhs->name < rhs->name;
    }
};

enum DISampleCompoundAdductQuantType{None=0, IdentifiedCompound=1, Reextraction=2};

//Issue 363
struct DISampleCompoundAdductQuant {

    map<string,float> quantValues{};
    map<string, float> fractionValues{};

    string partitionGroupId;
    DISampleCompoundAdductQuantType diSampleCompoundAdductQuantType = DISampleCompoundAdductQuantType::None;

    DISampleCompoundAdductQuant(){}

    DISampleCompoundAdductQuant(string partitionGroupId, DISampleCompoundAdductQuantType diSampleCompoundAdductQuantType) {
        this->partitionGroupId = partitionGroupId;
        this->diSampleCompoundAdductQuantType = diSampleCompoundAdductQuantType;
    }

    void addQuantValue(string quantName, float quantVal) {
        if (quantValues.find(quantName) == quantValues.end()) {
            quantValues.insert(make_pair(quantName, -1.0f));
        }
        quantValues[quantName] = quantVal;
    }

    float getQuantValue(string quantName) {
        if (quantValues.find(quantName) == quantValues.end()) {
            return -1.0f;
        }
        return quantValues[quantName];
    }

    void addFractionValue(string fractionName, float fractionVal) {
        if (fractionValues.find(fractionName) == fractionValues.end()) {
            fractionValues.insert(make_pair(fractionName, -1.0f));
        }
        fractionValues[fractionName] = fractionVal;
    }

    float getFractionValue(string fractionName) {
        if (fractionValues.find(fractionName) == fractionValues.end()) {
            return -1.0f;
        }
        return fractionValues[fractionName];
    }

    void addMs1ScanIntensity(float intensity) {addQuantValue("ms1_scan_intensity", intensity);}
    void addNearestScanNormalizedIntensity(float intensity) {addQuantValue("ms1_intensity_is_nearest_scan_normalized", intensity);}
    void addMs1PartitionFraction(float fraction){addFractionValue("ms1_partition_fraction", fraction);}
    void addMs1PartitionFractionSAF(float fraction){addFractionValue("ms1_partition_fraction_SAF", fraction);}
    void addMs1AdductFraction(float fraction){addFractionValue("ms1_adduct_fraction", fraction);}

    float getMs1ScanIntensity() {return getQuantValue("ms1_scan_intensity");}
    float getNearestScanNormalizedIntensity() {return getQuantValue("ms1_intensity_is_nearest_scan_normalized");}
    float getMs1PartitionFraction(){return getFractionValue("ms1_partition_fraction");}
    float getMs1PartitionFractionSAF(){return getFractionValue("ms1_partition_fraction_SAF");}
    float getMs1AdductFraction(){return getFractionValue("ms1_adduct_fraction");}

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
    vector<Scan*> validMs1Scans{};
    Fragment* ms1Fragment = nullptr;
    map<int, vector<Scan*>> ms2ScansByPrecursorRangeId{};
    map<pair<int, int>, vector<Scan*>> validMs1ScansByMzRange{};

    //Pipeline Search Results

    //is results
    long searchNumOutputRows = 0;
    map<int, DirectInfusionAnnotation*> searchAnnotationsByPrecursorRangeId{};

    //search results
    long isNumOutputRows = 0;
    map<int, DirectInfusionAnnotation*> isAnnotationsByPrecursorRangeId{};

    //Issue 319, 363: store quant measurements for different adduct forms of identified compounds
    map<Compound*, map<Adduct*, DISampleCompoundAdductQuant, adduct_less>, compound_less> compoundQuantByAdduct{};

    //Issue 365: intensity value = sum of all (SAF_partition * ms1_scan_intensity) for all IDed adducts
    map<Compound*, float, compound_less> identifiedAdductsCompoundQuant{};

    // m/z reex key
    // <int, int>
    // <mzUtils::intKeyToMz(mz, 4), mzUtils::intKeyToMz(intensity, 10)>
    static const long MZ_REEX_MZ_MULT_FACTOR = 10L;
    static const long MZ_REEX_INTENSITY_MULT_FACTOR = 1L;

    //Issue 365: mz reex key convenient functions
    static pair<long, long> getMzReexKey(float mz, float intensity) {
        long mzKey = mzUtils::mzToIntKey(static_cast<double>(mz), MZ_REEX_MZ_MULT_FACTOR);
        long intensityKey = mzUtils::mzToIntKey(static_cast<double>(intensity), MZ_REEX_INTENSITY_MULT_FACTOR);
        return make_pair(mzKey, intensityKey);
    }

    static pair<float, float> getMzIntensityFromReexKey(pair<int, int> mzReexKey) {
        float mz = static_cast<float>(mzUtils::intKeyToMz(mzReexKey.first, MZ_REEX_MZ_MULT_FACTOR));
        float intensity = static_cast<float>(mzUtils::intKeyToMz(mzReexKey.second, MZ_REEX_INTENSITY_MULT_FACTOR));
        return make_pair(mz, intensity);
    }

    //Issue 363: key = m/z reex key of identified compounds.
    //used when considering adduct table.
    set<pair<long, long>> identifiedMzs{};

    //Issue 365: key = m/z reex key, value = all identified compounds suggesting this m/z as a reextraction
    map<pair<long, long>, set<Compound*>> reextractedMzToIdentifiedCompounds{};

    //Issue 365: key = m/z reex key, value = unidentified <Compound*, Adduct*> pair
    map<pair<long, long>, vector<pair<Compound*, Adduct*>>> reextractedMzToUnidentifiedCompoundAdduct{};

    //IS normalization info (additions @ Issue 480)
    //key = <lipidClass, lipidAdduct>
    map<pair<string, string>, float> precursorQuantNormalizationIntensityMap{};
    map<pair<string, string>, float> precursorQuantNormalizationMzMap{};
    map<pair<string, string>, float> precursorQuantRelativeFractions{};
    map<pair<string, string>, int> precursorQuantAdductIntensityRank{};

    // <lipidClass,  adducts>
    map<string, vector<string>> internalStandardAdducts{};

    //IS quant for adduct table
    map<pair<string, string>, float> precursorQuantByAdductNormalizationIntensityMap{};

    //fragments
    map<tuple<string, string, string>, float> fragmentQuantNormalizationMap{};

    //Issue 487:  add some fragment information
    // Stores the IS quant information for the chosen diagnostic fragment
    //  <compoundName, adductName, fragmentLabel>, IS observed intensity
    map<tuple<string, string, string>, float> fragmentQuantMaxDiagnosticObservedISIntensityMap{};

    //Stores the normalized diagnostic intensity information
    //    <lipidClass, adductName>, normalized fragment intensity
    map<pair<string, string>, float> fragmentQuantMaxDiagnosticNormalizedIntensityMap{};

    //    <compoundName, adductName>
    map<pair<string, string>, float> fragmentQuantAcylSumNormalizedIntensityMap{};
    map<pair<string, string>, float> fragmentQuantAcylSumSAFNormalizedIntensityMap{};
    map<pair<string, string>, float> fragmentQuantDiagnosticSumNormalizedIntensityMap{};
    map<pair<string, string>, float> fragmentQuantDiagnosticSumSAFNormalizedIntensityMap{};

    //    <compoundName, adductName>
    map<pair<string, string>, float> nearestScanNormalizedIntensityMap{};

    //Issue 492: needed for normalization downstream
    //  <lipidClass, adductName>
    map<pair<string, string>, float> diagnosticFragmentSumISMap{};
    map<pair<string, string>, float> acylChainFragmentSumISMap{};

    constexpr static const float MAP_NO_VALUE = -1.0f;

    //Issue 523: fall back to checking alternative key in some cases
    float getClassAdductMapValue(pair<string, string>& class_adduct_key,
                            map<pair<string, string>, float>& class_adduct_map,
                            shared_ptr<DirectInfusionSearchParameters> params,
                            bool debug) {

        if (class_adduct_map.find(class_adduct_key) != class_adduct_map.end()) {
            return class_adduct_map.at(class_adduct_key);
        }

        string lipidClass = class_adduct_key.first;
        string adductName = class_adduct_key.second;

        if (params->normClassMap.find(lipidClass) != params->normClassMap.end()) {

            string substituteLipidClass = params->normClassMap.at(lipidClass);
            pair<string, string> alt_class_adduct_key = make_pair(substituteLipidClass, adductName);

            if (class_adduct_map.find(alt_class_adduct_key) != class_adduct_map.end()) {
                if (debug){
                    cout << "DIPipelineSampleData::getClassAdductMapValue(): "
                         << "original key: (" << lipidClass << ", " << adductName << ") "
                         << "subsitute key: (" << substituteLipidClass << ", " << adductName << ") "
                         << "original key missing, subsitute key found with "
                         << "intensity=" << class_adduct_map.at(alt_class_adduct_key)
                         << endl;
                }
                return class_adduct_map.at(alt_class_adduct_key);
            }

        }

        return MAP_NO_VALUE;
    }

    float getClassAdductFragmentMapValue(
            tuple<string, string, string>& class_adduct_fragment_key,
            map<tuple<string, string, string>, float>& class_adduct_fragment_map,
            shared_ptr<DirectInfusionSearchParameters> params,
            bool debug) {

        if (class_adduct_fragment_map.find(class_adduct_fragment_key) != class_adduct_fragment_map.end()) {
            return class_adduct_fragment_map.at(class_adduct_fragment_key);
        }

        string lipidClass = get<0>(class_adduct_fragment_key);
        string adductName = get<1>(class_adduct_fragment_key);
        string fragment = get<2>(class_adduct_fragment_key);

        if (params->normClassMap.find(lipidClass) != params->normClassMap.end()) {

            string substituteLipidClass = params->normClassMap.at(lipidClass);
            tuple<string, string, string> alt_class_adduct_fragment_key = make_tuple(substituteLipidClass, adductName, fragment);

            if (class_adduct_fragment_map.find(alt_class_adduct_fragment_key) != class_adduct_fragment_map.end()) {
                if (debug){
                    cout << "DIPipelineSampleData::getClassAdductFragmentMapValue(): "
                         << "original key: (" << lipidClass << ", " << adductName << ", " << fragment << ") "
                         << "subsitute key: (" << substituteLipidClass << ", " << adductName << ", " << fragment << ") "
                         << "original key missing, subsitute key found with "
                         << "intensity=" << class_adduct_fragment_map.at(alt_class_adduct_fragment_key)
                         << endl;
                }
                return class_adduct_fragment_map.at(alt_class_adduct_fragment_key);
            }

        }

        return MAP_NO_VALUE;
    }

};

enum ScanIntensityType{STANDARD=0, QUERY=1};

struct ScanIntensity;
struct NearestScanIntensityPair;

struct ScanIntensity {

    Scan* scan = nullptr;
    float intensity = 0;
    ScanIntensityType scanIntensityType = ScanIntensityType::QUERY;
    int scanWidth = 0;

    ScanIntensity(Scan* scanVal, float intensityVal, ScanIntensityType scanIntensityTypeVal){
        scan = scanVal;
        intensity = intensityVal;
        scanIntensityType = scanIntensityTypeVal;
        if (scanVal) {
            scanWidth = static_cast<int>(round(scan->upperLimitMz - scan->lowerLimitMz));
        }
    }

    ScanIntensity(){}

    //assumes that queryScans and standardScans are sorted in increasing order by scannum
    static vector<NearestScanIntensityPair> matchStandardScanIntensitiesToQueryScanIntensities(
            vector<ScanIntensity> queryScans,
            vector<ScanIntensity> standardScans,
            shared_ptr<DirectInfusionSearchParameters> params,
            bool debug=false);
};

struct ScanIntensityMatch {
    ScanIntensity standardScan;
    ScanIntensity queryScan;
    int dist;
};

struct NearestScanIntensityPair {

    ScanIntensity standardScan;
    ScanIntensity queryScan;
    int dist = -1;

    NearestScanIntensityPair(ScanIntensity standardScanVal, ScanIntensity queryScanVal) {
        standardScan = standardScanVal;
        queryScan = queryScanVal;
        if (standardScanVal.scan && queryScanVal.scan) {
            dist = abs(standardScanVal.scan->scannum - queryScanVal.scan->scannum);
        }
    }

    float getIntensity(){return queryScan.intensity/standardScan.intensity;}
};
