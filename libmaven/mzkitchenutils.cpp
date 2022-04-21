#include "mzSample.h"

/**
 * @brief LCLipidProcessor::matchLipids
 * @param groups
 * @param compounds
 * @param params
 * @param debug
 *
 * This function assumes that all peak groups contain MS2s.
 *
 * This function assumes that the compounds are lipids, they have precursorMz values,
 * they have labeled fragments, and the labels follow these reserved character rules:
 *
 * Lipid Fragment Labels:
 * sn1: @
 * sn2: $
 * sn3: !
 * sn4: ^
 * diagnostic: *
 * oxidation: &
 */
void MzKitchenProcessor::matchLipids_LC(
        vector<PeakGroup>& groups,
        vector<Compound*>& compounds,
        shared_ptr<LCLipidSearchParameters> params,
        bool debug){

    for (auto& group : groups) {

        float minMz = group.meanMz - (group.meanMz*params->ms1PpmTolr/1000000);
        float maxMz = group.meanMz + (group.meanMz*params->ms1PpmTolr/1000000);

        auto lb = lower_bound(compounds.begin(), compounds.end(), minMz, [](const Compound* lhs, const float& rhs){
            return lhs->precursorMz < rhs;
        });

        if (group.fragmentationPattern.mzs.empty()) {

            //TODO: more flexibility around consensus spectrum formation
            group.computeFragPattern(params->ms2PpmTolr);
        }

        vector<pair<Compound*, FragmentationMatchScore>> scores{};

        for (unsigned int pos = lb - compounds.begin(); pos < compounds.size(); pos++){

            Compound *compound = compounds[pos];
            float precMz = compound->precursorMz;

            //stop searching when the maxMz has been exceeded.
            if (precMz > maxMz) {
                break;
            }

            Fragment library;
            library.precursorMz = static_cast<double>(precMz);
            library.mzs = compound->fragment_mzs;
            library.intensity_array = compound->fragment_intensity;
            library.fragment_labels = compound->fragment_labels;

            Fragment observed = group.fragmentationPattern;
            float maxDeltaMz = (params->ms2PpmTolr * precMz)/ 1000000;

            FragmentationMatchScore s;

            vector<int> ranks = Fragment::findFragPairsGreedyMz(&library, &observed, maxDeltaMz);

            for (unsigned int i = 0; i < ranks.size(); i++) {

                int observedIndex = ranks[i];

                if (observedIndex != -1) {

                    s.numMatches++;

                    string compoundLabel = compound->fragment_labels[i];

                    //check labels
                    for (char c : compoundLabel){

                        if (c == '*'){
                            s.numDiagnosticMatches++;
                        } else if (c == '@') {
                            s.numSn1Matches++;
                        } else if (c == '$'){
                            s.numSn2Matches++;
                        } else if (c == '!') {
                            s.numSn3Matches++;
                        } else if (c == '^') {
                            s.numSn4Matches++;
                        } else if (c == '&') {
                            s.numOxidations++;
                        } else {
                            //finished with special labels for this fragment
                            break;
                        }

                    }
                }
            }

            if (s.numMatches < params->ms2MinNumMatches) continue;
            if (s.numDiagnosticMatches < params->ms2MinNumDiagnosticMatches) continue;

            s.hypergeomScore = Fragment::SHP(static_cast<int>(s.numMatches),
                                             static_cast<int>(library.mzs.size()),
                                             static_cast<int>(observed.nobs()),
                                             100000);

            s.dotProduct = library.dotProduct(&observed);

            s.fractionMatched = s.numMatches/library.mzs.size();

            scores.push_back(make_pair(compound, s));
        }

        //based on scores, determine a result
        //TODO: compare scores, summarization, etc

        //this is a dummy for testing
        float maxScore = -1.0f;
        pair<Compound*, FragmentationMatchScore> bestPair;
        for (auto score : scores) {
            if (score.second.hypergeomScore > maxScore) {
                bestPair = score;
                maxScore = score.second.hypergeomScore;
            }
        }

        if (maxScore > -1.0f) {
            group.compound = bestPair.first;
            group.fragMatchScore = bestPair.second;
            group.fragMatchScore.mergedScore = bestPair.second.hypergeomScore;
        }
    }
}

void MzKitchenProcessor::matchMetabolites(
        vector<PeakGroup>& groups,
        vector<Compound*>& compounds,
        shared_ptr<MzkitchenMetaboliteSearchParameters> params,
        bool debug){

}

string LCLipidSearchParameters::encodeParams() {
    string encodedParams;

    //START SearchParameters

    //program level
    encodedParams = encodedParams + "searchVersion" + "=" + searchVersion + ";";

    //scan filter params
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

    //consensus spectrum params
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

    // ms1 matching
    encodedParams = encodedParams + "ms1PpmTolr" + "=" + to_string(ms1PpmTolr) + ";";

    // ms2 search
    encodedParams = encodedParams + "ms2MinNumMatches" + "=" + to_string(ms2MinNumMatches) + ";";
    encodedParams = encodedParams + "ms2MinNumDiagnosticMatches" + "=" + to_string(ms2MinNumDiagnosticMatches) + ";";
    encodedParams = encodedParams + "ms2MinNumUniqueMatches" + "=" + to_string(ms2MinNumUniqueMatches) + ";";
    encodedParams = encodedParams + "ms2PpmTolr" + "=" + to_string(ms2PpmTolr) + ";";
    encodedParams = encodedParams + "ms2MinIntensity" + "=" + to_string(ms2MinIntensity) + ";";

    // END SearchParameters
    // START LCLipidSearchParameters

    encodedParams = encodedParams + "mspFilePath" + "=" + mspFilePath + ";";

    // END LCLipidSearchParameters

    return encodedParams;
}

shared_ptr<LCLipidSearchParameters> LCLipidSearchParameters::decode(string encodedParams){
    shared_ptr<LCLipidSearchParameters> lipidSearchParameters = shared_ptr<LCLipidSearchParameters>(new LCLipidSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    //START SearchParameters

    //program level
    if (decodedMap.find("searchVersion") != decodedMap.end()) {
        lipidSearchParameters->searchVersion = decodedMap["searchVersion"];
    }

    //scan filter params
    if (decodedMap.find("scanFilterMinFracIntensity") != decodedMap.end()){
        lipidSearchParameters->scanFilterMinFracIntensity = stof(decodedMap["scanFilterMinFracIntensity"]);
    }
    if (decodedMap.find("scanFilterMinSNRatio") != decodedMap.end()){
        lipidSearchParameters->scanFilterMinSNRatio = stof(decodedMap["scanFilterMinSNRatio"]);
    }
    if (decodedMap.find("scanFilterMaxNumberOfFragments") != decodedMap.end()) {
        lipidSearchParameters->scanFilterMaxNumberOfFragments = stoi(decodedMap["scanFilterMaxNumberOfFragments"]);
    }
    if (decodedMap.find("scanFilterBaseLinePercentile") != decodedMap.end()) {
        lipidSearchParameters->scanFilterBaseLinePercentile = stoi(decodedMap["scanFilterBaseLinePercentile"]);
    }
    if (decodedMap.find("scanFilterIsRetainFragmentsAbovePrecursorMz") != decodedMap.end()) {
        lipidSearchParameters->scanFilterIsRetainFragmentsAbovePrecursorMz = decodedMap["scanFilterIsRetainFragmentsAbovePrecursorMz"] == "1";
    }
    if (decodedMap.find("scanFilterPrecursorPurityPpm") != decodedMap.end()){
        lipidSearchParameters->scanFilterPrecursorPurityPpm = stof(decodedMap["scanFilterPrecursorPurityPpm"]);
    }
    if (decodedMap.find("scanFilterMinIntensity") != decodedMap.end()){
        lipidSearchParameters->scanFilterMinIntensity = stof(decodedMap["scanFilterMinIntensity"]);
    }

    //scan filter for MS1 scans
    if (decodedMap.find("scanFilterMs1MinRt") != decodedMap.end()) {
        lipidSearchParameters->scanFilterMs1MinRt = stof(decodedMap["scanFilterMs1MinRt"]);
    }
    if (decodedMap.find("scanFilterMs1MaxRt") != decodedMap.end()) {
        lipidSearchParameters->scanFilterMs1MaxRt = stof(decodedMap["scanFilterMs1MaxRt"]);
    }

    //scan filter for MS2 scans
    if (decodedMap.find("scanFilterMs2MinRt") != decodedMap.end()) {
        lipidSearchParameters->scanFilterMs2MinRt = stof(decodedMap["scanFilterMs2MinRt"]);
    }
    if (decodedMap.find("scanFilterMs2MaxRt") != decodedMap.end()) {
        lipidSearchParameters->scanFilterMs2MaxRt = stof(decodedMap["scanFilterMs2MaxRt"]);
    }

    //consensus spectrum params (all ms levels)

    if (decodedMap.find("consensusIsIntensityAvgByObserved") != decodedMap.end()){
        lipidSearchParameters->consensusIsIntensityAvgByObserved = decodedMap["consensusIsIntensityAvgByObserved"] == "1";
    }
    if (decodedMap.find("consensusIsNormalizeTo10K") != decodedMap.end()){
        lipidSearchParameters->consensusIsNormalizeTo10K = decodedMap["consensusIsNormalizeTo10K"] == "1";
    }
    if (decodedMap.find("consensusIntensityAgglomerationType") != decodedMap.end()) {
        string consensusIntensityAgglomerationTypeStr = decodedMap["consensusIntensityAgglomerationType"];
        if (consensusIntensityAgglomerationTypeStr == "MEAN") {
            lipidSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
        } else if (consensusIntensityAgglomerationTypeStr == "MEDIAN") {
            lipidSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
        }
    }

    //ms1 consensus spectrum params
    if (decodedMap.find("consensusMs1PpmTolr") != decodedMap.end()){
        lipidSearchParameters->consensusMs1PpmTolr = stof(decodedMap["consensusMs1PpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs1Scans") != decodedMap.end()){
        lipidSearchParameters->consensusMinNumMs1Scans = stoi(decodedMap["consensusMinNumMs1Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs1Scans") != decodedMap.end()){
        lipidSearchParameters->consensusMinFractionMs1Scans = stof(decodedMap["consensusMinFractionMs1Scans"]);
    }

    //ms2 consensus spectrum params
    if (decodedMap.find("consensusPpmTolr") != decodedMap.end()){
        lipidSearchParameters->consensusPpmTolr = stof(decodedMap["consensusPpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs2Scans") != decodedMap.end()){
        lipidSearchParameters->consensusMinNumMs2Scans = stoi(decodedMap["consensusMinNumMs2Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs2Scans") != decodedMap.end()){
        lipidSearchParameters->consensusMinFractionMs2Scans = stof(decodedMap["consensusMinFractionMs2Scans"]);
    }

    // ms1 matching
    if (decodedMap.find("ms1PpmTolr") != decodedMap.end()) {
        lipidSearchParameters->ms1PpmTolr = stof(decodedMap["ms1PpmTolr"]);
    }

    // ms2 search
    if (decodedMap.find("ms2MinNumMatches") != decodedMap.end()) {
        lipidSearchParameters->ms2MinNumMatches = stoi(decodedMap["ms2MinNumMatches"]);
    }
    if (decodedMap.find("ms2MinNumDiagnosticMatches") != decodedMap.end()) {
        lipidSearchParameters->ms2MinNumDiagnosticMatches = stoi(decodedMap["ms2MinNumDiagnosticMatches"]);
    }
    if (decodedMap.find("ms2MinNumUniqueMatches") != decodedMap.end()) {
        lipidSearchParameters->ms2MinNumUniqueMatches = stoi(decodedMap["ms2MinNumUniqueMatches"]);
    }
    if (decodedMap.find("ms2PpmTolr") != decodedMap.end()) {
        lipidSearchParameters->ms2PpmTolr = stof(decodedMap["ms2PpmTolr"]);
    }
    if (decodedMap.find("ms2MinIntensity") != decodedMap.end()) {
        lipidSearchParameters->ms2MinIntensity = stof(decodedMap["ms2MinIntensity"]);
    }

    // END SearchParameters
    // START LCLipidSearchParameters

    if (decodedMap.find("mspFilePath") != decodedMap.end()) {
        lipidSearchParameters->mspFilePath = decodedMap["mspFilePath"];
    }

    // END LCLipidSearchParameters

    return lipidSearchParameters;
}

string MzkitchenMetaboliteSearchParameters::encodeParams() {

    string encodedParams;

    //START SearchParameters

    //program level
    encodedParams = encodedParams + "searchVersion" + "=" + searchVersion + ";";

    //scan filter params
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

    //consensus spectrum params
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

    // ms1 matching
    encodedParams = encodedParams + "ms1PpmTolr" + "=" + to_string(ms1PpmTolr) + ";";

    // ms2 search
    encodedParams = encodedParams + "ms2MinNumMatches" + "=" + to_string(ms2MinNumMatches) + ";";
    encodedParams = encodedParams + "ms2MinNumDiagnosticMatches" + "=" + to_string(ms2MinNumDiagnosticMatches) + ";";
    encodedParams = encodedParams + "ms2MinNumUniqueMatches" + "=" + to_string(ms2MinNumUniqueMatches) + ";";
    encodedParams = encodedParams + "ms2PpmTolr" + "=" + to_string(ms2PpmTolr) + ";";
    encodedParams = encodedParams + "ms2MinIntensity" + "=" + to_string(ms2MinIntensity) + ";";

    // END SearchParameters
    // START MzkitchenMetaboliteSearchParameters

    encodedParams = encodedParams + "mspFilePath" + "=" + mspFilePath + ";";
    encodedParams = encodedParams + "rtIsRequireRtMatch" + "=" + to_string(rtIsRequireRtMatch) + ";";
    encodedParams = encodedParams + "rtMatchTolerance" + "=" + to_string(rtMatchTolerance) + ";";

    // END MzkitchenMetaboliteSearchParameters

    return encodedParams;
}

shared_ptr<MzkitchenMetaboliteSearchParameters> MzkitchenMetaboliteSearchParameters::decode(string encodedParams){

    shared_ptr<MzkitchenMetaboliteSearchParameters> metaboliteSearchParameters = shared_ptr<MzkitchenMetaboliteSearchParameters>(new MzkitchenMetaboliteSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    //START SearchParameters

    //program level
    if (decodedMap.find("searchVersion") != decodedMap.end()) {
        metaboliteSearchParameters->searchVersion = decodedMap["searchVersion"];
    }

    //scan filter params
    if (decodedMap.find("scanFilterMinFracIntensity") != decodedMap.end()){
        metaboliteSearchParameters->scanFilterMinFracIntensity = stof(decodedMap["scanFilterMinFracIntensity"]);
    }
    if (decodedMap.find("scanFilterMinSNRatio") != decodedMap.end()){
        metaboliteSearchParameters->scanFilterMinSNRatio = stof(decodedMap["scanFilterMinSNRatio"]);
    }
    if (decodedMap.find("scanFilterMaxNumberOfFragments") != decodedMap.end()) {
        metaboliteSearchParameters->scanFilterMaxNumberOfFragments = stoi(decodedMap["scanFilterMaxNumberOfFragments"]);
    }
    if (decodedMap.find("scanFilterBaseLinePercentile") != decodedMap.end()) {
        metaboliteSearchParameters->scanFilterBaseLinePercentile = stoi(decodedMap["scanFilterBaseLinePercentile"]);
    }
    if (decodedMap.find("scanFilterIsRetainFragmentsAbovePrecursorMz") != decodedMap.end()) {
        metaboliteSearchParameters->scanFilterIsRetainFragmentsAbovePrecursorMz = decodedMap["scanFilterIsRetainFragmentsAbovePrecursorMz"] == "1";
    }
    if (decodedMap.find("scanFilterPrecursorPurityPpm") != decodedMap.end()){
        metaboliteSearchParameters->scanFilterPrecursorPurityPpm = stof(decodedMap["scanFilterPrecursorPurityPpm"]);
    }
    if (decodedMap.find("scanFilterMinIntensity") != decodedMap.end()){
        metaboliteSearchParameters->scanFilterMinIntensity = stof(decodedMap["scanFilterMinIntensity"]);
    }

    //scan filter for MS1 scans
    if (decodedMap.find("scanFilterMs1MinRt") != decodedMap.end()) {
        metaboliteSearchParameters->scanFilterMs1MinRt = stof(decodedMap["scanFilterMs1MinRt"]);
    }
    if (decodedMap.find("scanFilterMs1MaxRt") != decodedMap.end()) {
        metaboliteSearchParameters->scanFilterMs1MaxRt = stof(decodedMap["scanFilterMs1MaxRt"]);
    }

    //scan filter for MS2 scans
    if (decodedMap.find("scanFilterMs2MinRt") != decodedMap.end()) {
        metaboliteSearchParameters->scanFilterMs2MinRt = stof(decodedMap["scanFilterMs2MinRt"]);
    }
    if (decodedMap.find("scanFilterMs2MaxRt") != decodedMap.end()) {
        metaboliteSearchParameters->scanFilterMs2MaxRt = stof(decodedMap["scanFilterMs2MaxRt"]);
    }

    //consensus spectrum params (all ms levels)

    if (decodedMap.find("consensusIsIntensityAvgByObserved") != decodedMap.end()){
        metaboliteSearchParameters->consensusIsIntensityAvgByObserved = decodedMap["consensusIsIntensityAvgByObserved"] == "1";
    }
    if (decodedMap.find("consensusIsNormalizeTo10K") != decodedMap.end()){
        metaboliteSearchParameters->consensusIsNormalizeTo10K = decodedMap["consensusIsNormalizeTo10K"] == "1";
    }
    if (decodedMap.find("consensusIntensityAgglomerationType") != decodedMap.end()) {
        string consensusIntensityAgglomerationTypeStr = decodedMap["consensusIntensityAgglomerationType"];
        if (consensusIntensityAgglomerationTypeStr == "MEAN") {
            metaboliteSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
        } else if (consensusIntensityAgglomerationTypeStr == "MEDIAN") {
            metaboliteSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
        }
    }

    //ms1 consensus spectrum params
    if (decodedMap.find("consensusMs1PpmTolr") != decodedMap.end()){
        metaboliteSearchParameters->consensusMs1PpmTolr = stof(decodedMap["consensusMs1PpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs1Scans") != decodedMap.end()){
        metaboliteSearchParameters->consensusMinNumMs1Scans = stoi(decodedMap["consensusMinNumMs1Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs1Scans") != decodedMap.end()){
        metaboliteSearchParameters->consensusMinFractionMs1Scans = stof(decodedMap["consensusMinFractionMs1Scans"]);
    }

    //ms2 consensus spectrum params
    if (decodedMap.find("consensusPpmTolr") != decodedMap.end()){
        metaboliteSearchParameters->consensusPpmTolr = stof(decodedMap["consensusPpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs2Scans") != decodedMap.end()){
        metaboliteSearchParameters->consensusMinNumMs2Scans = stoi(decodedMap["consensusMinNumMs2Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs2Scans") != decodedMap.end()){
        metaboliteSearchParameters->consensusMinFractionMs2Scans = stof(decodedMap["consensusMinFractionMs2Scans"]);
    }

    // ms1 matching
    if (decodedMap.find("ms1PpmTolr") != decodedMap.end()) {
        metaboliteSearchParameters->ms1PpmTolr = stof(decodedMap["ms1PpmTolr"]);
    }

    // ms2 search
    if (decodedMap.find("ms2MinNumMatches") != decodedMap.end()) {
        metaboliteSearchParameters->ms2MinNumMatches = stoi(decodedMap["ms2MinNumMatches"]);
    }
    if (decodedMap.find("ms2MinNumDiagnosticMatches") != decodedMap.end()) {
        metaboliteSearchParameters->ms2MinNumDiagnosticMatches = stoi(decodedMap["ms2MinNumDiagnosticMatches"]);
    }
    if (decodedMap.find("ms2MinNumUniqueMatches") != decodedMap.end()) {
        metaboliteSearchParameters->ms2MinNumUniqueMatches = stoi(decodedMap["ms2MinNumUniqueMatches"]);
    }
    if (decodedMap.find("ms2PpmTolr") != decodedMap.end()) {
        metaboliteSearchParameters->ms2PpmTolr = stof(decodedMap["ms2PpmTolr"]);
    }
    if (decodedMap.find("ms2MinIntensity") != decodedMap.end()) {
        metaboliteSearchParameters->ms2MinIntensity = stof(decodedMap["ms2MinIntensity"]);
    }

    // END SearchParameters
    // START MzkitchenMetaboliteSearchParameters

    if (decodedMap.find("mspFilePath") != decodedMap.end()) {
        metaboliteSearchParameters->mspFilePath = decodedMap["mspFilePath"];
    }
    if (decodedMap.find("rtIsRequireRtMatch") != decodedMap.end()) {
        metaboliteSearchParameters->rtIsRequireRtMatch = decodedMap["rtIsRequireRtMatch"]=="1";
    }
    if (decodedMap.find("rtMatchTolerance") != decodedMap.end()) {
        metaboliteSearchParameters->rtMatchTolerance = stof(decodedMap["rtMatchTolerance"]);
    }

    // END MzkitchenMetaboliteSearchParameters

    return metaboliteSearchParameters;

}

MzkitchenMspSearchParameters::~MzkitchenMspSearchParameters() = default;
