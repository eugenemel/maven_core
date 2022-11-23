#include "mzSample.h"
#include "lipidsummarizationutils.h"

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
            group.computeFragPattern(params.get());
        }

        vector<pair<Compound*, FragmentationMatchScore>> scores{};

        if (debug) {
            cout << group.meanMz << "@" << group.meanRt << ":\n";
        }

        for (long pos = lb - compounds.begin(); pos < static_cast<long>(compounds.size()); pos++){

            Compound *compound = compounds[static_cast<unsigned long>(pos)];
            float precMz = compound->precursorMz;

            //stop searching when the maxMz has been exceeded.
            if (precMz > maxMz) {
                break;
            }

            if (debug) {

                cout << "(" << minMz << ", " << maxMz << ")\n";

                if (pos >= 2){
                    cout << "compounds[" << (pos-2) << "]: " << compounds[pos-2]->precursorMz << "\n";
                    cout << "compounds[" << (pos-1) << "]: " << compounds[pos-1]->precursorMz << "\n";
                }

                cout << "compounds[" << (pos) << "]: " << precMz << " <--> " << compound->id << "\n";

                if (pos <= compounds.size()-2) {
                    cout << "compounds[" << (pos+1) << "]: " << compounds[pos+1]->precursorMz << "\n";
                    cout << "compounds[" << (pos+2) << "]: " << compounds[pos+1]->precursorMz << "\n";
                }

                cout << "\n\n";

            }

            Fragment library;
            library.precursorMz = static_cast<double>(precMz);
            library.mzs = compound->fragment_mzs;
            library.intensity_array = compound->fragment_intensity;
            library.fragment_labels = compound->fragment_labels;

            //lazily enumerate summarization info, cache for future comparisons
            if (compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) == compound->metaDataMap.end()) {
                LipidNameComponents lipidNameComponents = LipidSummarizationUtils::getNameComponents(compound->name);
                compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getLipidClassSummaryKey(), lipidNameComponents.lipidClass));
            }

            //skip entries when the RT is out of range
            string lipidClass = compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];
            if (params->lipidClassToRtRange.find(lipidClass) != params->lipidClassToRtRange.end()) {
                pair<float, float> validRtRange = params->lipidClassToRtRange[lipidClass];

                if (group.medianRt() < validRtRange.first || group.medianRt() > validRtRange.second) continue;
            }

            Fragment observed = group.fragmentationPattern;
            float maxDeltaMz = (params->ms2PpmTolr * precMz)/ 1000000;

            FragmentationMatchScore s;

            vector<int> ranks = Fragment::findFragPairsGreedyMz(&library, &observed, maxDeltaMz);

            for (unsigned int i = 0; i < ranks.size(); i++) {

                int observedIndex = ranks[i];

                if (observedIndex != -1) {

                    s.numMatches++;

                    string compoundLabel = compound->fragment_labels[i];

                    s.addLabelSpecificMatches(compoundLabel);
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

        if (debug) {
            cout << "\n";
        }

        //based on scores, determine a result
        //TODO: compare scores, summarization, etc

        //this is a dummy for testing
        float maxScore = -1.0f;
        pair<Compound*, FragmentationMatchScore> bestPair;
        for (auto score : scores) {

            if (debug) {
                cout << score.first->id << ", score=" << score.second.hypergeomScore << "\n";
            }

            if (score.second.hypergeomScore > maxScore) {
                bestPair = score;
                maxScore = score.second.hypergeomScore;
            }
        }

        if (maxScore > -1.0f) {
            group.compound = bestPair.first;
            group.fragMatchScore = bestPair.second;
            group.fragMatchScore.mergedScore = bestPair.second.hypergeomScore;

            if (debug) {
                cout << "MATCH: " << group.meanMz << "@" << group.meanRt  << " <--> " << group.compound->id << "\n" << endl;
            }
        }
    }
}

/**
 * @brief MzKitchenProcessor::matchMetabolites
 * @param groups
 * @param compounds
 * @param params
 * @param debug
 */
void MzKitchenProcessor::matchMetabolites(
        vector<PeakGroup>& groups,
        vector<Compound*>& compounds,
        shared_ptr<MzkitchenMetaboliteSearchParameters> params,
        bool debug){

    //Issue 559: Arrange groups in order for debugging
    if (debug) {
        sort(groups.begin(), groups.end(), [](const PeakGroup& lhs, const PeakGroup& rhs){
            if (lhs.meanMz == rhs.meanMz){
                return lhs.meanRt < rhs.meanRt;
            } else {
                return lhs.meanMz < rhs.meanMz;
            }
        });
    }

    for (auto& group : groups) {

        float minMz = group.meanMz - (group.meanMz*params->ms1PpmTolr/1000000);
        float maxMz = group.meanMz + (group.meanMz*params->ms1PpmTolr/1000000);

        auto lb = lower_bound(compounds.begin(), compounds.end(), minMz, [](const Compound* lhs, const float& rhs){
            return lhs->precursorMz < rhs;
        });

        if (group.fragmentationPattern.mzs.empty()) {
            group.computeFragPattern(params.get());
        }

        vector<pair<Compound*, FragmentationMatchScore>> scores{};

        for (long pos = lb - compounds.begin(); pos < static_cast<long>(compounds.size()); pos++){

            Compound *compound = compounds[static_cast<unsigned long>(pos)];
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

            //skip entries when the RT is required, and out of range
            float rtDiff = abs(group.medianRt() - compound->expectedRt);
            if (params->rtIsRequireRtMatch & (rtDiff > params->rtMatchTolerance)) {
                continue;
            }

            FragmentationMatchScore s = library.scoreMatch(&(group.fragmentationPattern), params->ms2PpmTolr);

            if (s.numMatches < params->ms2MinNumMatches) continue;

            double normCosineScore = Fragment::normCosineScore(&library, &(group.fragmentationPattern), s.ranks);
            s.dotProduct = normCosineScore;

            scores.push_back(make_pair(compound, s));
        }

        //Issue 559: print empty groups
        if (debug) {
            cout << group.meanMz << "@" << group.medianRt() << " #ms2s=" << group.ms2EventCount << " :" << endl;
        }

        if (scores.empty()) continue;

        group.compounds = scores;

        float maxScore = -1.0f;
        pair<Compound*, FragmentationMatchScore> bestPair;
        for (auto score : scores) {
            if (score.second.dotProduct > maxScore) {
                bestPair = score;
                maxScore = score.second.dotProduct;
            }
        }

        if (maxScore > -1.0f) {
            group.compound = bestPair.first;
            group.fragMatchScore = bestPair.second;
            group.fragMatchScore.mergedScore = bestPair.second.dotProduct;
        }

        //Issue 546: debugging
        if (debug) {
            //cout << group.meanMz << "@" << group.medianRt() << ":" << endl;

            for (auto pair : group.compounds) {
                cout << "\t" << pair.first->name << " "
                     << pair.second.numMatches << " "
                     << pair.second.fractionMatched << " "
                     << pair.second.dotProduct  << " "
                     << pair.second.hypergeomScore << " "
                     << endl;
            }
        }
    }
}

string LCLipidSearchParameters::encodeParams() {

    string encodedParams = baseParams();

    encodedParams = encodedParams + "mspFilePath" + "=" + mspFilePath + ";";

    //Issue 552
    encodedParams = encodedParams + "lipidClassToRtRange" + "=" + "{";

    for (auto it = lipidClassToRtRange.begin(); it != lipidClassToRtRange.end(); ++it) {
        string key = it->first;
        string value = to_string(it->second.first) + TUPLE_MAP_KEY_DELIMITER + to_string(it->second.second);
        encodedParams = encodedParams + key + "=" + value + INTERNAL_MAP_DELIMITER;
    }

    encodedParams = encodedParams + "};";

    //Issue 586
    encodedParams = encodedParams + "ms2MinNumAcylMatches" + "=" + to_string(ms2MinNumAcylMatches) + ";";

    encodedParams = encodedParams + "ms2MinNumAcylMatchesByLipidClassAndAdduct" + "=" + encodeByLipidToClassAndAdductToIntMap(
                ms2MinNumAcylMatchesByLipidClassAndAdduct,
                TUPLE_MAP_KEY_DELIMITER,
                INTERNAL_MAP_DELIMITER);

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
    if (decodedMap.find("lipidClassToRtRange") != decodedMap.end()) {

        string encodedLipidClassToRtRange = decodedMap["lipidClassToRtRange"];
        unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedLipidClassToRtRange, INTERNAL_MAP_DELIMITER);

        for (auto it = decodedMap.begin(); it != decodedMap.end(); ++it){
            string key = it->first;
            string valueEncoded = it->second;

            auto pos = valueEncoded.find(TUPLE_MAP_KEY_DELIMITER);

            string startRtEncoded = valueEncoded.substr(0, pos);
            string endRtEncoded = valueEncoded.substr(pos+1, valueEncoded.size());

            try {
                float startRt = stof(startRtEncoded);
                float endRt = stof(endRtEncoded);

                if (startRt < endRt) {
                    lipidSearchParameters->lipidClassToRtRange.insert(make_pair(key, make_pair(startRt, endRt)));
                }

            } catch (...) {
                //skip this entry - something wrong with the encoding
            }
        }
    }

    if (decodedMap.find("ms2MinNumAcylMatches") != decodedMap.end()) {
        lipidSearchParameters->ms2MinNumAcylMatches = stoi(decodedMap["ms2MinNumAcylMatches"]);
    }

    if (decodedMap.find("ms2MinNumAcylMatchesByLipidClassAndAdduct") != decodedMap.end()) {
        string encodedMs2MinNumAcylMatchesByLipidClassAndAdduct = decodedMap["ms2MinNumAcylMatchesByLipidClassAndAdduct"];
        //TODO: decoding
    }

    // END LCLipidSearchParameters

    return lipidSearchParameters;
}

string MzkitchenMetaboliteSearchParameters::encodeParams() {

    string encodedParams = baseParams();

    encodedParams = encodedParams + "mspFilePath" + "=" + mspFilePath + ";";
    encodedParams = encodedParams + "rtIsRequireRtMatch" + "=" + to_string(rtIsRequireRtMatch) + ";";
    encodedParams = encodedParams + "rtMatchTolerance" + "=" + to_string(rtMatchTolerance) + ";";

    //Issue 547
    string matchingPolicyStr = "SINGLE_TOP_HIT";
    if (matchingPolicy == PeakGroupCompoundMatchingPolicy::ALL_MATCHES) {
        matchingPolicyStr = "ALL_MATCHES";
    } else if (matchingPolicy == PeakGroupCompoundMatchingPolicy::TOP_SCORE_HITS) {
        matchingPolicyStr = "TOP_SCORE_HITS";
    }
    encodedParams = encodedParams + "matchingPolicy" + "=" + matchingPolicyStr + ";";
    encodedParams = encodedParams + "isComputeAllFragScores" + "=" + to_string(isComputeAllFragScores) + ";";

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

    //Issue 547
    if (decodedMap.find("matchingPolicy") != decodedMap.end()) {
        if (decodedMap["matchingPolicy"] == "SINGLE_TOP_HIT") {
            metaboliteSearchParameters->matchingPolicy = PeakGroupCompoundMatchingPolicy::SINGLE_TOP_HIT;
        } else if (decodedMap["matchingPolicy"] == "ALL_MATCHES") {
            metaboliteSearchParameters->matchingPolicy = PeakGroupCompoundMatchingPolicy::ALL_MATCHES;
        } else if (decodedMap["matchingPolicy"] == "TOP_SCORE_HITS") {
            metaboliteSearchParameters->matchingPolicy = PeakGroupCompoundMatchingPolicy::TOP_SCORE_HITS;
        }
    }
    if (decodedMap.find("isComputeAllFragScores") != decodedMap.end()) {
        metaboliteSearchParameters->isComputeAllFragScores = decodedMap["isComputeAllFragScores"] == "1";
    }

    // END MzkitchenMetaboliteSearchParameters

    return metaboliteSearchParameters;

}

void MzkitchenMspSearchParameters::setLegacyPeakGroupParameters(){

    //Fragment() constructor
    scanFilterMinFracIntensity = 0.01f;
    scanFilterMinSNRatio = 1.0f;
    scanFilterMaxNumberOfFragments = 1024;
    scanFilterBaseLinePercentile = 5;
    scanFilterIsRetainFragmentsAbovePrecursorMz = false;
    scanFilterPrecursorPurityPpm=10;
    scanFilterMinIntensity=0;

    //Fragment::buildConsensus()
    consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
    consensusIsIntensityAvgByObserved = false;
    consensusIsNormalizeTo10K = true;
    consensusMinNumMs2Scans=0;
    consensusMinFractionMs2Scans=0.0f;
    consensusIsRetainOriginalScanIntensities=false;
}

MzkitchenMspSearchParameters::~MzkitchenMspSearchParameters() = default;
