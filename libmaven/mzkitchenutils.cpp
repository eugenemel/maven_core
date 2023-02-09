#include "mzSample.h"
#include "lipidsummarizationutils.h"

///**
// * @brief LCLipidProcessor::matchLipids
// * @param groups
// * @param compounds
// * @param params
// * @param debug
// *
// * This function assumes that all peak groups contain MS2s.
// *
// * This function assumes that the compounds are lipids, they have precursorMz values,
// * they have labeled fragments, and the labels follow these reserved character rules:
// *
// * Lipid Fragment Labels:
// * sn1: @
// * sn2: $
// * sn3: !
// * sn4: ^
// * diagnostic: *
// * oxidation: &
// */
//void MzKitchenProcessor::matchLipids_LC(
//        vector<PeakGroup>& groups,
//        vector<Compound*>& compounds,
//        shared_ptr<LCLipidSearchParameters> params,
//        bool debug){

//    for (auto& group : groups) {

//        float minMz = group.meanMz - (group.meanMz*params->ms1PpmTolr/1000000);
//        float maxMz = group.meanMz + (group.meanMz*params->ms1PpmTolr/1000000);
//        float deltaMz = group.meanMz*params->ms1PpmTolr/1000000;

//        auto lb = lower_bound(compounds.begin(), compounds.end(), minMz, [](const Compound* lhs, const float& rhs){
//            return lhs->precursorMz < rhs;
//        });

//        if (group.fragmentationPattern.mzs.empty()) {
//            group.computeFragPattern(params.get());
//        }

//        vector<pair<Compound*, FragmentationMatchScore>> scores{};

//        if (debug) {
//            cout << group.meanMz << "@" << group.meanRt << ":\n";
//            cout << "tol: " << params->ms1PpmTolr
//                 << " ppm, deltaMz=" << deltaMz
//                 <<  ", search range: ["
//                 << minMz << " - " << maxMz << "]\n";
//        }

//        for (long pos = lb - compounds.begin(); pos < static_cast<long>(compounds.size()); pos++){

//            Compound *compound = compounds[static_cast<unsigned long>(pos)];
//            float precMz = compound->precursorMz;

//            if (debug) {
//                cout << compound->name << " " << compound->adductString << ": " << precMz << "\n";
//            }

//            //stop searching when the maxMz has been exceeded.
//            if (precMz > maxMz) {
//                break;
//            }

//            if (debug) {

//                cout << "(minMz=" << minMz << ", maxMz=" << maxMz << "):\n";

//                if (pos >= 2){
//                    cout << "compounds[" << (pos-2) << "]: " << compounds[pos-2]->precursorMz << "\n";
//                    cout << "compounds[" << (pos-1) << "]: " << compounds[pos-1]->precursorMz << "\n";
//                }

//                cout << "compounds[" << (pos) << "]: " << precMz << " <--> " << compound->id << "\n";

//                if (pos <= compounds.size()-2) {
//                    cout << "compounds[" << (pos+1) << "]: " << compounds[pos+1]->precursorMz << "\n";
//                    cout << "compounds[" << (pos+2) << "]: " << compounds[pos+1]->precursorMz << "\n";
//                }

//                cout << "\n";

//            }

//            Fragment library;
//            library.precursorMz = static_cast<double>(precMz);
//            library.mzs = compound->fragment_mzs;
//            library.intensity_array = compound->fragment_intensity;
//            library.fragment_labels = compound->fragment_labels;

//            //lazily enumerate summarization info, cache for future comparisons
//            if (compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) == compound->metaDataMap.end()) {
//                LipidNameComponents lipidNameComponents = LipidSummarizationUtils::getNameComponents(compound->name);
//                compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getLipidClassSummaryKey(), lipidNameComponents.lipidClass));
//            }

//            string lipidClass = compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];

//            //Issue 588: skip entries unless the (class, adduct) is explicitly permitted, or no entries provided.
//            bool isValidClassAdduct = params->validClassAdducts.empty();
//            for (auto pair : params->validClassAdducts) {
//                if (pair.first == lipidClass && (pair.second == "*" || pair.second == compound->adductString)) {
//                    isValidClassAdduct = true;
//                    break;
//                }
//            }

//            if (!isValidClassAdduct) continue;

//            //skip entries when the RT is out of range
//            if (params->lipidClassToRtRange.find(lipidClass) != params->lipidClassToRtRange.end()) {
//                pair<float, float> validRtRange = params->lipidClassToRtRange[lipidClass];

//                if (group.medianRt() < validRtRange.first || group.medianRt() > validRtRange.second) continue;
//            }

//            Fragment observed = group.fragmentationPattern;
//            float maxDeltaMz = (params->ms2PpmTolr * precMz)/ 1000000;

//            FragmentationMatchScore s;

//            vector<int> ranks = Fragment::findFragPairsGreedyMz(&library, &observed, maxDeltaMz);

//            if (debug) {
//                for (unsigned int i = 0; i < library.mzs.size(); i++) {
//                    cout << "i=" << i << ": " << library.fragment_labels[i] << " " << library.mzs[i];
//                    int observedIndex = ranks.at(i);
//                    if (observedIndex >= 0) {
//                        cout << " <--> " << observed.mzs.at(static_cast<unsigned int>(observedIndex));
//                    }
//                    cout << "\n";
//                }
//            }

//            for (unsigned int i = 0; i < ranks.size(); i++) {

//                int observedIndex = ranks[i];

//                if (observedIndex != -1) {

//                    s.numMatches++;

//                    string compoundLabel = library.fragment_labels[i];

//                    s.addLabelSpecificMatches(compoundLabel, debug);
//                }
//            }

//            //Issue 606: START OLD APPROACH
////            if (s.numMatches < params->ms2MinNumMatches) continue;
////            if (s.numDiagnosticMatches < params->ms2MinNumDiagnosticMatches) continue;
//            //Issue 606: END OLD APPROACH

//            //Issue 606: NEW APPROACH more complex parameter matching
//            if (compound->adductString == "") continue;
//            if (!params->isMatchPassesLCLipidSearchThresholds(s, lipidClass, compound->adductString)) continue;
//            //Issue 606: END NEW APPROACH

//            s.hypergeomScore = Fragment::SHP(static_cast<int>(s.numMatches),
//                                             static_cast<int>(library.mzs.size()),
//                                             static_cast<int>(observed.nobs()),
//                                             100000);

//            s.dotProduct = Fragment::normCosineScore(&library, &observed, ranks);

//            s.fractionMatched = s.numMatches/library.mzs.size();
//            s.ppmError = static_cast<double>(mzUtils::ppmDist(compound->precursorMz, group.meanMz));

//            //debugging
//            if (debug) {
//                cout << "Candidate Score: " << compound->name << " " << compound->adductString <<":\n";
//                cout << "numMatches= " << s.numMatches
//                     << ", numDiagnosticMatches= " << s.numDiagnosticMatches
//                     << ", numAcylMatches= " << s.numAcylChainMatches
//                     << ", hyperGeometricScore= " << s.hypergeomScore
//                     << ", cosineScore= " << s.dotProduct
//                     << "\n\n\n";
//            }

//            //Issue 606: relocated this call here for debugging
//            // if (!params->isMatchPassesLCLipidSearchThresholds(s, lipidClass, compound->adductString)) continue;

//            scores.push_back(make_pair(compound, s));
//        }

//        if (debug) {
//            cout << endl;
//        }

//        if (!scores.empty()) {

//            //Issue 606: Pass along lipid scores values as additional table.
//            group.compounds = scores;

//            //Issue 593: guarantee non-determinism
//            sort(scores.begin(), scores.end(), [](pair<Compound*, FragmentationMatchScore>& lhs, pair<Compound*, FragmentationMatchScore>& rhs){
//                if (lhs.second.hypergeomScore != rhs.second.hypergeomScore) {
//                    return lhs.second.hypergeomScore > rhs.second.hypergeomScore;
//                }

//                if (lhs.second.dotProduct != rhs.second.dotProduct) {
//                    return lhs.second.dotProduct > rhs.second.dotProduct;
//                }

//                return lhs.first->name < rhs.first->name;
//            });

//            if (debug) {
//                cout << "Sorted scores:\n";
//                for (auto score : scores) {
//                    cout << score.first->name << " " << score.first->adductString << ":\n";
//                    cout << "numMatches= " << score.second.numMatches
//                         << ", numDiagnosticMatches= " << score.second.numDiagnosticMatches
//                         << ", numAcylMatches= " << score.second.numAcylChainMatches
//                         << ", hyperGeometricScore= " << score.second.hypergeomScore
//                         << ", cosineScore= " << score.second.dotProduct
//                         << endl;
//                }
//            }

//            pair<Compound*, FragmentationMatchScore> bestPair = scores[0];

//            group.compound = bestPair.first;
//            group.fragMatchScore = bestPair.second;
//            group.fragMatchScore.mergedScore = bestPair.second.hypergeomScore;

//            if (debug) {
//                cout << "MATCH: " << group.meanMz << "@" << group.meanRt  << " <--> " << group.compound->id << "\n" << endl;
//            }
//        }
//    }
//}

void MzKitchenProcessor::matchLipids_LC(
        vector<PeakGroup>& groups,
        vector<Compound*>& compounds,
        shared_ptr<LCLipidSearchParameters> params,
        bool debug){

    auto compoundIons = vector<CompoundIon>(compounds.size());
    for (unsigned int i = 0; i < compounds.size(); i++) {
        compoundIons[i] = CompoundIon(compounds[i]);
    }

    for (auto& g : groups) {
        MzKitchenProcessor::assignBestLipidToGroup(
                    &g,
                    compoundIons,
                    params,
                    debug);
    }
}

/**
 * @brief MzKitchenProcessor::assignBestLipidToGroup
 * @param g
 * @param compounds
 * @param params
 * @param debug
 *
 */
void MzKitchenProcessor::assignBestLipidToGroup(
        PeakGroup *g,
        vector<CompoundIon>& compounds,
        shared_ptr<LCLipidSearchParameters> params,
        bool debug){

    float minMz = g->meanMz - (g->meanMz*params->ms1PpmTolr/1000000);
    float maxMz = g->meanMz + (g->meanMz*params->ms1PpmTolr/1000000);
    float deltaMz = g->meanMz*params->ms1PpmTolr/1000000;

    auto lb = lower_bound(compounds.begin(), compounds.end(), minMz, [](const CompoundIon& lhs, const float& rhs){
        return lhs.precursorMz < rhs;
    });

    if (g->fragmentationPattern.mzs.empty()) {
        g->computeFragPattern(params.get());
    }

    vector<pair<Compound*, FragmentationMatchScore>> scores{};

    if (debug) {
        cout << g->meanMz << "@" << g->meanRt << ":\n";
        cout << "tol: " << params->ms1PpmTolr
             << " ppm, deltaMz=" << deltaMz
             <<  ", search range: ["
             << minMz << " - " << maxMz << "]\n";
    }

    for (long pos = lb - compounds.begin(); pos < static_cast<long>(compounds.size()); pos++){

        CompoundIon ion = compounds[static_cast<unsigned long>(pos)];
        Compound *compound = compounds[static_cast<unsigned long>(pos)].compound;

        float precMz = ion.precursorMz;
        string adductName = ion.getAdductName();

        if (debug) {
            cout << compound->name << " " << adductName << ": " << precMz << "\n";
        }

        //stop searching when the maxMz has been exceeded.
        if (precMz > maxMz) {
            break;
        }

        if (debug) {

            cout << "(minMz=" << minMz << ", maxMz=" << maxMz << "):\n";

            if (pos >= 2){
                cout << "compounds[" << (pos-2) << "]: " << compounds[pos-2].precursorMz << "\n";
                cout << "compounds[" << (pos-1) << "]: " << compounds[pos-1].precursorMz << "\n";
            }

            cout << "compounds[" << (pos) << "]: " << precMz << " <--> " << compound->id << "\n";

            if (pos <= static_cast<long>(compounds.size()-2)) {
                cout << "compounds[" << (pos+1) << "]: " << compounds[pos+1].precursorMz << "\n";
                cout << "compounds[" << (pos+2) << "]: " << compounds[pos+1].precursorMz << "\n";
            }

            cout << "\n";

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

        string lipidClass = compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];

        //Issue 588: skip entries unless the (class, adduct) is explicitly permitted, or no entries provided.
        bool isValidClassAdduct = params->validClassAdducts.empty();
        for (auto pair : params->validClassAdducts) {
            if (pair.first == lipidClass && (pair.second == "*" || pair.second == adductName)) {
                isValidClassAdduct = true;
                break;
            }
        }

        if (!isValidClassAdduct) continue;

        //skip entries when the RT is out of range
        if (params->lipidClassToRtRange.find(lipidClass) != params->lipidClassToRtRange.end()) {
            pair<float, float> validRtRange = params->lipidClassToRtRange[lipidClass];

            if (g->medianRt() < validRtRange.first || g->medianRt() > validRtRange.second) continue;
        }

        Fragment observed = g->fragmentationPattern;
        float maxDeltaMz = (params->ms2PpmTolr * precMz)/ 1000000;

        FragmentationMatchScore s;

        vector<int> ranks = Fragment::findFragPairsGreedyMz(&library, &observed, maxDeltaMz);

        if (debug) {
            for (unsigned int i = 0; i < library.mzs.size(); i++) {
                cout << "i=" << i << ": " << library.fragment_labels[i] << " " << library.mzs[i];
                int observedIndex = ranks.at(i);
                if (observedIndex >= 0) {
                    cout << " <--> " << observed.mzs.at(static_cast<unsigned int>(observedIndex));
                }
                cout << "\n";
            }
        }

        for (unsigned int i = 0; i < ranks.size(); i++) {

            int observedIndex = ranks[i];

            if (observedIndex != -1) {

                s.numMatches++;

                string compoundLabel = library.fragment_labels[i];

                s.addLabelSpecificMatches(compoundLabel, debug);
            }
        }

        if (adductName == "") continue;
        if (!params->isMatchPassesLCLipidSearchThresholds(s, lipidClass, adductName)) continue;

        s.hypergeomScore = Fragment::SHP(static_cast<int>(s.numMatches),
                                         static_cast<int>(library.mzs.size()),
                                         static_cast<int>(observed.nobs()),
                                         100000);

        s.dotProduct = Fragment::normCosineScore(&library, &observed, ranks);

        s.fractionMatched = s.numMatches/library.mzs.size();
        s.ppmError = static_cast<double>(mzUtils::ppmDist(compound->precursorMz, g->meanMz));

        //debugging
        if (debug) {
            cout << "Candidate Score: " << compound->name << " " << adductName <<":\n";
            cout << "numMatches= " << s.numMatches
                 << ", numDiagnosticMatches= " << s.numDiagnosticMatches
                 << ", numAcylMatches= " << s.numAcylChainMatches
                 << ", hyperGeometricScore= " << s.hypergeomScore
                 << ", cosineScore= " << s.dotProduct
                 << "\n\n\n";
        }

        scores.push_back(make_pair(compound, s));
    }

    if (debug) {
        cout << endl;
    }

    if (!scores.empty()) {

        //Issue 606: Pass along lipid scores values as additional table.
        g->compounds = scores;

        //Issue 593: guarantee non-determinism
        sort(scores.begin(), scores.end(), [](pair<Compound*, FragmentationMatchScore>& lhs, pair<Compound*, FragmentationMatchScore>& rhs){
            if (lhs.second.hypergeomScore != rhs.second.hypergeomScore) {
                return lhs.second.hypergeomScore > rhs.second.hypergeomScore;
            }

            if (lhs.second.dotProduct != rhs.second.dotProduct) {
                return lhs.second.dotProduct > rhs.second.dotProduct;
            }

            return lhs.first->name < rhs.first->name;
        });

        if (debug) {
            cout << "Sorted scores:\n";
            for (auto score : scores) {
                cout << score.first->name << " " << score.first->adductString << ":\n";
                cout << "numMatches= " << score.second.numMatches
                     << ", numDiagnosticMatches= " << score.second.numDiagnosticMatches
                     << ", numAcylMatches= " << score.second.numAcylChainMatches
                     << ", hyperGeometricScore= " << score.second.hypergeomScore
                     << ", cosineScore= " << score.second.dotProduct
                     << endl;
            }
        }

        pair<Compound*, FragmentationMatchScore> bestPair = scores[0];

        g->compound = bestPair.first;
        g->fragMatchScore = bestPair.second;
        g->fragMatchScore.mergedScore = bestPair.second.hypergeomScore;

        if (debug) {
            cout << "MATCH: " << g->meanMz << "@" << g->meanRt  << " <--> " << g->compound->id << "\n" << endl;
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
    encodedParams = encodedParams + getEncodedLipidParameters(TUPLE_MAP_KEY_DELIMITER, INTERNAL_MAP_DELIMITER);

    return encodedParams;
}

shared_ptr<LCLipidSearchParameters> LCLipidSearchParameters::decode(string encodedParams){
    shared_ptr<LCLipidSearchParameters> lipidSearchParameters = shared_ptr<LCLipidSearchParameters>(new LCLipidSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    lipidSearchParameters->fillInBaseParams(decodedMap);
    lipidSearchParameters->fillInLipidParameters(decodedMap, TUPLE_MAP_KEY_DELIMITER, INTERNAL_MAP_DELIMITER);

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
    // END LCLipidSearchParameters

    return lipidSearchParameters;
}

bool LCLipidSearchParameters::isMatchPassesLCLipidSearchThresholds(
        FragmentationMatchScore &s,
        string lipidClass,
        string adductName,
        bool debug){

    return isMatchPassesLipidSearchThresholds(
                s,
                lipidClass,
                adductName,
                ms2MinNumMatches,
                ms2MinNumDiagnosticMatches,
                ms2IsRequirePrecursorMatch,
                debug);
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

    metaboliteSearchParameters->fillInBaseParams(decodedMap);

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

