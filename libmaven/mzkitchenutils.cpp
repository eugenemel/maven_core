#include "mzSample.h"
#include "lipidsummarizationutils.h"

/**
 * @brief LCLipidProcessor::matchLipids_LC
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

    if (!g) return;

    float peakGroupMz = g->meanMz;
    float peakGroupRt = g->meanRt;

    if (params->isUseGroupMaxPeakVals) {
        peakGroupMz = g->maxPeakMzVal;
        peakGroupRt = g->maxPeakRtVal;
    }

    float minMz = peakGroupMz - (peakGroupMz*params->ms1PpmTolr/1000000);
    float maxMz = peakGroupMz + (peakGroupMz*params->ms1PpmTolr/1000000);
    float deltaMz = peakGroupMz*params->ms1PpmTolr/1000000;


    auto lb = lower_bound(compounds.begin(), compounds.end(), minMz, [](const CompoundIon& lhs, const float& rhs){
        return lhs.precursorMz < rhs;
    });

    if (g->fragmentationPattern.mzs.empty()) {
        g->computeFragPattern(params.get());
    }

    vector<tuple<CompoundIon, FragmentationMatchScore, RtAgreementState>> scores{};

    if (debug) {
        cout << peakGroupMz << "@" << peakGroupRt << ":\n";
        cout << "tol: " << params->ms1PpmTolr
             << " ppm, deltaMz=" << deltaMz
             <<  ", search range: ["
             << minMz << " - " << maxMz << "]\n";
    }

    for (long pos = lb - compounds.begin(); pos < static_cast<long>(compounds.size()); pos++){

        CompoundIon ion = compounds[static_cast<unsigned long>(pos)];
        Compound *compound = ion.compound;
        Adduct *adduct = ion.adduct;

        if (!compound) continue;

        bool isMatchingAdduct = adduct && compound->adductString == adduct->name;
        if (params->IDisRequireMatchingAdduct && !isMatchingAdduct) {
            continue;
        }

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

        // Assess RTAgreementState

        RtAgreementState rtAgreementState = RtAgreementState::COMPOUND_MISSING_RT;

        //skip entries when the RT is out of range.
        if (params->lipidClassToRtRange.find(lipidClass) != params->lipidClassToRtRange.end()) {
            pair<float, float> validRtRange = params->lipidClassToRtRange[lipidClass];

            if (g->medianRt() < validRtRange.first || g->medianRt() > validRtRange.second) {
                rtAgreementState = RtAgreementState::RT_DISAGREEMENT;
            }
        }

        //skip entries when the RT is out of range.
        if (rtAgreementState == RtAgreementState::RT_DISAGREEMENT) {
            continue;
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
        s.ppmError = static_cast<double>(mzUtils::ppmDist(compound->precursorMz, peakGroupMz));

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

        //Issue 777: MS2 score is not controlled by (class, adduct) settings for lipids
        if (s.hypergeomScore < params->ms2MinScore) continue;

        scores.push_back(make_tuple(ion, s, rtAgreementState));
    }

    if (debug) {
        cout << endl;
    }

    if (!scores.empty()) {

        //Issue 606: Pass along lipid scores values as additional table.
        g->compounds = scores;

        //Issue 593: guarantee non-determinism
        sort(scores.begin(), scores.end(), [](tuple<CompoundIon, FragmentationMatchScore, RtAgreementState>& lhs, tuple<CompoundIon, FragmentationMatchScore, RtAgreementState>& rhs){

            if (get<1>(lhs).hypergeomScore != get<1>(rhs).hypergeomScore) {
                return get<1>(lhs).hypergeomScore > get<1>(rhs).hypergeomScore;
            }

            if (get<1>(lhs).dotProduct != get<1>(rhs).dotProduct) {
                return get<1>(lhs).dotProduct > get<1>(rhs).dotProduct;
            }

            return get<0>(lhs).compound->name < get<0>(rhs).compound->name;
        });

        if (debug) {
            cout << "Sorted scores:\n";
            for (auto score : scores) {
                Compound *compound = get<0>(score).compound;
                if (compound) {
                    cout << compound->name << " " << get<0>(score).getAdductName() << ":\n";
                }
                FragmentationMatchScore fragmentationMatchScore = get<1>(score);
                cout << "numMatches= " << fragmentationMatchScore.numMatches
                     << ", numDiagnosticMatches= " << fragmentationMatchScore.numDiagnosticMatches
                     << ", numAcylMatches= " << fragmentationMatchScore.numAcylChainMatches
                     << ", hyperGeometricScore= " << fragmentationMatchScore.hypergeomScore
                     << ", cosineScore= " << fragmentationMatchScore.dotProduct
                     << endl;
            }
        }

        tuple<CompoundIon, FragmentationMatchScore, RtAgreementState> bestScore = scores[0];

        Adduct *adduct = get<0>(bestScore).adduct;
        FragmentationMatchScore fragmentationMatchScore = get<1>(bestScore);

        g->compound = get<0>(bestScore).compound;
        if (adduct) g->adduct = adduct;
        g->fragMatchScore = fragmentationMatchScore;
        g->fragMatchScore.mergedScore = fragmentationMatchScore.hypergeomScore;

        //Issue 769
        MzKitchenProcessor::labelRtAgreement(g, 'l', debug);

        //Issue 835
        if (g->compound->metaDataMap.find(MzKitchenProcessor::MZKITCHEN_NOTES_KEY) != g->compound->metaDataMap.end()) {
            g->notes = g->compound->metaDataMap.at(MzKitchenProcessor::MZKITCHEN_NOTES_KEY);
        }

        if (debug) {
            cout << "MATCH: " << peakGroupMz << "@" << peakGroupRt  << " <--> " << g->compound->id << "\n" << endl;
        }
    }

}

void MzKitchenProcessor::labelRtAgreement(PeakGroup *g, char rtMatchLabel, bool debug){

    if (!g) return;
    if (!g->compound) return;
    if (g->compound->expectedRtMin < 0) return;
    if (g->compound->expectedRtMax < 0) return;

    if (debug) {
        cout << "MzKitchenProcessor::labelRtAgreement(): Group: " << g->meanMz
             << "@" << g->meanRt
             << " (" << g->compound->name
             << "): Compound RT range [" << g->compound->expectedRtMin
             << ", " << g->compound->expectedRtMax << "]"
             << endl;
    }

    if (g->meanRt >= g->compound->expectedRtMin && g->meanRt <= g->compound->expectedRtMax) {
        g->addLabel(rtMatchLabel);
    }
}

/**
 * @brief MzKitchenProcessor::assessRtAgreement
 * Determine if a peakgroup and compound are close enough in RT space to constitute a possible match.
 * This function will intelligently determine if the compound has an RT value, if it has an RT value
 * and no RT_min/Rt_max, or if it has an RT value, RT_min and RT_max, and will carry out
 * the comparison appropriately.
 *
 * @param group
 * @param compound
 * @param ms1RtTolr: fallback tolerance for RT matching
 */
RtAgreementState MzKitchenProcessor::assessRtAgreement(float groupRtVal, Compound *compound, float ms1RtTolr, bool debug) {
    if (!compound) return RtAgreementState::COMPOUND_MISSING;

    //Issue 792: If the compound is missing a valid RT value,
    //return 'true', indicating that there is no disagreement between
    //the compound's stated RT and the peak group's measured RT.
    //This was an intentional behavior change introduced in this case.
    if (compound->expectedRt < 0) return RtAgreementState::COMPOUND_MISSING_RT;

    // Case: compounds expected RT range is known
    if (compound->expectedRtMin > 0 && compound->expectedRtMax > 0) {
        if (debug) cout
                << compound->name
                 << " RT Range: ["
                << compound->expectedRtMin
                << " "
                << compound->expectedRt
                << " "
                << compound->expectedRtMax
                << "] vs. "
                << groupRtVal
                << endl;
        return groupRtVal >= compound->expectedRtMin && groupRtVal <= compound->expectedRtMax ? RtAgreementState::RT_AGREEMENT : RtAgreementState::RT_DISAGREEMENT;
    } else {
        //Case: compounds expected Rt range is not provided
        if (debug) cout
                << compound->name
                << " RT : "
                << compound->expectedRt
                << " vs. "
                << groupRtVal
                << endl;
        float rtDiff = abs(groupRtVal - compound->expectedRt);
        return rtDiff <= ms1RtTolr ? RtAgreementState::RT_AGREEMENT : RtAgreementState::RT_DISAGREEMENT;
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

            cout << "COMPOUNDS:\n";
            for (auto compound : compounds) {
                cout << compound->id << ": precMz=" << compound->precursorMz << "\n";
            }
            cout << endl;
        }

        auto compoundIons = vector<CompoundIon>(compounds.size());
        for (unsigned int i = 0; i < compounds.size(); i++) {
            compoundIons[i] = CompoundIon(compounds[i]);
        }

        for (auto& g : groups) {
            assignBestMetaboliteToGroup(
                        &g,
                        compoundIons,
                        params,
                        debug);
        }

}

void MzKitchenProcessor::assignBestMetaboliteToGroup(
        PeakGroup* g,
        vector<CompoundIon>& compounds,
        shared_ptr<MzkitchenMetaboliteSearchParameters> params,
        bool debug){

    if (!g) return;

    float peakGroupMz = g->meanMz;
    float peakGroupRt = g->meanRt;

    if (params->isUseGroupMaxPeakVals) {
        peakGroupMz = g->maxPeakMzVal;
        peakGroupRt = g->maxPeakRtVal;
    }

    float minMz = peakGroupMz - (peakGroupMz*params->ms1PpmTolr/1000000);
    float maxMz = peakGroupMz + (peakGroupMz*params->ms1PpmTolr/1000000);

    if (debug) {
        cout << "[minMz, maxMz] = [" << minMz << ", " << maxMz << "]" << endl;
    }

    auto lb = lower_bound(compounds.begin(), compounds.end(), minMz, [](const CompoundIon& lhs, const float& rhs){
        return lhs.precursorMz < rhs;
    });

    if (g->fragmentationPattern.mzs.empty()) {
        g->computeFragPattern(params.get(), debug);
    }

    vector<tuple<CompoundIon, FragmentationMatchScore, RtAgreementState>> scores{};

    for (long pos = lb - compounds.begin(); pos < static_cast<long>(compounds.size()); pos++){

        if (debug) cout << "pos=" << pos << endl;

        CompoundIon ion = compounds[static_cast<unsigned long>(pos)];
        Compound* compound = ion.compound;
        Adduct *adduct = ion.adduct;

        if (!compound) continue;

        if (debug) cout << "compound: " << compound->name << endl;

        bool isMatchingAdduct = adduct && compound->adductString == adduct->name;
        if (params->IDisRequireMatchingAdduct && !isMatchingAdduct) {
            continue;
        }

        float precMz = ion.precursorMz;

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

        //Issue 792: Altered logic around RT Agreement
        //Issue 816: Expanded options around RT matching to more properly deal with code paths.
        RtAgreementState rtAgreementState = MzKitchenProcessor::assessRtAgreement(peakGroupRt, compound, params->rtMatchTolerance, debug);

        // matches are only skipped if there is an explicit disagreement between compound and measured RT.
        // If insufficient information is available to make this comparison, the compound passes.
        if (params->rtIsRequireRtMatch && (rtAgreementState == RtAgreementState::RT_DISAGREEMENT)){
            continue;
        }

        FragmentationMatchScore s = library.scoreMatch(&(g->fragmentationPattern), params->ms2PpmTolr);

        double normCosineScore = Fragment::normCosineScore(&library, &(g->fragmentationPattern), s.ranks);
        s.dotProduct = normCosineScore;

        //debugging
        if (debug) {
            cout << "Peak group fragmentation pattern has " << g->fragmentationPattern.nobs() << " peaks.\n";
            cout << "Candidate Score: " << compound->id << ":\n";
            cout << "numMatches= " << s.numMatches
                 << ", hyperGeometricScore= " << s.hypergeomScore
                 << ", cosineScore= " << s.dotProduct
                 << " VS params->ms2MinNumMatches = " << params->ms2MinNumMatches
                 << " " << (s.numMatches >= params->ms2MinNumMatches ? "yes" : "no")
                 << " VS params->ms2MinScore = " << params->ms2MinScore
                 << " " << (s.dotProduct >= params->ms2MinScore ? "yes" : "no")
                 << "\n\n\n";
        }

        if (s.numMatches < params->ms2MinNumMatches) continue;

        //Issue 752: valid scores must meet score threshold
        if (s.dotProduct < params->ms2MinScore) continue;


        scores.push_back(make_tuple(ion, s, rtAgreementState));
    }

    //Issue 559: print empty groups
    if (debug) {
        cout << peakGroupMz << "@" << peakGroupRt << " #ms2s=" << g->ms2EventCount << " :" << endl;
    }

    if (scores.empty()) return;

    g->compounds = scores;

    float maxScore = -1.0f;
    tuple<CompoundIon, FragmentationMatchScore, RtAgreementState> bestScore;
    for (auto score : scores) {
        float dotProduct = get<1>(score).dotProduct;
        if (get<1>(score).dotProduct > maxScore) {
            bestScore = score;
            maxScore = dotProduct;
        }
    }

    if (maxScore > -1.0f) {

        Compound *compound = get<0>(bestScore).compound;
        Adduct *adduct =get<0>(bestScore).adduct;
        FragmentationMatchScore fragMatchScore = get<1>(bestScore);
        RtAgreementState rtAgreementState = get<2>(bestScore);

        g->compound = compound;

        if (adduct) g->adduct = adduct;
        g->fragMatchScore = fragMatchScore;
        g->fragMatchScore.mergedScore = fragMatchScore.dotProduct;

        if (compound->expectedRt > 0) {
            g->expectedRtDiff = abs(g->compound->expectedRt - peakGroupRt);
        }

        // Issue 816: Add an explicit label for RT agreement.
        // compounds missing RT values do not receive the label
        // Issue 828: label is only added once group is finalized, to prevent bug
        // where RT label might be added from a candidate that is ultimately filtered out
        // or usurped by a higher-scoring candidate.
        if (rtAgreementState == RtAgreementState::RT_AGREEMENT) {
            g->addLabel('l');
        }

        //Issue 835
        if (g->compound->metaDataMap.find(MzKitchenProcessor::MZKITCHEN_NOTES_KEY) != g->compound->metaDataMap.end()) {
            g->notes = g->compound->metaDataMap.at(MzKitchenProcessor::MZKITCHEN_NOTES_KEY);
        }
    }

    //Issue 546: debugging
    if (debug) {
        for (auto score : g->compounds) {
            FragmentationMatchScore fragmentationMatchScore = get<1>(score);
            cout << "\t" << get<0>(score).compound->name << " "
                 << fragmentationMatchScore.numMatches << " "
                 << fragmentationMatchScore.fractionMatched << " "
                 << fragmentationMatchScore.dotProduct  << " "
                 << fragmentationMatchScore.hypergeomScore << " "
                 << endl;
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

    string peakPickingEncodedParams = peakPickingAndGroupingParameters->getEncodedPeakParameters();

    encodedParams = encodedParams + peakPickingEncodedParams;

    //Issue 586
    encodedParams = encodedParams + getEncodedLipidParameters(TUPLE_MAP_KEY_DELIMITER, INTERNAL_MAP_DELIMITER);

    //Issue 792
    encodedParams = encodedParams + "isUseGroupMaxPeakVals" + "=" + to_string(isUseGroupMaxPeakVals) + ";";

    return encodedParams;
}

shared_ptr<LCLipidSearchParameters> LCLipidSearchParameters::decode(string encodedParams){
    shared_ptr<LCLipidSearchParameters> lipidSearchParameters = shared_ptr<LCLipidSearchParameters>(new LCLipidSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    lipidSearchParameters->fillInBaseParams(decodedMap);

    lipidSearchParameters->peakPickingAndGroupingParameters = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());
    lipidSearchParameters->peakPickingAndGroupingParameters->fillInPeakParameters(decodedMap);

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

    //Isuse 792
    if (decodedMap.find("isUseGroupMaxPeakVals") != decodedMap.end()) {
        lipidSearchParameters->isUseGroupMaxPeakVals = decodedMap["isUseGroupMaxPeakVals"] == "1";
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

    //Issue 792
    encodedParams = encodedParams + "isUseGroupMaxPeakVals" + "=" + to_string(isUseGroupMaxPeakVals) + ";";

    string peakPickingEncodedParams = peakPickingAndGroupingParameters->getEncodedPeakParameters();

    encodedParams = encodedParams + peakPickingEncodedParams;

    return encodedParams;
}

shared_ptr<MzkitchenMetaboliteSearchParameters> MzkitchenMetaboliteSearchParameters::decode(string encodedParams){

    shared_ptr<MzkitchenMetaboliteSearchParameters> metaboliteSearchParameters = shared_ptr<MzkitchenMetaboliteSearchParameters>(new MzkitchenMetaboliteSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    metaboliteSearchParameters->fillInBaseParams(decodedMap);

    metaboliteSearchParameters->peakPickingAndGroupingParameters = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());
    metaboliteSearchParameters->peakPickingAndGroupingParameters->fillInPeakParameters(decodedMap);

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

    //Isue 792
    if (decodedMap.find("isUseGroupMaxPeakVals") != decodedMap.end()) {
        metaboliteSearchParameters->isUseGroupMaxPeakVals = decodedMap["isUseGroupMaxPeakVals"] == "1";
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

// new parameters for encoding/decoding

string HRMSQCSearchParameters::encodeParams() {

    string encodedParams = baseParams();

    encodedParams = encodedParams + "ms1PpmTolr" + "=" + to_string(this->ms1PpmTolr) + ";";
    encodedParams = encodedParams + "rtTol" + "=" + to_string(this->rtTol) + ";";

    string peakPickingEncodedParams = peakPickingAndGroupingParameters->getEncodedPeakParameters();

    encodedParams = encodedParams + peakPickingEncodedParams;

    return encodedParams;
}

shared_ptr<HRMSQCSearchParameters> HRMSQCSearchParameters::decode(string encodedParams){

    shared_ptr<HRMSQCSearchParameters> hrmsQcSearchParameters = shared_ptr<HRMSQCSearchParameters>(new HRMSQCSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    hrmsQcSearchParameters->fillInBaseParams(decodedMap);

    hrmsQcSearchParameters->peakPickingAndGroupingParameters = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());
    hrmsQcSearchParameters->peakPickingAndGroupingParameters->fillInPeakParameters(decodedMap);

    if (decodedMap.find("ms1PpmTolr") != decodedMap.end()) {
        hrmsQcSearchParameters->ms1PpmTolr = stof(decodedMap["ms1PpmTolr"]);
    }
    if (decodedMap.find("rtTol") != decodedMap.end()) {
        hrmsQcSearchParameters->rtTol = stof(decodedMap["rtTol"]);
    }

    return hrmsQcSearchParameters;
}
