#include "mzSample.h"
#include <regex>


vector<SRMTransition*> QQQProcessor::getSRMTransitions(
        vector<mzSample*>& samples,
        shared_ptr<QQQSearchParameters> params,
        vector<Compound*>& compounds,
        vector<Adduct*>& adducts,
        bool debug) {

    //+118.001@cid34.00 [57.500-58.500]
    //+ c ESI SRM ms2 102.000@cid19.00 [57.500-58.500]
    //-87.000 [42.500-43.500]
    //- c ESI SRM ms2 159.000 [113.500-114.500]

    //regexes for precursorMz, productMz
    regex rx1a("[+/-](\\d+\\.\\d+)");
    regex rx1b("ms2\\s*(\\d+\\.\\d+)");
    regex rx2("(\\d+\\.\\d+)-(\\d+\\.\\d+)");

    //regex for transition name (Sciex 7500)
    regex rx3("name=.*$");

    //regex for collision energy (Sciex 7500)
    regex rx4("ce=[^\\s|$]+");

    float amuQ1 = params->amuQ1;
    float amuQ3 = params->amuQ3;

    //Issue 563: SRM transitions now include a transition ID string to distinguish data
    map<tuple<long, long, string>, SRMTransition*> srmTransitions{};

    for(unsigned int i=0; i < samples.size(); i++ ) {
        mzSample* sample = samples[i];

        //Issue 658: Within a sample, only need to detect a single SRM ID string.
        //However, this should not be used across samples, as samples in different
        //plates may have identical SRM ID strings.
        set<string> srms{};

        for(unsigned int j=0; j < sample->scans.size(); j++ ) {
            Scan* scan = sample->getScan(j);
            if (!scan) continue;

            // -------------------------- //
            // Extract SRM data from scan //
            // -------------------------- //

            string filterLine(scan->filterLine.c_str());
            if (filterLine.empty()) continue;

            if (srms.find(filterLine) != srms.end()) continue;
            srms.insert(filterLine);

            float precursorMz = scan->precursorMz;
            float productMz   = scan->productMz;
            string transitionName = "";
            float collisionEnergy = 0;

            int  polarity= scan->getPolarity();
            if (polarity==0) filterLine[0] == '+' ? polarity=1 : polarity =-1;

            if ( precursorMz < 1e-6f ) {

                smatch m1, m2;
                regex_search(filterLine, m1, rx1a);
                regex_search(filterLine, m2, rx1b);

                if (m1.size() > 0) {
                    precursorMz = stof(m1[0]);
                } else if (m2.size() > 0){
                    precursorMz = stof(m2[0]);
                }
            }

            if ( productMz < 1e-6f ) {

                smatch m3;
                regex_search(filterLine, m3, rx2);

                if (m3.size() > 2) {
                    float lb = stof(m3[1]);
                    float ub = stof(m3[2]);
                    productMz = lb+(0.5f*(ub-lb));
                }
            }

            smatch m4;
            regex_search(filterLine, m4, rx3);
            if (!m4.empty()) {
                transitionName = string(m4[0]);
                transitionName = transitionName.substr(5);
            }

            smatch m5;
            regex_search(filterLine, m5, rx4);
            if (!m5.empty()) {
                string ce = string(m5[0]);
                ce = ce.substr(3);
                collisionEnergy = stof(ce);
            }

            // ------------- //
            // SRMTransition //
            // ------------- //

            // Issue 578: silently ignore any SRM scans where the precursorMz, productMz are invalid.
            if (precursorMz <= 0 || productMz <= 0) continue;

            tuple<long, long, string> srmKey = make_tuple(
                        mzUtils::mzToIntKey(static_cast<double>(precursorMz)),
                        mzUtils::mzToIntKey(static_cast<double>(productMz)),
                        transitionName);

            SRMTransition *srmTransition;
            if (srmTransitions.find(srmKey) != srmTransitions.end()) {
                srmTransition = srmTransitions[srmKey];
            } else {
                srmTransition = new SRMTransition();
                srmTransition->precursorMz = precursorMz;
                srmTransition->productMz = productMz;
                srmTransition->name = transitionName;
                srmTransition->collisionEnergy = collisionEnergy;
            }

            if (debug) {
                srmTransition->printKey();
                cout << ": " << endl;
            }

            //Associated compounds with SRMTransition, if applicable
            if (precursorMz > 0 && productMz > 0 ) {
                for (auto db_compound : compounds) {

                    //only consider SRM compounds
                    if (db_compound->precursorMz <= 0 || db_compound->productMz <= 0) continue;

                    //compounds must be within tolerance
                    auto q1_dist = abs(db_compound->precursorMz-precursorMz) ;
                    if (q1_dist > amuQ1) continue;

                    auto q3_dist = abs(db_compound->productMz-productMz);
                    if (q3_dist > amuQ3) continue;

                    // Issue 563/564: match transition_id filter
                    // If a transition ID is known, it must be present in the srmId, for this compound to be correct.
                    // If no transitionId is provided, match only on precursorMz and productMz.
                    // This allows that the compound name and transition name be different strings, though by default
                    // the transition_id will usually be the same as the compound name.
                    if (!db_compound->srmId.empty()) {

                        string compoundTransitionId = db_compound->srmId;

                        //compounds should match, except for whitespace differences, and removing special characters
                        string compoundTransitionIdForComparison = compoundTransitionId;
                        string transitionNameForComparison = transitionName;

                        //remove certain tricky characters, white space
                        mzUtils::replaceAll(transitionNameForComparison, "(±)", "");
                        mzUtils::replaceAll(transitionNameForComparison, " ", "");

                        mzUtils::replaceAll(compoundTransitionIdForComparison, "(±)", "");
                        mzUtils::replaceAll(compoundTransitionIdForComparison, " ", "");

                        if (compoundTransitionIdForComparison != transitionNameForComparison) continue;
                    }

                    Adduct* db_adduct = nullptr;
                    if (!db_compound->adductString.empty()) {
                        for (auto availableAdduct : adducts) {
                            if (availableAdduct->name == db_compound->adductString) {
                                db_adduct = availableAdduct;
                                break;
                            }
                        }
                    }

                    if (debug) {
                        cout << "    " << db_compound->id << endl;
                    }

                    srmTransition->addCompound(db_compound, db_adduct);
                }
            }

            if (srmTransition->srmIdBySample.find(sample) == srmTransition->srmIdBySample.end()) {
                srmTransition->srmIdBySample.insert(make_pair(sample, set<string>{}));
            }
            srmTransition->srmIdBySample.at(sample).insert(filterLine);
            srmTransition->srmIds.insert(filterLine);

            srmTransitions[srmKey] = srmTransition;
        }
    }

    vector<SRMTransition*> transitionsZeroCompounds{};
    vector<SRMTransition*> transitionsOneCompound{};
    vector<SRMTransition*> transitionsMultipleCompounds{};

    for (auto it = srmTransitions.begin(); it != srmTransitions.end(); ++it) {
        unsigned long numCompounds = it->second->compounds.size();
        if (numCompounds == 0) {
            transitionsZeroCompounds.push_back(it->second);
        } else if (numCompounds == 1) {
            transitionsOneCompound.push_back(it->second);
        } else {
            transitionsMultipleCompounds.push_back(it->second);
        }
    }

    vector<SRMTransition*> srmTransitionsAsVector = transitionsOneCompound;

    if (params->transitionCompoundMappingPolicy == REQUIRE_ALL_TRANSITIONS_EXACTLY_ONE_COMPOUND) {
        if (!transitionsZeroCompounds.empty() || !transitionsMultipleCompounds.empty()) {
            cout << "=====================================\n";
            cout << "ERROR: The current transition to compound mapping policy is set to require that all transitions map to exactly one compound." << endl;

            cout <<  transitionsZeroCompounds.size() << " transitions did not map to any compounds:" << endl;
            for (auto transition : transitionsZeroCompounds) {
                cout << transition->getKey() << endl;
            }
            cout << endl;

            cout << transitionsMultipleCompounds.size() << " transitions mapped to more than one compound:" << endl;

            for (auto transition : transitionsMultipleCompounds) {
                cout << transition->getKey() << endl;
                for (auto p : transition->compounds) {
                    cout << "   " << p.first->id << endl;
                }

            }

            cout << "=====================================\n";
            cout << endl;
            exit(-1);
        }
    } else if (params->transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::REQUIRE_ALL_TRANSITIONS_TO_ONE_OR_MORE_COMPOUNDS) {
        if (!transitionsZeroCompounds.empty()) {
            cout << "=====================================\n";
            cout << "ERROR: The current transition to compound mapping policy is set to require that all transitions map to at least one compound." << endl;

            cout <<  transitionsZeroCompounds.size() << " transitions did not map to any compounds:" << endl;
            for (auto transition : transitionsZeroCompounds) {
                cout << transition->getKey() << endl;
            }

            cout << "=====================================\n";
            cout << endl;
            exit(-1);
        } else {
            srmTransitionsAsVector.insert(srmTransitionsAsVector.end(), transitionsMultipleCompounds.begin(), transitionsMultipleCompounds.end());
        }
    } else if (params->transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_ONE_OR_MORE_COMPOUNDS) {
        srmTransitionsAsVector.insert(srmTransitionsAsVector.end(), transitionsMultipleCompounds.begin(), transitionsMultipleCompounds.end());
    } else if (params->transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_EXACTLY_ONE_COMPOUND) {
        //done during initialization, leaving here for readability
        //srmTransitionsAsVector = transitionsOneCompound;
    } else { //if an unfamiliar policy is detected, same as RETAIN_ALL_TRANSITIONS
        srmTransitionsAsVector.insert(srmTransitionsAsVector.end(), transitionsMultipleCompounds.begin(), transitionsMultipleCompounds.end());
        srmTransitionsAsVector.insert(srmTransitionsAsVector.end(), transitionsZeroCompounds.begin(), transitionsZeroCompounds.end());
    }

    sort(srmTransitionsAsVector.begin(), srmTransitionsAsVector.end(), [](SRMTransition* lhs, SRMTransition* rhs){
       return lhs->getKey() < rhs->getKey();
    });

    if (debug) {
        cout << "# SRMTransitions: " << srmTransitionsAsVector.size() << endl;
    }

    return srmTransitionsAsVector;
}

string QQQSearchParameters::encodeParams(){

    string encodedParams = baseParams();

    encodedParams = encodedParams + "amuQ1" + "=" + to_string(amuQ1) + ";";
    encodedParams = encodedParams + "amuQ3" + "=" + to_string(amuQ3) + ";";
    encodedParams = encodedParams + "transitionListFilePath" + "=" + transitionListFilePath + ";";

    encodedParams = encodedParams + "transitionCompoundMappingPolicy" + "=";
    if (transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::REQUIRE_ALL_TRANSITIONS_EXACTLY_ONE_COMPOUND) {
        encodedParams = encodedParams + "REQUIRE_ALL_TRANSITIONS_EXACTLY_ONE_COMPOUND";
    } else if (transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::REQUIRE_ALL_TRANSITIONS_TO_ONE_OR_MORE_COMPOUNDS) {
        encodedParams = encodedParams + "REQUIRE_ALL_TRANSITIONS_TO_ONE_OR_MORE_COMPOUNDS";
    } else if (transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_ONE_OR_MORE_COMPOUNDS) {
        encodedParams = encodedParams + "RETAIN_TRANSITIONS_ONE_OR_MORE_COMPOUNDS";
    } else if (transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_EXACTLY_ONE_COMPOUND) {
        encodedParams = encodedParams + "RETAIN_TRANSITIONS_EXACTLY_ONE_COMPOUND";
    } else if (transitionCompoundMappingPolicy == QQQTransitionCompoundMappingPolicy::RETAIN_ALL_TRANSITIONS) {
        encodedParams = encodedParams + "RETAIN_ALL_TRANSITIONS";
    } else {
        encodedParams = encodedParams + "UNKNOWN";
    }
    encodedParams = encodedParams + ";";

    encodedParams = encodedParams + "rollUpRtTolerance" + "=" + to_string(rollUpRtTolerance) + ";";

    encodedParams = encodedParams + "qqqFilterMinSignalBlankRatio" + "=" + to_string(qqqFilterMinSignalBlankRatio) + ";";
    encodedParams = encodedParams + "qqqFilterMinPeakIntensityGroupBackgroundRatio" + "=" + to_string(qqqFilterMinPeakIntensityGroupBackgroundRatio) + ";";
    encodedParams = encodedParams + "qqqFilterIsRetainOnlyPassingPeaks" + "=" + to_string(qqqFilterIsRetainOnlyPassingPeaks) + ";";

    string peakPickingEncodedParams = peakPickingAndGroupingParameters->getEncodedPeakParameters();

    encodedParams = encodedParams + peakPickingEncodedParams;

    return encodedParams;
}

shared_ptr<QQQSearchParameters> QQQSearchParameters::decode(string encodedParams) {

    shared_ptr<QQQSearchParameters> qqqSearchParameters = shared_ptr<QQQSearchParameters>(new QQQSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    qqqSearchParameters->fillInBaseParams(decodedMap);

    qqqSearchParameters->peakPickingAndGroupingParameters = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());
    qqqSearchParameters->peakPickingAndGroupingParameters->fillInPeakParameters(decodedMap);

    // START QQQSearchParameters

    if (decodedMap.find("amuQ1") != decodedMap.end()) {
        qqqSearchParameters->amuQ1 = stof(decodedMap["amuQ1"]);
    }
    if (decodedMap.find("amuQ3") != decodedMap.end()) {
        qqqSearchParameters->amuQ3 = stof(decodedMap["amuQ3"]);
    }
    if (decodedMap.find("transitionListFilePath") != decodedMap.end()) {
        qqqSearchParameters->transitionListFilePath = decodedMap["transitionListFilePath"];
    }
    if (decodedMap.find("transitionCompoundMappingPolicy") != decodedMap.end()) {
        string transitionCompoundMappingPolicyStr = decodedMap["transitionCompoundMappingPolicy"];
        if (transitionCompoundMappingPolicyStr == "REQUIRE_ALL_TRANSITIONS_EXACTLY_ONE_COMPOUND") {
            qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::REQUIRE_ALL_TRANSITIONS_EXACTLY_ONE_COMPOUND;
        } else if (transitionCompoundMappingPolicyStr == "REQUIRE_ALL_TRANSITIONS_TO_ONE_OR_MORE_COMPOUNDS") {
            qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::REQUIRE_ALL_TRANSITIONS_TO_ONE_OR_MORE_COMPOUNDS;
        } else if (transitionCompoundMappingPolicyStr == "RETAIN_TRANSITIONS_ONE_OR_MORE_COMPOUNDS") {
            qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_ONE_OR_MORE_COMPOUNDS;
        } else if (transitionCompoundMappingPolicyStr == "RETAIN_TRANSITIONS_EXACTLY_ONE_COMPOUND") {
            qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_EXACTLY_ONE_COMPOUND;
        } else if (transitionCompoundMappingPolicyStr == "RETAIN_ALL_TRANSITIONS") {
            qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::RETAIN_ALL_TRANSITIONS;
        }
    }
    if (decodedMap.find("rollUpRtTolerance") != decodedMap.end()) {
        qqqSearchParameters->rollUpRtTolerance = stof(decodedMap["rollUpRtTolerance"]);
    }
    if (decodedMap.find("qqqFilterMinSignalBlankRatio") != decodedMap.end()) {
        qqqSearchParameters->qqqFilterMinSignalBlankRatio = stof(decodedMap["qqqFilterMinSignalBlankRatio"]);
    }
    if (decodedMap.find("qqqFilterMinPeakIntensityGroupBackgroundRatio") != decodedMap.end()) {
        qqqSearchParameters->qqqFilterMinPeakIntensityGroupBackgroundRatio = stof(decodedMap["qqqFilterMinPeakIntensityGroupBackgroundRatio"]);
    }
    if (decodedMap.find("qqqFilterIsRetainOnlyPassingPeaks") != decodedMap.end()) {
        qqqSearchParameters->qqqFilterIsRetainOnlyPassingPeaks = decodedMap["qqqFilterIsRetainOnlyPassingPeaks"] == "1";
    }

    // END QQQSearchParameters

    return qqqSearchParameters;
}


//Issue 568: Create dedicated method to convert SRMTransitions to mzSlice
vector<mzSlice*> QQQProcessor::getMzSlices(
        vector<SRMTransition*>& transitions,
        bool isRequireCompound,
        bool debug) {

    vector<mzSlice*> slices{};

    for (SRMTransition *transition : transitions) {

        if (isRequireCompound && !transition->compound) {
            // pass
        } else {
            mzSlice *slice = new mzSlice(transition);
            slices.push_back(slice);

            if (debug) {
                transition->printKey();
                cout << ":\n";
            }

            sort(slice->compoundVector.begin(), slice->compoundVector.end(), [](Compound* lhs, Compound* rhs){
                return lhs->id < rhs->id;
            });

            if (debug) {
                for (auto p : slice->compoundVector) {
                    cout << "   " << p->id << endl;
                }
            }
        }
    }

    return slices;
}

/**
 * @brief QQQProcessor::printSRMIds
 * @param transitions
 *
 * Print all srmIds, limiting to one (Q1, Q3, transition_name)
 * (print some random sample id number with the srmId information)
 */
void QQQProcessor::printSRMIds(vector<SRMTransition*>& transitions) {
    for (SRMTransition *transition : transitions) {
        if (transition && !transition->srmIds.empty()) {
            string srmIdString = *(transition->srmIds.begin());
            cout << srmIdString << endl;
        }
    }
}

/**
 * @brief QQQProcessor::getSRMIds
 * @param transitions
 * @return
 *
 * Retrieve all srmId information, from every sample.
 */
set<string> QQQProcessor::getSRMIds(vector<SRMTransition*>& transitions) {
    set<string> srmIds{};
    for (SRMTransition *transition : transitions) {
        if (transition && !transition->srmIds.empty()) {
            for (auto srmId : transition->srmIds) {
                srmIds.insert(srmId);
            }
        }
    }
    return srmIds;
}

void QQQProcessor::rollUpToCompoundQuant(vector<PeakGroup>& peakgroups, shared_ptr<QQQSearchParameters> params, bool debug){
    if (debug) cout << "Start QQQProcessor::rollUpToCompoundQuant()" << endl;

    vector<PeakGroup*> references(peakgroups.size());
    transform(peakgroups.begin(), peakgroups.end(), references.begin(), [](PeakGroup& pg){return &pg;});

    map<string, vector<PeakGroup*>> groupsByCategory{};
    map<string, vector<PeakGroup*>> allGroupsByCategory{};

    for (auto pg : references) {
        if (pg->compound && !pg->compound->category.empty()) {

            //Issue 635: All peakgroups should transfer over the peakRank quant, for possible manual reassignments.
            string quantType = "smoothedPeakAreaCorrected";
            if (pg->compound->metaDataMap.find(QQQProcessor::getTransitionPreferredQuantTypeStringKey()) != pg->compound->metaDataMap.end()) {
                quantType = pg->compound->metaDataMap.at(QQQProcessor::getTransitionPreferredQuantTypeStringKey());
            }

            float maxPeakRank = 0.0f;
            for (auto & p : pg->peaks) {
                p.peakRank = p.getQuantByName(quantType);
                if (p.peakRank > maxPeakRank) {
                    maxPeakRank = p.peakRank;
                }
            }

            pg->groupRank = maxPeakRank;

            string category = pg->compound->category.at(0);

            if (allGroupsByCategory.find(category) == allGroupsByCategory.end()) {
                allGroupsByCategory.insert(make_pair(category, vector<PeakGroup*>{}));
            }
            allGroupsByCategory[category].push_back(pg);

            //compounds are used for quant only if they are explicitly designated as such
            if (pg->compound->metaDataMap.find(QQQProcessor::getTransitionIonTypeFilterStringKey()) != pg->compound->metaDataMap.end()) {
                string quantType = pg->compound->metaDataMap.at(QQQProcessor::getTransitionIonTypeFilterStringKey());
                if (quantType == "qq" || quantType == "quant") {
                    if (groupsByCategory.find(category) == groupsByCategory.end()) {
                        groupsByCategory.insert(make_pair(category, vector<PeakGroup*>{}));
                    }
                    groupsByCategory[category].push_back(pg);
                }
            }
        }
    }

    //Take the PG with the highest maxSmoothedIntensity that is within tolerance

    for (auto it = groupsByCategory.begin(); it != groupsByCategory.end(); ++it) {
        string category = it->first;
        vector<PeakGroup*> peakGroups = it->second;

        bool isCheckRt = params->rollUpRtTolerance > 0 && peakGroups.at(0)->compound->expectedRt > 0;

        Compound* representativeCompound = nullptr;
        PeakGroup* representative = nullptr;
        for (auto pg : peakGroups) {
           if (isCheckRt && abs(pg->maxPeakRt() - pg->compound->expectedRt) > params->rollUpRtTolerance) {
                continue;
           }

           if (!representative || pg->maxIntensity > representative->maxIntensity) {
               representative = pg;
               representativeCompound = pg->compound;
           }

        }

        if (representative) {
            representative->addLabel('q');

            //Issue 635: Organize non-representative peakgroups as children of representatives.
            if (allGroupsByCategory.find(category) != allGroupsByCategory.end()) {
                vector<PeakGroup*> allPeakGroups = allGroupsByCategory[category];
                for (auto pg : allPeakGroups){
                    if (pg != representative) {
                        representative->addChild(*pg);
                        pg->deletedFlag = true;
                    }
                }
            }

        }
    }

    if (debug) cout << "End QQQProcessor::rollUpToCompoundQuant()" << endl;
}

void QQQProcessor::labelInternalStandards(vector<PeakGroup>& peakgroups, shared_ptr<QQQSearchParameters> params, bool debug){
    if (debug) cout << "Start QQQProcessor::labelInternalStandards()" << endl;

    for (auto & pg : peakgroups) {
        if (pg.compound && pg.compound->metaDataMap.find(QQQProcessor::getTransitionIsInternalStandardStringKey()) != pg.compound->metaDataMap.end()){
            if (pg.compound->metaDataMap.at(QQQProcessor::getTransitionIsInternalStandardStringKey()) == "TRUE") {
                pg.addLabel('i');
                if (debug) cout << pg.compound->name << " is an internal standard." << endl;

                bool isBothRepresentativeAndIS = false;
                for (auto label : pg.labels) {
                    if (label == 'q') {
                        isBothRepresentativeAndIS = true;
                        break;
                    }
                }

                if (isBothRepresentativeAndIS) {
                    pg.addLabel('l');
                }
            }
        }

    }

    if (debug) cout << "End QQQProcessor::labelInternalStandards()" << endl;
}

bool QQQProcessor::isPassPeakGroupBackground(
    Peak& p,
    PeakGroup& pg,
    shared_ptr<QQQSearchParameters> params,
    float peakQuant,
    bool debug
    ) {

    if (params->peakPickingAndGroupingParameters->groupBackgroundType == PeakGroupBackgroundType::MAX_BLANK_INTENSITY) {
        return p.peakIntensity >= pg.blankMaxHeight * params->qqqFilterMinPeakIntensityGroupBackgroundRatio;
    } else if (params->peakPickingAndGroupingParameters->groupBackgroundType == PeakGroupBackgroundType::MEDIAN_BLANK_INTENSITY) {
        return p.peakIntensity >= pg.blankMedianHeight * params->qqqFilterMinPeakIntensityGroupBackgroundRatio;
    } else if (params->peakPickingAndGroupingParameters->groupBackgroundType == PeakGroupBackgroundType::PREFERRED_QUANT_TYPE_MERGED_EIC_BASELINE ||
               params->peakPickingAndGroupingParameters->groupBackgroundType == PeakGroupBackgroundType::PREFERRED_QUANT_TYPE_MAX_BLANK_SIGNAL) {
        return peakQuant >= pg.groupBackground * params->qqqFilterMinPeakIntensityGroupBackgroundRatio;
    }

    return true;
}

//DOES NOT apply to children peak groups.
vector<PeakGroup> QQQProcessor::filterPeakGroups(vector<PeakGroup>& peakgroups, shared_ptr<QQQSearchParameters> params, bool debug){
    if (debug) cout << "Start QQQProcessor::labelInternalStandards()" << endl;
    vector<PeakGroup> filteredGroups{};

    for (auto & pg : peakgroups) {
        if (pg.compound && pg.compound->metaDataMap.find(QQQProcessor::getTransitionIsInternalStandardStringKey()) != pg.compound->metaDataMap.end()){
            if (pg.compound->metaDataMap.at(QQQProcessor::getTransitionIsInternalStandardStringKey()) == "TRUE") {
                filteredGroups.push_back(pg);
                if (debug) cout << "( " << pg.meanMz << ", " << pg.medianRt() << "): " << pg.compound->name << " is IS" << endl;

                continue;
            }
        }

        string quantType = "smoothedPeakAreaCorrected";
        if (pg.compound && pg.compound->metaDataMap.find(QQQProcessor::getTransitionPreferredQuantTypeStringKey()) != pg.compound->metaDataMap.end()) {
            quantType = pg.compound->metaDataMap.at(QQQProcessor::getTransitionPreferredQuantTypeStringKey());
        }

        float maxBlankQuant = 0.0f;
        float maxSampleQuant = 0.0f;
        float maxPeakIntensity = 0.0f;
        Peak maxPeak;

        for (auto & p : pg.peaks) {

            float peakQuant = p.getQuantByName(quantType);

            if (p.fromBlankSample && peakQuant > maxBlankQuant) {
                maxBlankQuant = peakQuant;
            } else if (!p.fromBlankSample) {
                if (peakQuant > maxSampleQuant) {
                    maxSampleQuant = peakQuant;
                }
                if (p.peakIntensity > maxPeakIntensity) {
                    maxPeakIntensity = p.peakIntensity;
                    maxPeak = p;
                }
            }
            if (debug) cout << p.sample->sampleName.c_str() << ": (blank=" << (p.fromBlankSample ? "yes" : "no") << ") " << peakQuant << endl;
        }

        float maxPeakQuant = maxPeak.getQuantByName(quantType);

        if (debug && pg.compound) cout << pg.compound->name << " ";
        if (debug) cout << "( " << pg.meanMz << ", " << pg.medianRt() << "): "
                        << "quantType = " << quantType << ", "
                        << "(maxSampleQuant/maxBlankQuant) = (" << maxSampleQuant << "/" << maxBlankQuant << ") = "
                        << (maxSampleQuant/maxBlankQuant)  << " "
                        << (maxBlankQuant * params->qqqFilterMinSignalBlankRatio <= maxSampleQuant ? "keep" : "skip")
                        << endl;

        //Issue 664
        if (params->qqqFilterIsRetainOnlyPassingPeaks) {

            bool isKeepPeakGroup = false;

            vector<Peak> passingPeaks{};

            for (auto p : pg.peaks) {
                float peakQuant = p.getQuantByName(quantType);
                if (p.fromBlankSample ||
                        (
                            peakQuant >= maxBlankQuant * params->qqqFilterMinSignalBlankRatio &&
                            QQQProcessor::isPassPeakGroupBackground(p, pg, params, peakQuant, debug)
                        )
                    ) {
                    passingPeaks.push_back(p);

                    if (!p.fromBlankSample) isKeepPeakGroup = true;
                }
            }

            if (isKeepPeakGroup) {

                pg.peaks.clear();
                for (auto& p : passingPeaks) {
                    pg.addPeak(p);
                }

                pg.groupStatistics(true); // force recomputation

                filteredGroups.push_back(pg);
            }

        } else if (
           maxSampleQuant >= maxBlankQuant * params->qqqFilterMinSignalBlankRatio &&
            QQQProcessor::isPassPeakGroupBackground(maxPeak, pg, params, maxPeakQuant, debug) ){
            filteredGroups.push_back(pg);
        }
    }

    if (debug) cout << "End QQQProcessor::filterPeakGroups()" << endl;
    return filteredGroups;
}

void QQQProcessor::setPeakGroupBackground(
    vector<PeakGroup>& peakgroups,
    shared_ptr<QQQSearchParameters> params,
    bool debug
    ) {

    if (debug) {
        cout << "QQQProcessor::setPeakGroupBackground(): "
             << peakgroups.size() << " peakgroups."
             << endl;
    }

    for (auto & pg : peakgroups) {

        if (!pg.compound) continue;

        string quantType = "smoothedPeakAreaCorrected";
        if (pg.compound->metaDataMap.find(QQQProcessor::getTransitionPreferredQuantTypeStringKey()) != pg.compound->metaDataMap.end()) {
            quantType = pg.compound->metaDataMap.at(QQQProcessor::getTransitionPreferredQuantTypeStringKey());

            //Issue 677: Save the peak group background as the corresponding quant from the blank samples, unless
            //choosing the merged EIC signal.
            //Filtering done in QQQProcessor::isPassPeakGroupBackgroumd() always respects the appropriate type.
            //
            if (params->peakPickingAndGroupingParameters->groupBackgroundType == PeakGroupBackgroundType::PREFERRED_QUANT_TYPE_MERGED_EIC_BASELINE) {
                pg.groupBackground = pg.mergedEICSummaryData.getCorrespondingBaseline(quantType);
            } else {
                pg.groupBackground = pg.getMaxBlankCorrespondingQuant(quantType);
            }

            if (debug) {
                cout << pg.compound->name << "@ "<< pg.medianRt()
                     << " min: quantType=" << quantType
                     << ", groupBackground=" << pg.groupBackground
                     << endl;
            }
        }
    }
}
