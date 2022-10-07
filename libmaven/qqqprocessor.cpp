#include "mzSample.h"
#include <regex>

pair<vector<mzSlice*>, vector<SRMTransition*>> QQQProcessor::getSRMSlices(
        vector<mzSample*>& samples,
        shared_ptr<QQQSearchParameters> params,
        vector<Compound*>& compounds,
        vector<Adduct*>& adducts,
        bool debug){

    set<string> srms;
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
    map<tuple<float, float, string>, SRMTransition*> srmTransitions{};

    for(unsigned int i=0; i < samples.size(); i++ ) {
        mzSample* sample = samples[i];
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

            tuple<float, float, string> srmKey = make_tuple(precursorMz, productMz, transitionName);

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

    vector<mzSlice*> slices{};  //Build slices from completed srmTransitions
    vector<SRMTransition*> srmTransitionsAsVector(srmTransitions.size());

    unsigned long counter = 0;
    for (auto it = srmTransitions.begin(); it != srmTransitions.end(); ++it) {
        srmTransitionsAsVector[counter] = it->second;
        counter++;

        for (string srmId : it->second->srmIds) {
            for (auto p : it->second->compounds) {
                mzSlice* s = new mzSlice(0,0,0,0);
                s->srmId = srmId;
                s->compound = p.first;
                s->adduct = p.second;
                if (s->compound) {
                    s->rt = s->compound->expectedRt;
                }
                slices.push_back(s);
            }
            if (!it->second->compound) {
                mzSlice* s = new mzSlice(0,0,0,0);
                s->srmId = srmId;
                slices.push_back(s);
            }
        }
    }

    sort(srmTransitionsAsVector.begin(), srmTransitionsAsVector.end(), [](SRMTransition* lhs, SRMTransition* rhs){
       return lhs->getKey() < rhs->getKey();
    });

    if (debug) {
        cout << "# SRMTransitions: " << srmTransitions.size()
             << ", # SRMs: " << srms.size()
             << ", # Slices: " << slices.size();
    }

    return make_pair(slices, srmTransitionsAsVector);
}

string QQQSearchParameters::encodeParams(){

    string encodedParams = baseParams();

    encodedParams = encodedParams + "amuQ1" + "=" + to_string(amuQ1) + ";";
    encodedParams = encodedParams + "amuQ3" + "=" + to_string(amuQ3) + ";";
    encodedParams = encodedParams + "transitionListFilePath" + "=" + transitionListFilePath + ";";

    return encodedParams;
}

shared_ptr<QQQSearchParameters> QQQSearchParameters::decode(string encodedParams) {

    shared_ptr<QQQSearchParameters> qqqSearchParameters = shared_ptr<QQQSearchParameters>(new QQQSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    //START SearchParameters

    //program level
    if (decodedMap.find("searchVersion") != decodedMap.end()) {
        qqqSearchParameters->searchVersion = decodedMap["searchVersion"];
    }

    //scan filter params
    if (decodedMap.find("scanFilterMinFracIntensity") != decodedMap.end()){
        qqqSearchParameters->scanFilterMinFracIntensity = stof(decodedMap["scanFilterMinFracIntensity"]);
    }
    if (decodedMap.find("scanFilterMinSNRatio") != decodedMap.end()){
        qqqSearchParameters->scanFilterMinSNRatio = stof(decodedMap["scanFilterMinSNRatio"]);
    }
    if (decodedMap.find("scanFilterMaxNumberOfFragments") != decodedMap.end()) {
        qqqSearchParameters->scanFilterMaxNumberOfFragments = stoi(decodedMap["scanFilterMaxNumberOfFragments"]);
    }
    if (decodedMap.find("scanFilterBaseLinePercentile") != decodedMap.end()) {
        qqqSearchParameters->scanFilterBaseLinePercentile = stoi(decodedMap["scanFilterBaseLinePercentile"]);
    }
    if (decodedMap.find("scanFilterIsRetainFragmentsAbovePrecursorMz") != decodedMap.end()) {
        qqqSearchParameters->scanFilterIsRetainFragmentsAbovePrecursorMz = decodedMap["scanFilterIsRetainFragmentsAbovePrecursorMz"] == "1";
    }
    if (decodedMap.find("scanFilterPrecursorPurityPpm") != decodedMap.end()){
        qqqSearchParameters->scanFilterPrecursorPurityPpm = stof(decodedMap["scanFilterPrecursorPurityPpm"]);
    }
    if (decodedMap.find("scanFilterMinIntensity") != decodedMap.end()){
        qqqSearchParameters->scanFilterMinIntensity = stof(decodedMap["scanFilterMinIntensity"]);
    }

    //scan filter for MS1 scans
    if (decodedMap.find("scanFilterMs1MinRt") != decodedMap.end()) {
        qqqSearchParameters->scanFilterMs1MinRt = stof(decodedMap["scanFilterMs1MinRt"]);
    }
    if (decodedMap.find("scanFilterMs1MaxRt") != decodedMap.end()) {
        qqqSearchParameters->scanFilterMs1MaxRt = stof(decodedMap["scanFilterMs1MaxRt"]);
    }

    //scan filter for MS2 scans
    if (decodedMap.find("scanFilterMs2MinRt") != decodedMap.end()) {
        qqqSearchParameters->scanFilterMs2MinRt = stof(decodedMap["scanFilterMs2MinRt"]);
    }
    if (decodedMap.find("scanFilterMs2MaxRt") != decodedMap.end()) {
        qqqSearchParameters->scanFilterMs2MaxRt = stof(decodedMap["scanFilterMs2MaxRt"]);
    }

    //consensus spectrum params (all ms levels)

    if (decodedMap.find("consensusIsIntensityAvgByObserved") != decodedMap.end()){
        qqqSearchParameters->consensusIsIntensityAvgByObserved = decodedMap["consensusIsIntensityAvgByObserved"] == "1";
    }
    if (decodedMap.find("consensusIsNormalizeTo10K") != decodedMap.end()){
        qqqSearchParameters->consensusIsNormalizeTo10K = decodedMap["consensusIsNormalizeTo10K"] == "1";
    }
    if (decodedMap.find("consensusIntensityAgglomerationType") != decodedMap.end()) {
        string consensusIntensityAgglomerationTypeStr = decodedMap["consensusIntensityAgglomerationType"];
        if (consensusIntensityAgglomerationTypeStr == "MEAN") {
            qqqSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
        } else if (consensusIntensityAgglomerationTypeStr == "MEDIAN") {
            qqqSearchParameters->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
        }
    }

    //ms1 consensus spectrum params
    if (decodedMap.find("consensusMs1PpmTolr") != decodedMap.end()){
        qqqSearchParameters->consensusMs1PpmTolr = stof(decodedMap["consensusMs1PpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs1Scans") != decodedMap.end()){
        qqqSearchParameters->consensusMinNumMs1Scans = stoi(decodedMap["consensusMinNumMs1Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs1Scans") != decodedMap.end()){
        qqqSearchParameters->consensusMinFractionMs1Scans = stof(decodedMap["consensusMinFractionMs1Scans"]);
    }

    //ms2 consensus spectrum params
    if (decodedMap.find("consensusPpmTolr") != decodedMap.end()){
        qqqSearchParameters->consensusPpmTolr = stof(decodedMap["consensusPpmTolr"]);
    }
    if (decodedMap.find("consensusMinNumMs2Scans") != decodedMap.end()){
        qqqSearchParameters->consensusMinNumMs2Scans = stoi(decodedMap["consensusMinNumMs2Scans"]);
    }
    if (decodedMap.find("consensusMinFractionMs2Scans") != decodedMap.end()){
        qqqSearchParameters->consensusMinFractionMs2Scans = stof(decodedMap["consensusMinFractionMs2Scans"]);
    }

    // ms1 matching
    if (decodedMap.find("ms1PpmTolr") != decodedMap.end()) {
        qqqSearchParameters->ms1PpmTolr = stof(decodedMap["ms1PpmTolr"]);
    }

    // ms2 search
    if (decodedMap.find("ms2MinNumMatches") != decodedMap.end()) {
        qqqSearchParameters->ms2MinNumMatches = stoi(decodedMap["ms2MinNumMatches"]);
    }
    if (decodedMap.find("ms2MinNumDiagnosticMatches") != decodedMap.end()) {
        qqqSearchParameters->ms2MinNumDiagnosticMatches = stoi(decodedMap["ms2MinNumDiagnosticMatches"]);
    }
    if (decodedMap.find("ms2MinNumUniqueMatches") != decodedMap.end()) {
        qqqSearchParameters->ms2MinNumUniqueMatches = stoi(decodedMap["ms2MinNumUniqueMatches"]);
    }
    if (decodedMap.find("ms2PpmTolr") != decodedMap.end()) {
        qqqSearchParameters->ms2PpmTolr = stof(decodedMap["ms2PpmTolr"]);
    }
    if (decodedMap.find("ms2MinIntensity") != decodedMap.end()) {
        qqqSearchParameters->ms2MinIntensity = stof(decodedMap["ms2MinIntensity"]);
    }

    // END SearchParameters
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

    // END QQQSearchParameters

    return qqqSearchParameters;
}
