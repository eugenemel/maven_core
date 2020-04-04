#include "directinfusionprocessor.h"
#include "lipidsummarizationutils.h"

#include <chrono>

using namespace std;
using namespace mzUtils;

shared_ptr<DirectInfusionSearchSet> DirectInfusionProcessor::getSearchSet(mzSample* sample,
                                                                              const vector<Compound*>& compounds,
                                                                              const vector<Adduct*>& adducts,
                                                                              shared_ptr<DirectInfusionSearchParameters> params,
                                                                              bool debug) {

    shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet = shared_ptr<DirectInfusionSearchSet>(new DirectInfusionSearchSet());

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            int mapKey = static_cast<int>(round(scan->precursorMz+0.001f)); //round to nearest int

            if (directInfusionSearchSet->mzRangesByMapKey.find(mapKey) == directInfusionSearchSet->mzRangesByMapKey.end()) {
                float precMzMin = scan->precursorMz - 0.5f * scan->isolationWindow;
                float precMzMax = scan->precursorMz + 0.5f * scan->isolationWindow;

                directInfusionSearchSet->mzRangesByMapKey.insert(make_pair(mapKey, make_pair(precMzMin, precMzMax)));
            }

            directInfusionSearchSet->mapKeys.insert(mapKey);
        }
    }

    typedef map<int, pair<float, float>>::iterator mzRangeIterator;

    MassCalculator massCalc;

    if (debug) cerr << "Organizing database into map for fast lookup..." << endl;

    for (Compound *compound : compounds) {
        for (Adduct *adduct : adducts) {

            if (params->isRequireAdductPrecursorMatch){

                if (compound->adductString != adduct->name){
                    continue;
                }

                //TODO: this code works for compounds that do not have an 'adductString' set, but the adduct in the name.
                //delete this eventually.
//                if(compound->name.length() < adduct->name.length() ||
//                   compound->name.compare (compound->name.length() - adduct->name.length(), adduct->name.length(), adduct->name) != 0){
//                    continue;
//                }
            }

            if (SIGN(adduct->charge) != SIGN(compound->charge)) {
                continue;
            }

            float compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));

            //determine which map key to associate this compound, adduct with

            for (mzRangeIterator it = directInfusionSearchSet->mzRangesByMapKey.begin(); it != directInfusionSearchSet->mzRangesByMapKey.end(); ++it) {
                int mapKey = it->first;
                pair<float, float> mzRange = it->second;

                if (compoundMz > mzRange.first && compoundMz < mzRange.second) {
                    directInfusionSearchSet->compoundsByMapKey.insert(make_pair(mapKey, make_pair(compound, adduct)));
                    break;
                }
            }
        }
    }

    return directInfusionSearchSet;

}

map<int, DirectInfusionAnnotation*> DirectInfusionProcessor::processSingleSample(mzSample* sample,
                                                                              shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
                                                                              shared_ptr<DirectInfusionSearchParameters> params,
                                                                              bool debug) {

    MassCalculator massCalc;
    map<int, DirectInfusionAnnotation*> annotations = {};

    double totalTimeBuildConsensus = 0;
    double totalTimeScoringHits = 0;
    double totalTimeMatchingSpectra = 0;
    double totalTimeFindingMs1 = 0;

    if (debug) cerr << "Started DirectInfusionProcessor::processSingleSample()" << endl;

    //Organize all scans by common precursor m/z

    map<int, vector<Scan*>> ms2ScansByBlockNumber = {};
    vector<Scan*> validMs1Scans;

    typedef multimap<int, pair<Compound*, Adduct*>>::iterator compoundsIterator;

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            int mapKey = static_cast<int>(round(scan->precursorMz+0.001f)); //round to nearest int

            if (ms2ScansByBlockNumber.find(mapKey) == ms2ScansByBlockNumber.end()) {
                ms2ScansByBlockNumber.insert(make_pair(mapKey, vector<Scan*>()));
            }

            ms2ScansByBlockNumber[mapKey].push_back(scan);
        }
        if (params->isFindPrecursorIonInMS1Scan && scan->mslevel == 1 && scan->filterString.find(params->ms1ScanFilter) != string::npos) {
            validMs1Scans.push_back(scan);
        }
    }
    if (debug) cerr << "Performing search over map keys..." << endl;

    for (auto mapKey : directInfusionSearchSet->mapKeys){

        //need MS2 scans to identify matches
        if (ms2ScansByBlockNumber.find(mapKey) == ms2ScansByBlockNumber.end()) continue;

        pair<float,float> mzRange = directInfusionSearchSet->mzRangesByMapKey.at(mapKey);

        float precMzMin = mzRange.first;
        float precMzMax = mzRange.second;

        DirectInfusionAnnotation *directInfusionAnnotation = new DirectInfusionAnnotation();
        directInfusionAnnotation->precMzMin = precMzMin;
        directInfusionAnnotation->precMzMax = precMzMax;
        directInfusionAnnotation->sample = sample;

        vector<shared_ptr<DirectInfusionMatchData>> dIAnnotatedCompounds;

        if (debug) {
            cerr << "=========================================" << endl;
            cerr << "Investigating precMzRange = [" << precMzMin << " - " << precMzMax << "]" << endl;
        }

        vector<Scan*> scans = ms2ScansByBlockNumber[mapKey];

        Fragment *f = nullptr;
        int numScansPerPrecursorMz = 0;
        for (auto scan : scans) {
            if (numScansPerPrecursorMz == 0){
                f = new Fragment(scan, 0, 0, UINT_MAX);

                directInfusionAnnotation->scan = scan;

            } else {
                Fragment *brother = new Fragment(scan, 0, 0, UINT_MAX);
                f->addFragment(brother);
            }
            numScansPerPrecursorMz++;
        }

        if (!f) {
            delete(directInfusionAnnotation);
            continue;
        }

        //Issue 192: time building consensus
        auto startBuildConsensus = std::chrono::system_clock::now();

        f->buildConsensus(params->productPpmTolr); //TODO: a separate parameter?        
        f->consensus->sortByMz();

        //Issue 192: time building consensus
        auto stopBuildConsensus = std::chrono::system_clock::now();
        std::chrono::duration<double> buildConsensusTime = stopBuildConsensus-startBuildConsensus;
        totalTimeBuildConsensus += buildConsensusTime.count();

        directInfusionAnnotation->fragmentationPattern = f;

        pair<compoundsIterator, compoundsIterator> compoundMatches = directInfusionSearchSet->compoundsByMapKey.equal_range(mapKey);

        //check for ID bug
        if (debug) cerr << "Precursor m/z of fragment spectrum: " << f->consensus->precursorMz << endl;

        int compCounter = 0;
        int matchCounter = 0;
        for (compoundsIterator it = compoundMatches.first; it != compoundMatches.second; ++it){

            Compound* compound = it->second.first;
            Adduct *adduct = it->second.second;

            if (debug) cerr << "Scoring compound hit: " <<  compound->name << "<--> f=" << f << endl;

            //Issue 192: time scoring hits
            auto startScoringHit = std::chrono::system_clock::now();

            //Issue 192: avoid unused scoring metrics
            Fragment t;
            t.precursorMz = compound->precursorMz;
            t.mzs = compound->fragment_mzs;
            t.intensity_array = compound->fragment_intensity;
            t.fragment_labels = compound->fragment_labels;

            FragmentationMatchScore s;

            float maxDeltaMz = (params->productPpmTolr * static_cast<float>(t.precursorMz))/ 1000000;

            //Issue 192: time matching spectra
            auto startMatchingSpectra = std::chrono::system_clock::now();

            s.ranks = Fragment::findFragPairsGreedyMz(&t, f->consensus, maxDeltaMz);

            //Issue 192: time building consensus
            auto stopMatchingSpectra = std::chrono::system_clock::now();
            std::chrono::duration<double> matchingSpectraTime = stopMatchingSpectra-startMatchingSpectra;
            totalTimeMatchingSpectra += matchingSpectraTime.count();

            for(int rank: s.ranks) { if(rank != -1) s.numMatches++; }

            bool isHasLabels = compound->fragment_labels.size() == s.ranks.size();

            int numMatchAboveIntensityThreshold = 0;
            int numDiagnosticMatches = 0;
            for (int i=0; i < s.ranks.size(); i++) {

                int y = s.ranks[i];

                if (y != -1 && f->consensus->intensity_array[y] >= params->productMinIntensity) {
                    numMatchAboveIntensityThreshold++;

                    //Issue 187
                    if (isHasLabels && compound->fragment_labels[i].find("*") == 0) {
                        numDiagnosticMatches++;
                    }
                }
            }

            if (debug) cerr << "numMatchAboveIntensityThreshold=" << numMatchAboveIntensityThreshold << ", numDiagnosticMatches=" << numDiagnosticMatches << endl;

            bool isPassesMs1PrecursorRequirements = true;

            if (params->isFindPrecursorIonInMS1Scan) {

                auto startFindingPrecursor = std::chrono::system_clock::now();

                double precMz = compound->precursorMz;
                if (!params->isRequireAdductPrecursorMatch) {

                    //Compute this way instead of using compound->precursorMz to allow for possibility of matching compound to unexpected adduct
                    float compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));
                    precMz = adduct->computeAdductMass(compoundMz);

                }

                double minMz = precMz - precMz*params->parentPpmTolr/1e6;
                double maxMz = precMz + precMz*params->parentPpmTolr/1e6;

                isPassesMs1PrecursorRequirements = false;

                for (auto scan : validMs1Scans) {

                    vector<int> matchingMzs = scan->findMatchingMzs(minMz, maxMz);

                    for (auto x : matchingMzs) {
                        if (scan->intensity[x] >= params->parentMinIntensity) {
                            isPassesMs1PrecursorRequirements = true;
                            break;
                        }
                    }

                    //no need to check other MS1 scans once a valid precursor has been found.
                    if (isPassesMs1PrecursorRequirements) break;
                }

                auto stopFindingPrecursor = std::chrono::system_clock::now();

                std::chrono::duration<double> findMs1Time = stopFindingPrecursor-startFindingPrecursor;
                totalTimeFindingMs1 += findMs1Time.count();
            }

            //Issue 192: time scoring hits
            auto stopScoringHit = std::chrono::system_clock::now();
            std::chrono::duration<double> scoringHitTime = stopScoringHit-startScoringHit;
            totalTimeScoringHits += scoringHitTime.count();

            if (numMatchAboveIntensityThreshold >= params->minNumMatches && numDiagnosticMatches >= params->minNumDiagnosticFragments && isPassesMs1PrecursorRequirements) {

                if (debug) cerr << "Retain " << compound->name << ": " << s.numMatches << " matches." << endl;

                shared_ptr<DirectInfusionMatchData> directInfusionMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());
                directInfusionMatchData->compound = compound;
                directInfusionMatchData->adduct = it->second.second;
                directInfusionMatchData->fragmentationMatchScore = s;

                dIAnnotatedCompounds.push_back(directInfusionMatchData);

                matchCounter++;
            }

            compCounter++;
        }

        if (debug) {
            cerr << "Matched " << matchCounter << "/" << compCounter << " compounds." << endl;
            cerr << "=========================================" << endl;
        }

        if (matchCounter != 0){
            if (params->spectralCompositionAlgorithm == SpectralCompositionAlgorithm::ALL_CANDIDATES) {
                directInfusionAnnotation->compounds = dIAnnotatedCompounds;
            } else {
//                if (debug) cerr << "Calling DirectInfusionProcessor::determineComposition()" << endl;
                directInfusionAnnotation->compounds = DirectInfusionProcessor::determineComposition(dIAnnotatedCompounds, f->consensus, params, debug);
            }

            annotations.insert(make_pair(mapKey, directInfusionAnnotation));
        } else {
            delete(directInfusionAnnotation->fragmentationPattern);
            delete(directInfusionAnnotation);
        }
    }

    if (debug) cerr << "Finished DirectInfusionProcessor::processSingleSample()" << endl;

    if (debug){
        cerr << "=========================================\n"
             << "DirectInfusionProcessor::processSingleSample() performance stats:"
             << "\n\tConsensus Spectrum Formation: " << to_string(totalTimeBuildConsensus) << " s"
             << "\n\tScoring Spectral Hits: " << to_string(totalTimeScoringHits) << " s"
             << "\n\t\tMatching Spectra Time: " << to_string(totalTimeMatchingSpectra) << " s"
             << "\n\t\tFind Precursor in MS1 Scans Time: " << to_string(totalTimeFindingMs1) << " s"
             << "\n=========================================" << endl;
    }

    return annotations;

}

DirectInfusionAnnotation* DirectInfusionProcessor::processBlock(int blockNum,
                                       const vector<Scan*>& ms2Scans,
                                       const vector<Scan*>& ms1Scans,
                                       const vector<pair<Compound*, Adduct*>> library,
                                       shared_ptr<DirectInfusionSearchParameters> params,
                                       bool debug){
    //TODO
    return nullptr;
}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::getMatchInformation(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    if (debug) cerr << "DirectInfusionProcessor::getMatchInformation()" << endl;

    unique_ptr<DirectInfusionMatchInformation> matchInfo = unique_ptr<DirectInfusionMatchInformation>(new DirectInfusionMatchInformation());

    for (auto directInfusionMatchData : allCandidates) {

        Compound *compound = directInfusionMatchData->compound;
        FragmentationMatchScore fragmentationMatchScore = directInfusionMatchData->fragmentationMatchScore;

        vector<int> compoundFrags(static_cast<unsigned int>(fragmentationMatchScore.numMatches));

        unsigned int matchCounter = 0;
        for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

//            if (debug) cerr << "allCandidates [start] i=" << i << ", ranks=" << fragmentationMatchScore.ranks.size() << endl;

            //skip unmatched peaks
            if (fragmentationMatchScore.ranks.at(i) == -1) continue;

            int fragInt = mzToIntKey(compound->fragment_mzs.at(i), 1000000);

            compoundFrags.at(matchCounter) = fragInt;
            matchCounter++;

            pair<int, shared_ptr<DirectInfusionMatchData>> key = make_pair(fragInt, directInfusionMatchData);

            matchInfo->fragToTheoreticalIntensity.insert(make_pair(key, (compound->fragment_intensity.at(i))));

            int observedIndex = fragmentationMatchScore.ranks.at(i);

            matchInfo->fragToObservedIntensity.insert(make_pair(key, observedSpectrum->intensity_array.at(observedIndex)));

            fragToMatchDataIterator it = matchInfo->fragToMatchData.find(fragInt);

            if (it != matchInfo->fragToMatchData.end()) {
                matchInfo->fragToMatchData[fragInt].push_back(directInfusionMatchData);
            } else {
                vector<shared_ptr<DirectInfusionMatchData>> matchingCompounds(1);
                matchingCompounds.at(0) = directInfusionMatchData;
                matchInfo->fragToMatchData.insert(make_pair(fragInt, matchingCompounds));
            }

//            if (debug) cerr << "allCandidates [end] i=" << i << ", ranks=" << fragmentationMatchScore.ranks.size() << endl;
        }

        matchInfo->matchDataToFrags.insert(make_pair(directInfusionMatchData, compoundFrags));

    }

    //Identify all cases where matchData matches to identical fragments
    //and compounds can be naturally summarized to a higher level.
    for (unsigned int i = 0; i < allCandidates.size(); i++) {

        shared_ptr<DirectInfusionMatchData> iMatchData = allCandidates[i];

        for (unsigned int j = i+1; j < allCandidates.size(); j++) {

            shared_ptr<DirectInfusionMatchData> jMatchData = allCandidates[j];

            vector<int> iFrags = matchInfo->matchDataToFrags[iMatchData];
            vector<int> jFrags = matchInfo->matchDataToFrags[jMatchData];

            if (iFrags == jFrags && iMatchData->adduct->name == jMatchData->adduct->name) {

                string iChainLengthSummary;
                if (iMatchData->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()) != iMatchData->compound->metaDataMap.end()){
                    iChainLengthSummary = iMatchData->compound->metaDataMap[LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()];
                }

                string jChainLengthSummary;
                if (jMatchData->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()) != jMatchData->compound->metaDataMap.end()){
                    jChainLengthSummary = jMatchData->compound->metaDataMap[LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()];
                }

                if (!iChainLengthSummary.empty() && !jChainLengthSummary.empty() && iChainLengthSummary == jChainLengthSummary) {
                    if (matchInfo->chainLengthSummaries.find(iChainLengthSummary) != matchInfo->chainLengthSummaries.end()){
                        matchInfo->chainLengthSummaries[iChainLengthSummary].insert(iMatchData);
                        matchInfo->chainLengthSummaries[iChainLengthSummary].insert(jMatchData);
                    } else {
                        set<shared_ptr<DirectInfusionMatchData>> matchDataSet = set<shared_ptr<DirectInfusionMatchData>>();
                        matchDataSet.insert(iMatchData);
                        matchDataSet.insert(jMatchData);
                        matchInfo->chainLengthSummaries.insert(make_pair(iChainLengthSummary, matchDataSet));
                    }

                    matchInfo->originalMatchToSummaryString.insert(make_pair(iMatchData, iChainLengthSummary));
                    matchInfo->originalMatchToSummaryString.insert(make_pair(jMatchData, iChainLengthSummary));
                }

                string iCompositionSummary;
                if (iMatchData->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()) != iMatchData->compound->metaDataMap.end()) {
                     iCompositionSummary = iMatchData->compound->metaDataMap[LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()];
                }

                string jCompositionSummary;
                if (jMatchData->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()) != jMatchData->compound->metaDataMap.end()) {
                    jCompositionSummary = jMatchData->compound->metaDataMap[LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()];
                }

                if (!iCompositionSummary.empty() && !jCompositionSummary.empty() && iCompositionSummary == jCompositionSummary) {
                    if (matchInfo->compositionSummaries.find(iCompositionSummary) != matchInfo->compositionSummaries.end()) {
                        matchInfo->compositionSummaries[iCompositionSummary].insert(iMatchData);
                        matchInfo->compositionSummaries[iCompositionSummary].insert(jMatchData);
                    } else {
                        set<shared_ptr<DirectInfusionMatchData>> matchDataSet = set<shared_ptr<DirectInfusionMatchData>>();
                        matchDataSet.insert(iMatchData);
                        matchDataSet.insert(jMatchData);
                        matchInfo->compositionSummaries.insert(make_pair(iCompositionSummary, matchDataSet));
                    }

                    //acyl chain length summarization takes precedence over composition summarization
                    if (matchInfo->originalMatchToSummaryString.find(iMatchData) == matchInfo->originalMatchToSummaryString.end()) {
                        matchInfo->originalMatchToSummaryString.insert(make_pair(iMatchData, iCompositionSummary));
                    }

                    if (matchInfo->originalMatchToSummaryString.find(jMatchData) == matchInfo->originalMatchToSummaryString.end()) {
                        matchInfo->originalMatchToSummaryString.insert(make_pair(jMatchData, iCompositionSummary));
                    }
                }

            }
        }
    }

    //build new candidates list, with summarized candidates (if applicable)
    vector<shared_ptr<DirectInfusionMatchData>> summarizedCandidates;
    set<string> addedSummaries;

    for (auto candidate : allCandidates) {
        if (matchInfo->originalMatchToSummaryString.find(candidate) != matchInfo->originalMatchToSummaryString.end()) {

            string summarizedName = matchInfo->originalMatchToSummaryString[candidate];

            //only add each summarized compound one time.
            if (addedSummaries.find(summarizedName) != addedSummaries.end()) {
                continue;
            }
            addedSummaries.insert(summarizedName);

            vector<Compound*> compounds;

            //check for chain length summary
            if (matchInfo->chainLengthSummaries.find(summarizedName) != matchInfo->chainLengthSummaries.end()) {

                set<shared_ptr<DirectInfusionMatchData>> matches = matchInfo->chainLengthSummaries[summarizedName];
                compounds.resize(matches.size());

                unsigned int counter = 0;
                for (auto match : matches) {
                    compounds[counter] = match->compound;
                    counter++;
                }

            //check for composition summary
            } else if (matchInfo->compositionSummaries.find(summarizedName) != matchInfo->compositionSummaries.end()) {

                set<shared_ptr<DirectInfusionMatchData>> matches = matchInfo->compositionSummaries[summarizedName];
                compounds.resize(matches.size());

                unsigned int counter = 0;
                for (auto match : matches) {
                    compounds[counter] = match->compound;
                    counter++;
                }

            //problem case
            } else {
                cerr << "summarizedName=" << summarizedName << " Did not match to chain or composition summaries." << endl;
                abort();
            }

            //TODO: not ever deleted now
            SummarizedCompound *summarizedCompound = new SummarizedCompound(summarizedName, compounds);

            //TODO: this is guaranteed to be equal based on how the summarized compound data maps are built
            summarizedCompound->adductString = compounds.at(0)->adductString;

            /**
             * TODO: these are not guaranteed to be equal among all children compounds.
             * A sensible value should be selected for the maven gui to function properly.
             */
            summarizedCompound->formula = compounds.at(0)->getFormula();
            summarizedCompound->precursorMz = compounds.at(0)->precursorMz;

            summarizedCompound->computeFragments();

            shared_ptr<DirectInfusionMatchData> summarizedMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());
            summarizedMatchData->compound = summarizedCompound;
            summarizedMatchData->adduct = candidate->adduct;
            summarizedMatchData->fragmentationMatchScore = summarizedCompound->scoreCompoundHit(observedSpectrum, params->productPpmTolr, false);

            summarizedCandidates.push_back(summarizedMatchData);

        } else {
            summarizedCandidates.push_back(candidate);
        }
    }

    if (summarizedCandidates.size() == allCandidates.size()) { // no summarization occurred
        matchInfo->fragToMatchDataSummarized = matchInfo->fragToMatchData;
        matchInfo->matchDataToFragsSummarized = matchInfo->matchDataToFrags;
        matchInfo->fragToTheoreticalIntensitySummarized = matchInfo->fragToTheoreticalIntensity;
    } else {

        for (auto directInfusionMatchData : summarizedCandidates) {

            Compound *compound = directInfusionMatchData->compound;
            FragmentationMatchScore fragmentationMatchScore = directInfusionMatchData->fragmentationMatchScore;

            vector<int> compoundFrags(static_cast<unsigned int>(fragmentationMatchScore.numMatches));

            unsigned int matchCounter = 0;
            for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

                if (debug) cerr << "summarizedCandidates [start] i=" << i << ", ranks=" << fragmentationMatchScore.ranks.size() << endl;

                //skip unmatched peaks
                if (fragmentationMatchScore.ranks.at(i) == -1) continue;

                int fragInt = mzToIntKey(compound->fragment_mzs.at(i), 1000000);

                compoundFrags.at(matchCounter) = fragInt;
                matchCounter++;

                pair<int, shared_ptr<DirectInfusionMatchData>> key = make_pair(fragInt, directInfusionMatchData);

                matchInfo->fragToTheoreticalIntensitySummarized.insert(make_pair(key, (compound->fragment_intensity.at(i))));

                int observedIndex = fragmentationMatchScore.ranks.at(i);

                matchInfo->fragToObservedIntensity.insert(make_pair(key, observedSpectrum->intensity_array.at(observedIndex)));

                fragToMatchDataIterator it = matchInfo->fragToMatchDataSummarized.find(fragInt);

                if (it != matchInfo->fragToMatchDataSummarized.end()) {
                    matchInfo->fragToMatchDataSummarized[fragInt].push_back(directInfusionMatchData);
                } else {
                    vector<shared_ptr<DirectInfusionMatchData>> matchingCompounds(1);
                    matchingCompounds.at(0) = directInfusionMatchData;
                    matchInfo->fragToMatchDataSummarized.insert(make_pair(fragInt, matchingCompounds));
                }

                if (debug) cerr << "summarizedCandidates [end] i=" << i << ", ranks=" << fragmentationMatchScore.ranks.size() << endl;

            }

            matchInfo->matchDataToFragsSummarized.insert(make_pair(directInfusionMatchData, compoundFrags));

        }
    }

    if (debug) {
        cerr << "Fragments --> Compounds: (" << matchInfo->matchDataToFrags.size() << " passing compounds)" << endl;

        for (fragToMatchDataIterator iterator = matchInfo->fragToMatchData.begin(); iterator != matchInfo->fragToMatchData.end(); ++iterator) {
            int frag = iterator->first;
            vector<shared_ptr<DirectInfusionMatchData>> compounds = iterator->second;
            cerr<< "frag= " << intKeyToMz(frag, 1000000) << " m/z : ";
            for (auto matchData : compounds) {
                cerr << matchData->compound->name << "|" << matchData->compound->adductString << " ";
            }
            cerr << endl;
        }

        cerr << "Compounds --> Fragments: (" << matchInfo->fragToMatchData.size() << " matched fragments)" << endl;

        for (matchDataToFragIterator iterator = matchInfo->matchDataToFrags.begin(); iterator != matchInfo->matchDataToFrags.end(); ++iterator) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = iterator->first;
            vector<int> frags = iterator->second;

            cerr << "Compound= " << directInfusionMatchData->compound->name << "|" << directInfusionMatchData->compound->adductString << ": ";
            for (auto frag : frags){
                cerr << intKeyToMz(frag, 1000000) << " ";
            }
            cerr << endl;
        }

        cerr << "Chain Summaries --> Compounds:" << endl;

        for (stringToMatchDataIterator iterator = matchInfo->chainLengthSummaries.begin(); iterator != matchInfo->chainLengthSummaries.end(); ++iterator){
            string summary = iterator->first;
            set<shared_ptr<DirectInfusionMatchData>> chainLengthMatchDataSet = iterator->second;

            cerr << "Summary= " << summary << ": ";
            for (auto chainMatch : chainLengthMatchDataSet) {
                cerr << chainMatch->compound->name << "|" << chainMatch->compound->adductString << " ";
            }
            cerr << endl;
        }

        cerr << "Composition Summaries --> Compounds:" << endl;

        for (stringToMatchDataIterator iterator = matchInfo->compositionSummaries.begin(); iterator != matchInfo->compositionSummaries.end(); ++iterator){
            string summary = iterator->first;
            set<shared_ptr<DirectInfusionMatchData>> compositionMatchDataSet = iterator->second;

            cerr << "Summary= " << summary << ": ";
            for (auto compMatch : compositionMatchDataSet) {
                cerr << compMatch->compound->name << "|" << compMatch->compound->adductString << " ";
            }
            cerr << endl;
        }

        cerr << "Summarized Fragments --> Summarized Compounds: (" << matchInfo->matchDataToFragsSummarized.size() << " passing compounds)" << endl;

        for (fragToMatchDataIterator iterator = matchInfo->fragToMatchDataSummarized.begin(); iterator != matchInfo->fragToMatchDataSummarized.end(); ++iterator) {
            int frag = iterator->first;
            vector<shared_ptr<DirectInfusionMatchData>> compounds = iterator->second;
            cerr<< "frag= " << intKeyToMz(frag, 1000000) << " m/z : ";
            for (auto matchData : compounds) {
                cerr << matchData->compound->name << "|" << matchData->compound->adductString << " ";
            }
            cerr << endl;
        }

        cerr << "Summarized Compounds --> Summarized Fragments: (" << matchInfo->fragToMatchDataSummarized.size() << " matched fragments)" << endl;

        for (matchDataToFragIterator iterator = matchInfo->matchDataToFragsSummarized.begin(); iterator != matchInfo->matchDataToFragsSummarized.end(); ++iterator) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = iterator->first;
            vector<int> frags = iterator->second;

            cerr << "Compound= " << directInfusionMatchData->compound->name << "|" << directInfusionMatchData->compound->adductString << ": ";
            for (auto frag : frags){
                cerr << intKeyToMz(frag, 1000000) << " ";
            }
            cerr << endl;
        }
    }

    return matchInfo;
}

vector<shared_ptr<DirectInfusionMatchData>> DirectInfusionProcessor::determineComposition(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    enum SpectralCompositionAlgorithm algorithm = params->spectralCompositionAlgorithm;

    unique_ptr<DirectInfusionMatchInformation> matchInfo = DirectInfusionProcessor::getMatchInformation(allCandidates, observedSpectrum, params, debug);

    if (debug) {
        cerr << "matchInfo->fragToMatchDataSummarized: " << matchInfo->fragToMatchDataSummarized.size() << " entries." << endl;
        cerr << "matchInfo->matchDataToFragsSummarized: " << matchInfo->matchDataToFragsSummarized.size() << " entries." << endl;
    }

    //TODO: refactor into class, subclass, etc
    if (algorithm == SpectralCompositionAlgorithm::AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE) {

        map<shared_ptr<DirectInfusionMatchData>, vector<shared_ptr<DirectInfusionSinglePeakMatchData>>> compoundToUniqueFragmentIntensities = {};

        for (fragToMatchDataIterator iterator = matchInfo->fragToMatchDataSummarized.begin(); iterator != matchInfo->fragToMatchDataSummarized.end(); ++iterator) {
            if (iterator->second.size() == 1) { // unique fragment

                shared_ptr<DirectInfusionMatchData> compound = iterator->second.at(0);
                int fragId = iterator->first;

//                if (debug) cerr << "Found unique fragment for " << compound->compound->name << ": fragId=" << fragId << endl;

                shared_ptr<DirectInfusionSinglePeakMatchData> intensityData = matchInfo->getSinglePeakMatchData(fragId, compound);

//                if (debug) cerr << "Retrieved intensityData for " << compound->compound->name << ": fragId=" << fragId  << "." << endl;

                matchDataToFragIntensityIterator it = compoundToUniqueFragmentIntensities.find(compound);
                if (it != compoundToUniqueFragmentIntensities.end()) {
                    compoundToUniqueFragmentIntensities[compound].push_back(intensityData);
                } else {
                    vector<shared_ptr<DirectInfusionSinglePeakMatchData>> observedIntensities(1);
                    observedIntensities.at(0) = intensityData;
                    compoundToUniqueFragmentIntensities.insert(make_pair(compound, observedIntensities));
                }

            }
        }

        map<shared_ptr<DirectInfusionMatchData>, float> results = {};
        float sumIntensity = 0;

        for (matchDataToFragIntensityIterator iterator = compoundToUniqueFragmentIntensities.begin(); iterator != compoundToUniqueFragmentIntensities.end(); ++iterator){

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = iterator->first;
            vector<shared_ptr<DirectInfusionSinglePeakMatchData>> fragIntensityDataVector = iterator->second;

            float representativeIntensity = 0;
            float maxNormalizedTheoreticalIntensity = 0;
            for (auto intensityData : fragIntensityDataVector) {
                if (intensityData->normalizedTheoreticalIntensity > maxNormalizedTheoreticalIntensity){
                    maxNormalizedTheoreticalIntensity = intensityData->normalizedTheoreticalIntensity;
                    representativeIntensity = intensityData->getIntensityRatio();
                }
            }

            sumIntensity += representativeIntensity;
            results.insert(make_pair(directInfusionMatchData, representativeIntensity));

            if (debug) {
                cerr << "Compound= " << directInfusionMatchData->compound->name << "|" << directInfusionMatchData->compound->adductString << ": ";
                for (auto frag : fragIntensityDataVector){
                    cerr << frag << " ";
                }
                cerr << "median=" << representativeIntensity << endl;
            }

        }

        vector<shared_ptr<DirectInfusionMatchData>> passingMatchData;
        for (matchDataToFloatIterator iterator = results.begin(); iterator != results.end(); ++iterator) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = iterator->first;
            float intensity = iterator->second;

            double proportion = static_cast<double>(intensity/sumIntensity);

            directInfusionMatchData->proportion = proportion;
            passingMatchData.push_back(directInfusionMatchData);

            if (debug) {
                cerr << "Compound= " << directInfusionMatchData->compound->name << "|" << directInfusionMatchData->compound->adductString <<": " << proportion << endl;
            }
        }

        return passingMatchData;

    }

    return allCandidates;
}

void DirectInfusionGroupAnnotation::clean() {
    for (map<mzSample*, DirectInfusionAnnotation*>::iterator it = annotationBySample.begin(); it != annotationBySample.end(); ++it) {

        if (it->second->fragmentationPattern) delete(it->second->fragmentationPattern);
        for (auto matchData : it->second->compounds) {
             //SummarizedCompounds are created transiently by directinfusionprocessor, Compounds are retrieved from DB.compounds
//            if (SummarizedCompound* sc = dynamic_cast<SummarizedCompound*>(matchData->compound)){
//                delete(sc);
//            }
        }
        if (it->second) delete(it->second);
    }
    annotationBySample.clear();
}

DirectInfusionGroupAnnotation* DirectInfusionGroupAnnotation::createByAverageProportions(vector<DirectInfusionAnnotation*> crossSampleAnnotations, shared_ptr<DirectInfusionSearchParameters> params, bool debug) {

    DirectInfusionGroupAnnotation *directInfusionGroupAnnotation = new DirectInfusionGroupAnnotation();

    directInfusionGroupAnnotation->precMzMin = crossSampleAnnotations.at(0)->precMzMin;
    directInfusionGroupAnnotation->precMzMax = crossSampleAnnotations.at(0)->precMzMax;

    if (debug) {
        cerr << "=========================================" << endl;
        cerr << "Merging peak groups in precMzRange = [" << directInfusionGroupAnnotation->precMzMin << " - " <<  directInfusionGroupAnnotation->precMzMax << "]" << endl;
    }

    Fragment *f = nullptr;

//    map<shared_ptr<DirectInfusionMatchData>, double, DirectInfusionMatchDataCompare> proportionSums = {};
//    map<shared_ptr<DirectInfusionMatchData>, FragmentationMatchScore, DirectInfusionMatchDataCompare> bestFragMatch = {};

    map<shared_ptr<DirectInfusionMatchData>, double, DirectInfusionMatchDataCompareByNames> proportionSums = {};
    map<shared_ptr<DirectInfusionMatchData>, FragmentationMatchScore, DirectInfusionMatchDataCompareByNames> bestFragMatch = {};

    unsigned int compoundInSampleMatchCounter = 0;

    /**
     * If a sample contains no compounds, it should be excluded from calculation for cross-sample adjusted proportions.
     *
     * Individual proportions of compounds within a single sample all sum to 1, except when no compounds are found in the sample,
     * in which case the individual proportions all sum to 0 (and should be excluded from re-calculating cross-sample contributions)
     */
    unsigned int numContributingSamples = 0;

    for (auto directInfusionAnnotation : crossSampleAnnotations){
        directInfusionGroupAnnotation->annotationBySample.insert(
                    make_pair(directInfusionAnnotation->sample,
                              directInfusionAnnotation)
                    );
        if (!f){
            f = new Fragment(directInfusionAnnotation->scan, 0, 0, UINT_MAX);
        } else {
            Fragment *brother = new Fragment(directInfusionAnnotation->scan, 0, 0, UINT_MAX);
            f->addFragment(brother);
        }

        if (debug) {
            cerr << "sample=" << directInfusionAnnotation->sample->sampleName
                 << ": " << directInfusionAnnotation->compounds.size() << " compounds."
                 << endl;
        }

        if (directInfusionAnnotation->compounds.size() > 0) {
            numContributingSamples++;
        }

        for (auto matchData : directInfusionAnnotation->compounds){

            compoundInSampleMatchCounter++;

            double runningSum = matchData->proportion;

            if (proportionSums.find(matchData) != proportionSums.end()){
                runningSum += proportionSums.at(matchData);
            } else {
                proportionSums.insert(make_pair(matchData, 0.0));
            }

            proportionSums.at(matchData) = runningSum;

            if (debug) {
                cerr << "sample=" << directInfusionAnnotation->sample->sampleName
                     << "(" << matchData->compound->name
                     << ", " << matchData->adduct->name
                     << ", proportion=" << matchData->proportion
                     << "): runningSum=" << runningSum << endl;
            }

            FragmentationMatchScore bestMatch = matchData->fragmentationMatchScore;
            if (bestFragMatch.find(matchData) != bestFragMatch.end()) {
                FragmentationMatchScore previousBestMatch = bestFragMatch.at(matchData);

                //TODO: how to decide on best match?
                if (bestMatch.hypergeomScore >= previousBestMatch.hypergeomScore){
                    bestFragMatch.at(matchData) =  bestMatch;
                }

            } else {
                bestFragMatch.insert(make_pair(matchData, bestMatch));
            }

        }

    }

    if (debug) {
        cerr << "Identified " << proportionSums.size() << " unique and " << compoundInSampleMatchCounter << " total compound-adduct pairs in all samples." << endl;
    }

    f->buildConsensus(params->productPpmTolr);
    f->consensus->sortByMz();

    directInfusionGroupAnnotation->fragmentationPattern = f;

    directInfusionGroupAnnotation->compounds.resize(proportionSums.size());

    double numSamples = static_cast<double>(directInfusionGroupAnnotation->annotationBySample.size());

    unsigned int annotationMatchIndex = 0;

    for (auto matchDataPair : proportionSums) {
       shared_ptr<DirectInfusionMatchData> groupMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());

       shared_ptr<DirectInfusionMatchData> matchData = matchDataPair.first;

       groupMatchData->compound = matchData->compound;
       groupMatchData->adduct = matchData->adduct;
       groupMatchData->proportion = matchDataPair.second / numContributingSamples;
       groupMatchData->fragmentationMatchScore = bestFragMatch.at(matchData);

       directInfusionGroupAnnotation->compounds.at(annotationMatchIndex) = groupMatchData;

       if (debug) {
           cerr << "Compound: " << groupMatchData->compound->name
                << ", Adduct: "<< groupMatchData->adduct->name
                << ", numMatches: " << groupMatchData->fragmentationMatchScore.numMatches
                << ", Proportion: " << groupMatchData->proportion << endl;
       }

       annotationMatchIndex++;
    }

    if (debug) {
        cerr << "Determined cross-sample proportions for " << annotationMatchIndex << " compounds." << endl;
        cerr << "=========================================" << endl;
    }

    return directInfusionGroupAnnotation;
}

