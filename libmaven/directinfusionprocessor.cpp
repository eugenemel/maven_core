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
                float precMzMin = scan->getPrecMzMin();
                float precMzMax = scan->getPrecMzMax();

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

            if (SIGN(adduct->charge) != SIGN(compound->charge)) {
                continue;
            }

            float compoundMz = compound->precursorMz;

            if (params->ms1IsRequireAdductPrecursorMatch){

                if (compound->adductString != adduct->name){
                    continue;
                }

                //TODO: this code works for compounds that do not have an 'adductString' set, but the adduct in the name.
                //delete this eventually.
//                if(compound->name.length() < adduct->name.length() ||
//                   compound->name.compare (compound->name.length() - adduct->name.length(), adduct->name.length(), adduct->name) != 0){
//                    continue;
//                }
            } else {
                compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));
            }

            //determine which map key to associate this compound, adduct with

            for (mzRangeIterator it = directInfusionSearchSet->mzRangesByMapKey.begin(); it != directInfusionSearchSet->mzRangesByMapKey.end(); ++it) {
                int mapKey = it->first;
                pair<float, float> mzRange = it->second;

                if (compoundMz > mzRange.first && compoundMz < mzRange.second) {

                    if (directInfusionSearchSet->compoundsByMapKey.find(mapKey) == directInfusionSearchSet->compoundsByMapKey.end()){
                        directInfusionSearchSet->compoundsByMapKey.insert(make_pair(mapKey, vector<pair<Compound*, Adduct*>>()));
                    }

                    directInfusionSearchSet->compoundsByMapKey[mapKey].push_back(make_pair(compound, adduct));
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

    if (debug) cerr << "Started DirectInfusionProcessor::processSingleSample()" << endl;

    //Organize all scans by common precursor m/z

    map<int, vector<Scan*>> ms2ScansByBlockNumber = {};
    vector<Scan*> validMs1Scans;

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            int mapKey = static_cast<int>(round(scan->precursorMz+0.001f)); //round to nearest int

            if (ms2ScansByBlockNumber.find(mapKey) == ms2ScansByBlockNumber.end()) {
                ms2ScansByBlockNumber.insert(make_pair(mapKey, vector<Scan*>()));
            }

            ms2ScansByBlockNumber[mapKey].push_back(scan);
        }
        if (params->ms1IsFindPrecursorIon && scan->mslevel == 1 && scan->filterString.find(params->ms1ScanFilter) != string::npos) {
            validMs1Scans.push_back(scan);
        }
    }
    if (debug) cerr << "Performing search over map keys..." << endl;

    for (auto mapKey : directInfusionSearchSet->mapKeys){

        DirectInfusionAnnotation* directInfusionAnnotation = processBlock(mapKey,
                                                                          directInfusionSearchSet->mzRangesByMapKey[mapKey],
                                                                          sample,
                                                                          ms2ScansByBlockNumber[mapKey],
                                                                          validMs1Scans,
                                                                          directInfusionSearchSet->compoundsByMapKey[mapKey],
                                                                          params,
                                                                          debug);

        if (directInfusionAnnotation) annotations.insert(make_pair(mapKey, directInfusionAnnotation));

    }

    return annotations;

}

DirectInfusionAnnotation* DirectInfusionProcessor::processBlock(int blockNum,
                                       const pair<float, float>& mzRange,
                                       mzSample* sample,
                                       const vector<Scan*>& ms2Scans,
                                       const vector<Scan*>& ms1Scans,
                                       const vector<pair<Compound*, Adduct*>> library,
                                       const shared_ptr<DirectInfusionSearchParameters> params,
                                       const bool debug){

    //need MS2 scans and compounds to identify matches
    if (ms2Scans.empty()) return nullptr;
    if (library.empty()) return nullptr;

    //build search spectrum
    Fragment *f = nullptr;
    Scan* representativeScan = nullptr;
    for (auto& scan : ms2Scans) {
        if (!f){
            f = new Fragment(scan,
                             params->scanFilterMinFracIntensity,
                             params->scanFilterMinSNRatio,
                             params->scanFilterMaxNumberOfFragments,
                             params->scanFilterBaseLinePercentile,
                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                             params->scanFilterPrecursorPurityPpm,
                             params->scanFilterMinIntensity);
            representativeScan = scan;
        } else {
            Fragment *brother = new Fragment(scan,
                                             params->scanFilterMinFracIntensity,
                                             params->scanFilterMinSNRatio,
                                             params->scanFilterMaxNumberOfFragments,
                                             params->scanFilterBaseLinePercentile,
                                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                             params->scanFilterPrecursorPurityPpm,
                                             params->scanFilterMinIntensity);

            f->addFragment(brother);
        }
    }

    f->buildConsensus(params->consensusPpmTolr,
                      params->consensusIntensityAgglomerationType,
                      params->consensusIsIntensityAvgByObserved,
                      params->consensusIsNormalizeTo10K,
                      params->consensusMinNumMs2Scans,
                      params->consensusMinFractionMs2Scans
                      );

    f->consensus->sortByMz();

    vector<shared_ptr<DirectInfusionMatchData>> libraryMatches;

    //Compare to library
    for (auto libraryEntry : library){

        unique_ptr<DirectInfusionMatchAssessment> matchAssessment = assessMatch(f, ms1Scans, libraryEntry, params, debug);
        FragmentationMatchScore s = matchAssessment->fragmentationMatchScore;
        float fragmentMaxObservedIntensity = matchAssessment->fragmentMaxObservedIntensity;

        if (s.numMatches >= params->ms2MinNumMatches &&
                s.numDiagnosticMatches >= params->ms2MinNumDiagnosticMatches &&
                params->isDiagnosticFragmentMapAgreement(matchAssessment->diagnosticFragmentMatchMap)) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());

            directInfusionMatchData->compound = libraryEntry.first;
            directInfusionMatchData->adduct = libraryEntry.second;
            directInfusionMatchData->fragmentationMatchScore = s;
            directInfusionMatchData->fragmentMaxObservedIntensity = fragmentMaxObservedIntensity;

            libraryMatches.push_back(directInfusionMatchData);
        }
    }

    //agglomerate (if necessary), and return if valid matches exist.
    if (!libraryMatches.empty()){

        //output
        DirectInfusionAnnotation *directInfusionAnnotation = new DirectInfusionAnnotation();
        directInfusionAnnotation->precMzMin = mzRange.first;
        directInfusionAnnotation->precMzMax = mzRange.second;
        directInfusionAnnotation->sample = sample;
        directInfusionAnnotation->scan = representativeScan;
        directInfusionAnnotation->fragmentationPattern = f;

        //determine fragment match maps, and mutate compounds, if needed.
        unique_ptr<DirectInfusionMatchInformation> matchInfo = DirectInfusionProcessor::getMatchInformation(
                    libraryMatches,
                    f->consensus,
                    params,
                    debug);

        //Issue 210
        DirectInfusionProcessor::addBlockSpecificMatchInfo(
                 libraryMatches,
                 matchInfo.get(),
                 f->consensus,
                 params,
                 debug);

        if (params->spectralCompositionAlgorithm == SpectralCompositionAlgorithm::ALL_CANDIDATES) {
            directInfusionAnnotation->compounds = libraryMatches;
        } else {
            directInfusionAnnotation->compounds = DirectInfusionProcessor::determineComposition(
                        libraryMatches,
                        f->consensus,
                        params,
                        debug);
        }

        //TODO: more filtering based on unique fragments criteria (or other criteria)

        return directInfusionAnnotation;
    }

    return nullptr;
}

unique_ptr<DirectInfusionMatchAssessment> DirectInfusionProcessor::assessMatch(const Fragment *f,
                                                             const vector<Scan *> &ms1Scans,
                                                             const pair<Compound*, Adduct*>& libraryEntry,
                                                             const shared_ptr<DirectInfusionSearchParameters> params,
                                                             const bool debug){
    //Initialize output
    unique_ptr<DirectInfusionMatchAssessment> directInfusionMatchAssessment = unique_ptr<DirectInfusionMatchAssessment>(new DirectInfusionMatchAssessment());

    Compound* compound = libraryEntry.first;
    Adduct *adduct = libraryEntry.second;

    //=============================================== //
    //START COMPARE MS1
    //=============================================== //

    bool isPassesMs1PrecursorRequirements = true;

     if (params->ms1IsFindPrecursorIon) {

         double precMz = compound->precursorMz;
         if (!params->ms1IsRequireAdductPrecursorMatch) {

             //Compute this way instead of using compound->precursorMz to allow for possibility of matching compound to unexpected adduct
             MassCalculator massCalc;
             float compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));
             precMz = adduct->computeAdductMass(compoundMz);

         }

         double minMz = precMz - precMz*params->ms1PpmTolr/1e6;
         double maxMz = precMz + precMz*params->ms1PpmTolr/1e6;

         isPassesMs1PrecursorRequirements = false;

         for (auto& scan : ms1Scans) {

             vector<int> matchingMzs = scan->findMatchingMzs(minMz, maxMz);

             for (auto x : matchingMzs) {
                 if (scan->intensity[x] >= params->ms1MinIntensity) {
                     isPassesMs1PrecursorRequirements = true;
                     break;
                 }
             }

             //no need to check other MS1 scans once a valid precursor has been found.
             if (isPassesMs1PrecursorRequirements) break;
         }

     }

     if (!isPassesMs1PrecursorRequirements) return directInfusionMatchAssessment; // will return with no matching fragments, 0 for every score

    //=============================================== //
    //END COMPARE MS1
    //=============================================== //

    //=============================================== //
    //START COMPARE MS2
    //=============================================== //

    Fragment t;
    t.precursorMz = compound->precursorMz;
    t.mzs = compound->fragment_mzs;
    t.intensity_array = compound->fragment_intensity;
    t.fragment_labels = compound->fragment_labels;

    float maxDeltaMz = (params->ms1PpmTolr * static_cast<float>(t.precursorMz))/ 1000000;
    directInfusionMatchAssessment->fragmentationMatchScore.ranks = Fragment::findFragPairsGreedyMz(&t, f->consensus, maxDeltaMz);

    bool isHasLabels = compound->fragment_labels.size() == directInfusionMatchAssessment->fragmentationMatchScore.ranks.size();

    float fragmentMaxObservedIntensity = 0;

    map<string, int> diagnosticMatchesMap = {};
    for (auto it = params->ms2MinNumDiagnosticMatchesMap.begin(); it != params->ms2MinNumDiagnosticMatchesMap.end(); ++it){
        diagnosticMatchesMap.insert(make_pair(it->first, 0));
    }

    for (unsigned long i=0; i < directInfusionMatchAssessment->fragmentationMatchScore.ranks.size(); i++) {

        int y = directInfusionMatchAssessment->fragmentationMatchScore.ranks[i];

        if (y != -1 && f->consensus->intensity_array[y] >= params->ms2MinIntensity) {

            float fragmentObservedIntensity = f->consensus->intensity_array[y];

            if (fragmentObservedIntensity > fragmentMaxObservedIntensity) {
                fragmentMaxObservedIntensity = fragmentObservedIntensity;
            }

            directInfusionMatchAssessment->fragmentationMatchScore.numMatches++;

            if (!isHasLabels) continue;

            if (compound->fragment_labels[i].find("*") == 0) {
                directInfusionMatchAssessment->fragmentationMatchScore.numDiagnosticMatches++;
            }

            for (auto it = params->ms2MinNumDiagnosticMatchesMap.begin(); it != params->ms2MinNumDiagnosticMatchesMap.end(); ++it){
                string diagnosticFragLabel = it->first;
                if (compound->fragment_labels[i].find(diagnosticFragLabel) == 0) {
                    diagnosticMatchesMap[diagnosticFragLabel]++;
                }
            }

        }
    }

    directInfusionMatchAssessment->diagnosticFragmentMatchMap = diagnosticMatchesMap;
    directInfusionMatchAssessment->fragmentMaxObservedIntensity = fragmentMaxObservedIntensity;

    //=============================================== //
    //END COMPARE MS2
    //=============================================== //

    return directInfusionMatchAssessment;
}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::getFragmentMatchMaps(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

       if (debug) cerr << "DirectInfusionProcessor::getFragmentMatchMaps()" << endl;

       unique_ptr<DirectInfusionMatchInformation> matchInfo = unique_ptr<DirectInfusionMatchInformation>(new DirectInfusionMatchInformation());

       for (auto directInfusionMatchData : allCandidates) {

           Compound *compound = directInfusionMatchData->compound;
           FragmentationMatchScore fragmentationMatchScore = directInfusionMatchData->fragmentationMatchScore;

           vector<int> compoundFrags(static_cast<unsigned int>(fragmentationMatchScore.numMatches));

           unsigned int matchCounter = 0;
           for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

               int observedIndex = fragmentationMatchScore.ranks[i];

               //Issue 209: peaks may be unmatched based on intensity as well as ranks[] position
               if (observedIndex == -1 || observedSpectrum->intensity_array[observedIndex] < params->ms2MinIntensity) continue;

               if (debug) cerr << "allCandidates[" << i << "]: " << compound->name << "|" << compound->adductString << " observedIndex=" << observedIndex << endl;

               int fragInt = mzToIntKey(compound->fragment_mzs[i], 1000000);

               compoundFrags[matchCounter] = fragInt;
               matchCounter++;

               pair<int, shared_ptr<DirectInfusionMatchData>> key = make_pair(fragInt, directInfusionMatchData);

               matchInfo->fragToTheoreticalIntensity.insert(make_pair(key, (compound->fragment_intensity[i])));

               matchInfo->fragToObservedIntensity.insert(make_pair(key, observedSpectrum->intensity_array[observedIndex]));

               fragToMatchDataIterator it = matchInfo->fragToMatchData.find(fragInt);

               if (it != matchInfo->fragToMatchData.end()) {
                   matchInfo->fragToMatchData[fragInt].push_back(directInfusionMatchData);
               } else {
                   vector<shared_ptr<DirectInfusionMatchData>> matchingCompounds(1);
                   matchingCompounds[0] = directInfusionMatchData;
                   matchInfo->fragToMatchData.insert(make_pair(fragInt, matchingCompounds));
               }

   //            if (debug) cerr << "allCandidates [end] i=" << i << ", ranks=" << fragmentationMatchScore.ranks.size() << endl;
           }

           matchInfo->matchDataToFrags.insert(make_pair(directInfusionMatchData, compoundFrags));

       }

       return matchInfo;

}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::summarizeByAcylChainsAndSumComposition(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        unique_ptr<DirectInfusionMatchInformation> matchInfo,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug) {

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
            summarizedCompound->setExactMass(compounds.at(0)->getExactMass());
            summarizedCompound->charge = compounds.at(0)->charge;
            summarizedCompound->db = "summarized";
            summarizedCompound->id = summarizedCompound->name + summarizedCompound->adductString;

            summarizedCompound->computeSummarizedData();

            shared_ptr<DirectInfusionMatchData> summarizedMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());
            summarizedMatchData->compound = summarizedCompound;
            summarizedMatchData->adduct = candidate->adduct;
            summarizedMatchData->fragmentationMatchScore = summarizedCompound->scoreCompoundHit(observedSpectrum, params->ms1PpmTolr, false);

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

                int observedIndex = fragmentationMatchScore.ranks[i];

                //Issue 209: peaks may be unmatched based on intensity as well as ranks[] position
                if (observedIndex == -1 || observedSpectrum->intensity_array[observedIndex] < params->ms2MinIntensity) continue;

                if (debug) cerr << "allCandidates[" << i << "]: " << compound->name << "|" << compound->adductString << " observedIndex=" << observedIndex << endl;

                int fragInt = mzToIntKey(compound->fragment_mzs[i], 1000000);

                compoundFrags[matchCounter] = fragInt;
                matchCounter++;

                pair<int, shared_ptr<DirectInfusionMatchData>> key = make_pair(fragInt, directInfusionMatchData);

                matchInfo->fragToTheoreticalIntensitySummarized.insert(make_pair(key, (compound->fragment_intensity[i])));

                matchInfo->fragToObservedIntensity.insert(make_pair(key, observedSpectrum->intensity_array[observedIndex]));

                fragToMatchDataIterator it = matchInfo->fragToMatchDataSummarized.find(fragInt);

                if (it != matchInfo->fragToMatchDataSummarized.end()) {
                    matchInfo->fragToMatchDataSummarized[fragInt].push_back(directInfusionMatchData);
                } else {
                    vector<shared_ptr<DirectInfusionMatchData>> matchingCompounds(1);
                    matchingCompounds[0] = directInfusionMatchData;
                    matchInfo->fragToMatchDataSummarized.insert(make_pair(fragInt, matchingCompounds));
                }
            }

            matchInfo->matchDataToFragsSummarized.insert(make_pair(directInfusionMatchData, compoundFrags));

        }
    }

    return matchInfo;
}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::getMatchInformation(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    if (debug) cerr << "DirectInfusionProcessor::getMatchInformation()" << endl;

    unique_ptr<DirectInfusionMatchInformation> matchInfo = getFragmentMatchMaps(allCandidates, observedSpectrum, params, debug);

    matchInfo = summarizeByAcylChainsAndSumComposition(allCandidates, move(matchInfo), observedSpectrum, params, debug);

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

                shared_ptr<DirectInfusionMatchData> compound = iterator->second[0];
                int fragId = iterator->first;

//                if (debug) cerr << "Found unique fragment for " << compound->compound->name << ": fragId=" << fragId << endl;

                shared_ptr<DirectInfusionSinglePeakMatchData> intensityData = matchInfo->getSinglePeakMatchData(fragId, compound);

//                if (debug) cerr << "Retrieved intensityData for " << compound->compound->name << ": fragId=" << fragId  << "." << endl;

                matchDataToFragIntensityIterator it = compoundToUniqueFragmentIntensities.find(compound);
                if (it != compoundToUniqueFragmentIntensities.end()) {
                    compoundToUniqueFragmentIntensities[compound].push_back(intensityData);
                } else {
                    vector<shared_ptr<DirectInfusionSinglePeakMatchData>> observedIntensities(1);
                    observedIntensities[0] = intensityData;
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

void DirectInfusionProcessor::addBlockSpecificMatchInfo(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        DirectInfusionMatchInformation *matchInfo,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){
    //TODO
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

        //Issue 218
        if (!f){
            f = new Fragment(directInfusionAnnotation->scan,
                             params->scanFilterMinFracIntensity,
                             params->scanFilterMinSNRatio,
                             params->scanFilterMaxNumberOfFragments,
                             params->scanFilterBaseLinePercentile,
                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                             params->scanFilterPrecursorPurityPpm,
                             params->scanFilterMinIntensity);
        } else {
            Fragment *brother = new Fragment(directInfusionAnnotation->scan,
                                             params->scanFilterMinFracIntensity,
                                             params->scanFilterMinSNRatio,
                                             params->scanFilterMaxNumberOfFragments,
                                             params->scanFilterBaseLinePercentile,
                                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                             params->scanFilterPrecursorPurityPpm,
                                             params->scanFilterMinIntensity);

            f->addFragment(brother);
        }

        //Issue 218
        if (directInfusionAnnotation->fragmentationPattern) {
            for (auto fragment : directInfusionAnnotation->fragmentationPattern->brothers) {
                if (fragment) {
                    Fragment *brother = new Fragment(fragment);
                    f->addFragment(brother);
                }
            }
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

    //Issue 218
    f->buildConsensus(params->consensusPpmTolr,
                      params->consensusIntensityAgglomerationType,
                      params->consensusIsIntensityAvgByObserved,
                      params->consensusIsNormalizeTo10K,
                      params->consensusMinNumMs2Scans,
                      params->consensusMinFractionMs2Scans
                      );
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
       groupMatchData->fragmentMaxObservedIntensity = matchData->fragmentMaxObservedIntensity;

       directInfusionGroupAnnotation->compounds.at(annotationMatchIndex) = groupMatchData;

       if (debug) {
           cerr << "Compound: " << groupMatchData->compound->name
                << ", Adduct: "<< groupMatchData->adduct->name
                << ", numMatches: " << groupMatchData->fragmentationMatchScore.numMatches
                << ", Proportion: " << groupMatchData->proportion
                << ", Fragment Max Observed Intensity: " << groupMatchData->fragmentMaxObservedIntensity
                << endl;
       }

       annotationMatchIndex++;
    }

    if (debug) {
        cerr << "Determined cross-sample proportions for " << annotationMatchIndex << " compounds." << endl;
        cerr << "=========================================" << endl;
    }

    return directInfusionGroupAnnotation;
}

