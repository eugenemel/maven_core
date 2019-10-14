#include "directinfusionprocessor.h"


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

    if (debug) cerr << "Started DirectInfusionProcessor::processSingleSample()" << endl;

    //Organize all scans by common precursor m/z

    multimap<int, Scan*> scansByPrecursor = {};

    typedef multimap<int, Scan*>::iterator scanIterator;
    typedef multimap<int, pair<Compound*, Adduct*>>::iterator compoundsIterator;

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            int mapKey = static_cast<int>(round(scan->precursorMz+0.001f)); //round to nearest int

            scansByPrecursor.insert(make_pair(mapKey, scan));
        }
    }
    if (debug) cerr << "Performing search over map keys..." << endl;

    for (auto mapKey : directInfusionSearchSet->mapKeys){

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

        pair<scanIterator, scanIterator> scansAtKey = scansByPrecursor.equal_range(mapKey);

        Fragment *f = nullptr;
        int numScansPerPrecursorMz = 0;
        for (scanIterator it = scansAtKey.first; it != scansAtKey.second; ++it) {
            if (numScansPerPrecursorMz == 0){
                f = new Fragment(it->second, 0, 0, UINT_MAX);

                //TODO: allow for possibility of smarter agglomeration, instead of just taking first valid scan.
                directInfusionAnnotation->scan = it->second;

            } else {
                Fragment *brother = new Fragment(it->second, 0, 0, UINT_MAX);
                f->addFragment(brother);
            }
            numScansPerPrecursorMz++;
        }

        if (!f) {
            delete(directInfusionAnnotation);
            continue;
        }

        f->buildConsensus(params->productPpmTolr); //TODO: a separate parameter?
        f->consensus->sortByMz();

        directInfusionAnnotation->fragmentationPattern = f;

        pair<compoundsIterator, compoundsIterator> compoundMatches = directInfusionSearchSet->compoundsByMapKey.equal_range(mapKey);

        int compCounter = 0;
        int matchCounter = 0;
        for (compoundsIterator it = compoundMatches.first; it != compoundMatches.second; ++it){

            Compound* compound = it->second.first;

            FragmentationMatchScore s = compound->scoreCompoundHit(f->consensus, params->productPpmTolr, false);

            if (s.numMatches >= params->minNumMatches) {
                if (debug) cerr << compound->name << ": " << s.numMatches << endl;

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
                directInfusionAnnotation->compounds = DirectInfusionProcessor::determineComposition(dIAnnotatedCompounds, f->consensus, params->spectralCompositionAlgorithm, true);
            }

            annotations.insert(make_pair(mapKey, directInfusionAnnotation));
        } else {
            delete(directInfusionAnnotation->fragmentationPattern);
            delete(directInfusionAnnotation);
        }
    }

    if (debug) cerr << "Finished DirectInfusionProcessor::processSingleSample()" << endl;

    return annotations;

}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::getMatchInformation(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        bool debug){

    unique_ptr<DirectInfusionMatchInformation> matchInfo = unique_ptr<DirectInfusionMatchInformation>(new DirectInfusionMatchInformation());

    for (auto directInfusionMatchData : allCandidates) {

        Compound *compound = directInfusionMatchData->compound;
        FragmentationMatchScore fragmentationMatchScore = directInfusionMatchData->fragmentationMatchScore;

        vector<int> compoundFrags(static_cast<unsigned int>(fragmentationMatchScore.numMatches));

        unsigned int matchCounter = 0;
        for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

            //skip unmatched peaks
            if (fragmentationMatchScore.ranks.at(i) == -1) continue;

            int fragInt = mzToIntKey(compound->fragment_mzs.at(i), 1000);

            compoundFrags.at(matchCounter) = fragInt;
            matchCounter++;

            pair<int, Compound*> key = make_pair(fragInt, compound);

            matchInfo->fragToTheoreticalIntensity.insert(make_pair(key, (compound->fragment_intensity.at(i))));

            int observedIndex = fragmentationMatchScore.ranks.at(i);

            matchInfo->fragToObservedIntensity.insert(make_pair(key, observedSpectrum->intensity_array.at(observedIndex)));

            fragToCompoundIterator it = matchInfo->fragToCompounds.find(fragInt);

            if (it != matchInfo->fragToCompounds.end()) {
                matchInfo->fragToCompounds[fragInt].push_back(compound);
            } else {
                vector<Compound*> matchingCompounds(1);
                matchingCompounds.at(0) = compound;
                matchInfo->fragToCompounds.insert(make_pair(fragInt, matchingCompounds));
            }
        }

        matchInfo->compoundToFrags.insert(make_pair(compound, compoundFrags));

    }

    if (debug) {
        cerr << "Fragments --> Compounds: (" << matchInfo->compoundToFrags.size() << " passing compounds)" << endl;

        for (fragToCompoundIterator iterator = matchInfo->fragToCompounds.begin(); iterator != matchInfo->fragToCompounds.end(); ++iterator) {
            int frag = iterator->first;
            vector<Compound*> compounds = iterator->second;
            cerr<< "frag= " << intKeyToMz(frag, 1000) << " m/z : ";
            for (auto compound : compounds) {
                cerr << compound->name << "|" << compound->adductString << " ";
            }
            cerr << endl;
        }

        cerr << "Compounds --> Fragments: (" << matchInfo->fragToCompounds.size() << " matched fragments)" << endl;

        for (compoundToFragIterator iterator = matchInfo->compoundToFrags.begin(); iterator != matchInfo->compoundToFrags.end(); ++iterator) {
            Compound* compound = iterator->first;
            vector<int> frags = iterator->second;
            cerr << "Compound= " << compound->name << "|" << compound->adductString << ": ";
            for (auto frag : frags){
                cerr << intKeyToMz(frag, 1000) << " ";
            }
            cerr << endl;
        }
    }

    return matchInfo;
}

vector<shared_ptr<DirectInfusionMatchData>> DirectInfusionProcessor::determineComposition(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        enum SpectralCompositionAlgorithm algorithm,
        bool debug){

    unique_ptr<DirectInfusionMatchInformation> matchInfo = DirectInfusionProcessor::getMatchInformation(allCandidates, observedSpectrum, debug);

    //TODO: refactor into class, subclass, etc
    if (algorithm == SpectralCompositionAlgorithm::MEDIAN_UNIQUE) {

        map<Compound*, vector<float>> compoundToUniqueFragmentIntensities = {};

        for (fragToCompoundIterator iterator = matchInfo->fragToCompounds.begin(); iterator != matchInfo->fragToCompounds.end(); ++iterator) {
            if (iterator->second.size() == 1) { // unique fragment

                Compound *compound = iterator->second.at(0);
                int fragId = iterator->first;

                float intensityRatio = matchInfo->getIntensityRatio(fragId, compound);

                compoundToFragIntensityIterator it = compoundToUniqueFragmentIntensities.find(compound);
                if (it != compoundToUniqueFragmentIntensities.end()) {
                    compoundToUniqueFragmentIntensities[compound].push_back(intensityRatio);
                } else {
                    vector<float> observedIntensities(1);
                    observedIntensities.at(0) = intensityRatio;
                    compoundToUniqueFragmentIntensities.insert(make_pair(compound, observedIntensities));
                }

            }
        }

        map<Compound*, float> results = {};
        float sumIntensity = 0;

        for (compoundToFragIntensityIterator iterator = compoundToUniqueFragmentIntensities.begin(); iterator != compoundToUniqueFragmentIntensities.end(); ++iterator){

            Compound* compound = iterator->first;
            vector<float> fragIntensityRatios = iterator->second;

            sort(fragIntensityRatios.begin(), fragIntensityRatios.end());

            float median = 0;
            unsigned long n = fragIntensityRatios.size();
            if (n % 2 != 0) {
                median = fragIntensityRatios.at(n / 2);
            } else {
                median = 0.5f*(fragIntensityRatios.at((n-1)/2) + fragIntensityRatios.at(n/2));
            }

            sumIntensity += median;
            results.insert(make_pair(compound, median));

            if (debug) {
                cerr << "Compound= " << compound->name << "|" << compound->adductString << ": ";
                for (auto frag : fragIntensityRatios){
                    cerr << frag << " ";
                }
                cerr << "median=" << median << endl;
            }

        }

        for (compoundToFloatIterator iterator = results.begin(); iterator != results.end(); ++iterator) {

            Compound *compound = iterator->first;
            float intensity = iterator->second;

            float proportion = intensity / sumIntensity;

            if (debug) {
                cerr << "Compound= " <<compound->name << "|" << compound->adductString <<": " << proportion << endl;
            }
        }

    }

    return allCandidates;
}


