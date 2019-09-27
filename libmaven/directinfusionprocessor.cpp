#include "directinfusionprocessor.h"


using namespace std;


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

            float compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));

            //determine which map key to associate this compound, adduct with

            for (mzRangeIterator it = directInfusionSearchSet->mzRangesByMapKey.begin(); it != directInfusionSearchSet->mzRangesByMapKey.end(); ++it) {
                int mapKey = it->first;
                pair<float, float> mzRange = it->second;

                if (compoundMz > mzRange.first && compoundMz < mzRange.second) {
                    directInfusionSearchSet->compoundsByMapKey.insert(make_pair(mapKey, make_pair(compound, adduct)));
                }
            }
        }
    }

    return directInfusionSearchSet;

}

vector<DirectInfusionAnnotation*> DirectInfusionProcessor::processSingleSample(mzSample* sample,
                                                                              shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
                                                                              shared_ptr<DirectInfusionSearchParameters> params,
                                                                              bool debug) {

    MassCalculator massCalc;
    vector<DirectInfusionAnnotation*> annotations;

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

        vector<tuple<Compound*, Adduct*, double, FragmentationMatchScore>> dIAnnotatedCompounds;

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

        f->buildConsensus(20); //TODO: refactor as parameter
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

                //TODO: this should be refined with a much better algorithm
                dIAnnotatedCompounds.push_back(make_tuple(compound, it->second.second, 0, s));

                matchCounter++;
            }

            compCounter++;
        }

        if (debug) {
            cerr << "Matched " << matchCounter << "/" << compCounter << " compounds." << endl;
            cerr << "=========================================" << endl;
        }

        if (matchCounter != 0){
            if (params->spectralCompositionAlgorithm == SpectralDeconvolutionAlgorithm::DO_NOTHING) {
                directInfusionAnnotation->compounds = dIAnnotatedCompounds;
            } else {
                //TODO: fancier algorithm here
            }

            annotations.push_back(directInfusionAnnotation);
        } else {
            delete(directInfusionAnnotation->fragmentationPattern);
            delete(directInfusionAnnotation);
        }
    }

    if (debug) cerr << "Finished DirectInfusionProcessor::processSingleSample()" << endl;

    return annotations;

}


