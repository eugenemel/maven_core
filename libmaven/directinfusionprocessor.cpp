#include "directinfusionprocessor.h"

using namespace std;

vector<DirectInfusionAnnotation*> DirectInfusionProcessor::processSingleSample(mzSample* sample,
                                                                              const vector<Compound*>& compounds,
                                                                              const vector<Adduct*>& adducts,
                                                                              bool debug) {

    MassCalculator massCalc;
    vector<DirectInfusionAnnotation*> annotations;

    float minFractionalIntensity = 0.000001f; //TODO: refactor as a parameter

    if (debug) cerr << "Started DirectInfusionProcessor::processSingleSample()" << endl;

    //Organize all scans by common precursor m/z

    multimap<int, Scan*> scansByPrecursor = {};
    map<int, pair<float,float>> mzRangeByPrecursor = {};
    set<int> mapKeys = {};

    multimap<int, pair<Compound*, Adduct*>> compoundsByMapKey = {};

    typedef multimap<int, Scan*>::iterator scanIterator;
    typedef map<int, pair<float, float>>::iterator mzRangeIterator;
    typedef multimap<int, pair<Compound*, Adduct*>>::iterator compoundsIterator;

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            int mapKey = static_cast<int>(round(scan->precursorMz+0.001f)); //round to nearest int

            if (mzRangeByPrecursor.find(mapKey) == mzRangeByPrecursor.end()) {
                float precMzMin = scan->precursorMz - 0.5f * scan->isolationWindow;
                float precMzMax = scan->precursorMz + 0.5f * scan->isolationWindow;

                mzRangeByPrecursor.insert(make_pair(mapKey, make_pair(precMzMin, precMzMax)));
            }

            mapKeys.insert(mapKey);
            scansByPrecursor.insert(make_pair(mapKey, scan));
        }
    }

    if (debug) cerr << "Organizing database into map for fast lookup..." << endl;

    //TODO: restructure this code to be processed separately
    //TODO: prefilter compound/adduct pairs? think about handling generally.

    for (Compound *compound : compounds) {
        for (Adduct *adduct : adducts) {

            //require adduct match to compound name, assumes that compounds have adduct as part of name
            if(compound->name.length() < adduct->name.length() ||
               compound->name.compare (compound->name.length() - adduct->name.length(), adduct->name.length(), adduct->name) != 0){
                continue;
            }

            float compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));

            //determine which map key to associate this compound, adduct with

            for (mzRangeIterator it = mzRangeByPrecursor.begin(); it != mzRangeByPrecursor.end(); ++it) {
                int mapKey = it->first;
                pair<float, float> mzRange = it->second;

                if (compoundMz > mzRange.first && compoundMz < mzRange.second) {
                    compoundsByMapKey.insert(make_pair(mapKey, make_pair(compound, adduct)));
                }
            }
        }
    }

    if (debug) cerr << "Performing search over map keys..." << endl;

    for (auto mapKey : mapKeys){

        pair<float,float> mzRange = mzRangeByPrecursor.at(mapKey);

        float precMzMin = mzRange.first;
        float precMzMax = mzRange.second;

        DirectInfusionAnnotation *directInfusionAnnotation = new DirectInfusionAnnotation();
        directInfusionAnnotation->precMzMin = precMzMin;
        directInfusionAnnotation->precMzMax = precMzMax;

        vector<tuple<Compound*, Adduct*, double>> dIAnnotatedCompounds;

        if (debug) {
            cerr << "=========================================" << endl;
            cerr << "Investigating precMzRange = [" << precMzMin << " - " << precMzMax << "]" << endl;
        }

        pair<scanIterator, scanIterator> scansAtKey = scansByPrecursor.equal_range(mapKey);

        Fragment *f;
        int numScansPerPrecursorMz = 0;
        for (scanIterator it = scansAtKey.first; it != scansAtKey.second; ++it) {
            if (numScansPerPrecursorMz == 0){
                f = new Fragment(it->second, minFractionalIntensity, 0, UINT_MAX);

                //TODO: allow for possibility of smarter agglomeration, instead of just taking first valid scan.
                directInfusionAnnotation->scan = it->second;

            } else {
                Fragment *brother = new Fragment(it->second, minFractionalIntensity, 0, UINT_MAX);
                f->addFragment(brother);
            }
            numScansPerPrecursorMz++;
        }

        f->buildConsensus(20); //TODO: refactor as parameter
        f->consensus->sortByMz();

        pair<compoundsIterator, compoundsIterator> compoundMatches = compoundsByMapKey.equal_range(mapKey);

        int compCounter = 0;
        int matchCounter = 0;
        for (compoundsIterator it = compoundMatches.first; it != compoundMatches.second; ++it){

            Compound* compound = it->second.first;

            FragmentationMatchScore s = compound->scoreCompoundHit(f->consensus, 20, false); //TODO: parameters

            if (s.numMatches > 3) {
                if (debug) cerr << compound->name << ": " << s.numMatches << endl;

                //TODO: this should be refined with a much better algorithm
                dIAnnotatedCompounds.push_back(make_tuple(compound, it->second.second, 0));

                matchCounter++;
            }

            compCounter++;
        }

        delete(f);

        if (debug) {
            cerr << "Matched " << matchCounter << "/" << compCounter << " compounds." << endl;
            cerr << "=========================================" << endl;
        }

        if (matchCounter != 0){
            directInfusionAnnotation->compounds = dIAnnotatedCompounds;
            annotations.push_back(directInfusionAnnotation);
        }
    }

    if (debug) cerr << "Finished DirectInfusionProcessor::processSingleSample()" << endl;

    return annotations;

}
