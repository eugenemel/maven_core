#include "directinfusionprocessor.h"

using namespace std;

vector<DirectInfusionAnnotation> DirectInfusionProcessor::processSingleSample(mzSample* sample, const vector<Compound*>& compounds, const vector<Adduct*>& adducts) {

    MassCalculator massCalc;
    vector<DirectInfusionAnnotation> annotations;

    float minFractionalIntensity = 0.000001f; //TODO: refactor as a parameter

    cerr << "Started DirectInfusionProcessor::processSingleSample()" << endl;

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

    cerr << "Organizing database into map for fast lookup..." << endl;

    for (Compound *compound : compounds) {
        for (Adduct *adduct : adducts) {

            float compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));

            //determine which map key to assocaite this compound, adduct with

            for (mzRangeIterator it = mzRangeByPrecursor.begin(); it != mzRangeByPrecursor.end(); ++it) {
                int mapKey = it->first;
                pair<float, float> mzRange = it->second;

                if (compoundMz > mzRange.first && compoundMz < mzRange.second) {
                    compoundsByMapKey.insert(make_pair(mapKey, make_pair(compound, adduct)));
                }
            }
        }
    }

    cerr << "Performing search over map keys..." << endl;

    for (auto mapKey : mapKeys){

        pair<float,float> mzRange = mzRangeByPrecursor.at(mapKey);

        float precMzMin = mzRange.first;
        float precMzMax = mzRange.second;

        cerr << "=========================================" << endl;
        cerr << "Investigating precMzRange = [" << precMzMin << " - " << precMzMax << "]" << endl;

        pair<scanIterator, scanIterator> scansAtKey = scansByPrecursor.equal_range(mapKey);

        Fragment *f;
        int numScansPerPrecursorMz = 0;
        for (scanIterator it = scansAtKey.first; it != scansAtKey.second; ++it) {
            if (numScansPerPrecursorMz == 0){
                f = new Fragment(it->second, minFractionalIntensity, 0, __FLT_MAX__);
            } else {
                Fragment *brother = new Fragment(it->second, minFractionalIntensity, 0, __FLT_MAX__);
                f->addFragment(brother);
            }
            numScansPerPrecursorMz++;
        }

        f->buildConsensus(20); //TODO: refactor as parameter
        f->consensus->sortByMz();

        pair<compoundsIterator, compoundsIterator> compoundMatches = compoundsByMapKey.equal_range(mapKey);

        for (compoundsIterator it = compoundMatches.first; it != compoundMatches.second; ++it){

            int mapKey = it->first;
            Compound* compound = it->second.first;

            FragmentationMatchScore s = compound->scoreCompoundHit(f->consensus, 20, false); //TODO: parameters

            if (s.numMatches > 3) {
                cerr << compound->name << ": " << s.numMatches << endl;
            }

        }

        delete(f);
        cerr << "=========================================" << endl;

        //cerr << "mapKey= " << mapKey << ": " << numScansPerPrecursorMz << " MS2 scans." << endl;
    }

    cerr << "Finished DirectInfusionProcessor::processSingleSample()" << endl;

    return annotations;

}
