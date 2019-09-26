#include "directinfusionprocessor.h"

using namespace std;

void DirectInfusionProcessor::processSingleSample(mzSample* sample, const vector<Compound*>& compounds) {

    double minFractionalIntensity = 0.000001;

    cerr << "Started DirectInfusionProcessor::processSingleSample()" << endl;

    //Organize all scans by common precursor m/z

    multimap<int, Scan*> scansByPrecursor = {};
    set<int> mapKeys = {};

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            int mapKey = static_cast<int>(round(scan->precursorMz+0.001)); //round to nearest int

            mapKeys.insert(mapKey);
            scansByPrecursor.insert(make_pair(mapKey, scan));
        }
    }

    typedef multimap<int, Scan*>::iterator scanIterator;

    for (auto mapKey : mapKeys){

        pair<scanIterator, scanIterator> scansAtKey = scansByPrecursor.equal_range(mapKey);

        Fragment f;
        int numScansPerPrecursorMz = 0;
        for (scanIterator it = scansAtKey.first; it != scansAtKey.second; ++it) {
            if (numScansPerPrecursorMz == 0){
                Fragment f(it->second, minFractionalIntensity, 0, FLT_MAX);
            } else {
                Fragment brother(it->second, minFractionalIntensity, 0, FLT_MAX);
                f.addFragment(&brother);
            }
            numScansPerPrecursorMz++;
        }

        f.buildConsensus(20); //TODO: refactor as parameter
        f.consensus->sortByMz();

        //TODO: actually get adductlist
        vector<Adduct*> adductList;
        adductList.push_back(MassCalculator::PlusHAdduct);

         MassCalculator massCalc;

         //TODO: finish this up, fewer hacks, think about how to communicate results back to Maven GUI.

        //TODO: this will be very slow, restructure as map based on precursor m/z
        for (Compound *compound : compounds) {
            for (Adduct *adduct : adductList) {
                double compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));
            }
        }

        //cerr << "mapKey= " << mapKey << ": " << numScansPerPrecursorMz << " MS2 scans." << endl;
    }

    cerr << "Finished DirectInfusionProcessor::processSingleSample()" << endl;

}
