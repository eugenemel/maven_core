#include "directinfusionprocessor.h"

using namespace std;

void DirectInfusionProcessor::processSingleSample(mzSample* sample, const vector<Compound*>& compounds) {

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

        int numScansPerPrecursorMz = 0;
        for (scanIterator it = scansAtKey.first; it != scansAtKey.second; ++it) {
            numScansPerPrecursorMz++;
        }

        cerr << "mapKey= " << mapKey << ": " << numScansPerPrecursorMz << " MS2 scans." << endl;
    }

    cerr << "Finished DirectInfusionProcessor::processSingleSample()" << endl;

}
