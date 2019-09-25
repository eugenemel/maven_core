#include "directinfusionprocessor.h"

using namespace std;

void DirectInfusionProcessor::processSingleSample(mzSample* sample, const vector<Compound*>& compounds) {

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            //TODO: do something
        }
    }
}
