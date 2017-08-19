#include "MSReader.h"
#include "mzSample.h"

using namespace MSToolkit;

int main(int argc, char** argv) {
    if (argc == 1) { cerr << "Specifiy massspec sample"; return -1; }

   char* filename = argv[1];
   mzSample* sample = new mzSample();
   sample->loadMsToolsSample(filename);

   for(Scan* scan: sample->scans) {
        cerr << scan->mslevel << "\t" << scan->totalIntensity() << endl;
   }

   return 0;
}

