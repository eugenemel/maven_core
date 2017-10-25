#include "mzSample.h"
//#include "MSReader.h"

//using namespace MSToolkit;

int main(int argc, char** argv) {
	    if (argc == 1) { cerr << "Specifiy massspec sample"; return -1; }
	
		for(int i=1; i<argc; i++ )  {
			char* filename = argv[i];
			cerr << "Testing " << filename << endl;
			// if using MSToolkit
			// mzSample* sample = new mzSample();
			// sample->loadMsToolsSample(filename);

			mzSample* sample = new mzSample();
			sample->loadSample(filename);

			unsigned int totalPeaks=0;
			long double totalInts=0;

			for(Scan* scan: sample->scans) {
				totalPeaks += scan->nobs();
				totalInts += scan->totalIntensity();
				//cerr << scan->mslevel << "\t" << scan->totalIntensity() << endl;
				for (int j=0; j < scan->nobs(); j++ ) {
					cerr << scan->mz[j] << "\t" << scan->intensity[j] << endl;
				}
			}

			cerr << filename << "\t" << sample->scans.size() << "\tpeaks=" << totalPeaks << "\ttotalInt=" << totalInts << endl;
			delete sample;
	}
   return 0;
}

