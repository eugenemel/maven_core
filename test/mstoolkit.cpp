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



			//example peak detection glucose.. 
			float mzmin = 186 - 186/1e6*10;
			float mzmax = 186 + 186/1e6*10;
			float rtmin = 5;
			float rtmax = 10;
			int smoothing_window=10; //smooth 10 scans
			float maxrtdiff = 3;     //grouping limits 
			//1. get chromatograms
			vector<mzSample*>samples;
			vector<EIC*>eics;
			for(mzSample* sample: samples) {
				EIC* eic = sample->getEIC(mzmin,mzmax,rtmin,rtmax,1);//get choromatogram
				eic->getPeakPositions(smoothing_window);//find peaks
				eics.push_back(eic);
			}

			//2. group peaks across samples.. 
			vector <PeakGroup> peakgroups = EIC::groupPeaks(eics,smoothing_window,maxrtdiff);

			//3. do some filtering
			for(PeakGroup& group: peakgroups) {
				if(group.maxIntensity > 1e6) { cerr << group.maxIntensity << endl; }

				//walk the peaks
				for(Peak& peak: group.peaks) {
					if(peak.peakArea > 1e5) { cerr << "Large peak"; }
				}
			}
	}
   return 0;
}

