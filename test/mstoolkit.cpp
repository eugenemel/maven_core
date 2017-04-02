#include "MSReader.h"
#include "mzSample.h"

using namespace MSToolkit;
Scan* getNextScan(MSReader& reader, int scanNum);
mzSample* loadSample(const char* filename);

int main(int argc, char** argv) {
    if (argc == 1) cerr << "Specifiy mzXML";

   char* mzxmlfilename = argv[1];
   mzSample* sample = loadSample(mzxmlfilename);

   return 0;
}

mzSample* loadSample(const char* filename) {

    mzSample* sample = new mzSample();
    string filenameString = string(filename);
    sample->sampleName = sample->cleanSampleName(filename);
    sample->fileName = filenameString;

    Spectrum spec;           // For holding spectrum.
    MSReader mstReader;

   int iFirstScan=0;
   mstReader.readFile(filename, spec, iFirstScan);

   MSSpectrumType filter = MS1;
   mstReader.setFilter(filter);

   int iFileLastScan = mstReader.getLastScan();
   cerr << sample->sampleName <<  " #scans=" << iFileLastScan << endl;

    for(int scanNum=0; scanNum<iFileLastScan; scanNum++) {

        spec.clearPeaks();
        spec.clearMZ();

        mstReader.readFile(NULL, spec,scanNum);

        //basic scan information
        Scan* scan = new Scan(sample, 
                scanNum,
                spec.getMsLevel(),
                spec.getRTime(), 0, 0);


        //precursor information
        if( spec.sizeMZ()) {
                scan->precursorMz = spec.getMZ(0);
                scan->precursorCharge = spec.getCharge();
        }

        //Activation method
        MSActivation act = spec.getActivationMethod();
        string actMethod;
        switch(act){
            case mstETD: actMethod="ETD"; break;
            case mstETDSA: actMethod="ETDSA"; break;
            case mstCID: actMethod="CID"; break;
            case mstECD: actMethod="ECD"; break;
            case mstPQD: actMethod = "PQD"; break;
            case mstHCD: actMethod = "HCD"; break;
            case mstNA: default: actMethod="UNKNOWN";
            break;
        }
        scan->activationMethod = actMethod;

        // actMethod.c_str();
        //
        vector<Peak_T>* vPeaks =  spec.getPeaks();

        if (vPeaks and vPeaks->size()) {
            for(unsigned int i=0; i<vPeaks->size(); i++) {
                scan->intensity.push_back( (vPeaks->at(i)).intensity );
                scan->mz.push_back( (vPeaks->at(i)).mz );
            }
        }

        cerr << scan->scannum 
            << "\t" << scan->mslevel
            << "\t" << scan->totalIntensity()
            << "\t" << scan->precursorMz
            << "\t" << scan->precursorCharge
            << "\t" << scan->activationMethod
            << endl;
            


        sample->addScan(scan);
    }

    //set min and max values for rt
    sample->calculateMzRtRange();

    //recalculate precursor masses
    cerr << "Recalculating Ms2 Precursor Masses" << endl;
    for(Scan* ms2scan: sample->scans) {
        if(ms2scan->mslevel==2) {
            ms2scan->precursorMz=sample->getMS1PrecurursorMass(ms2scan,20);
        }
    }

    if (mystrcasestr(filename,"blan") != NULL) {
        sample->isBlank = true;
        cerr << "Found Blank: " << filename << endl;
    }

    return sample;
}
