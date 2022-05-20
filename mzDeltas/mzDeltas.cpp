// dump basic scan information
#include "mzSample.h"
#include "mzUtils.h"
#include <iomanip>
#include <list>
#include <fstream>
#include <stdio.h>

using namespace std;

// default Conservative Params
double ppmWindow = 1.0;
double minIntensity = 1e3;
int   scanHistorySize = 10;
double minCorrelation = 0.9;
int  maxNumMzs = 1000;
int MZ_DELTA_PRECISION = 10;


struct mzPair {
    mzPair(float m1, float m2, float c, int fScan, int lScan, int lastFS) {
        mz1 = m1;
        mz2 = m2;
        cor = c;
        firstScan = fScan;
        lastScan  = lScan;
        lastFullScanIndex = lastFS;
    }

    bool withinMatchTolr(mzPair& a, double ppm, int minScanDist = 5) {
        if ( a.lastFullScanIndex > lastFullScanIndex
                && mzUtils::ppmDist(mz1, a.mz1) < ppm
                && mzUtils::ppmDist(mz2, a.mz2) < ppm
                && a.lastFullScanIndex - lastFullScanIndex < minScanDist) {
            return true;
        } else {
            return false;
        }
    }

    float mz1;
    float mz2;
    float cor;
    int firstScan;
    int lastScan;
    int lastFullScanIndex;

};

// list of observed pairs
vector<mzPair> mzPairCache;
vector<mzPair> mzPairSpillOverCache;


void spillOver(int lowestEndScanIndex) {

    vector<mzPair> newCache;
    for (mzPair& a : mzPairCache) {
        if (a.lastFullScanIndex < lowestEndScanIndex) {
            mzPairSpillOverCache.push_back(a);
        } else {
            newCache.push_back(a);
        }
    }
    mzPairCache = newCache;
}

int updateCache(mzPair a, double ppm) {
    for (mzPair& b : mzPairCache) {
        if (b.withinMatchTolr(a, 3)) {  // gets updated
            b.lastScan = a.lastScan;
            b.lastFullScanIndex = a.lastFullScanIndex;
            return 0;
        }
    }

    mzPairCache.push_back(a);
    return 1;
}

int writeMZOutput(string outputFileName, const vector<mzPair> &spillOverCache, const vector<mzPair> &cache, string sampleName) {

    // Add data from a single sample to the output file
    ofstream outputFile;
    outputFile.precision(MZ_DELTA_PRECISION);
    outputFile.setf(ios::fixed);
    outputFile.setf(ios::showpoint);
    outputFile.open(outputFileName, ios_base::app);
    for (const mzPair& a : spillOverCache) {
        outputFile << sampleName << "\t" <<  a.mz1 << "\t" << a.mz2 << "\t" << a.cor << "\t" << a.firstScan << "\t" << a.lastScan << "\n";
    }
    for (const mzPair& a : cache) {
        outputFile << sampleName << "\t" <<  a.mz1 << "\t" << a.mz2 << "\t" << a.cor << "\t" << a.firstScan << "\t" << a.lastScan << "\n";
    }
    outputFile.close();
    return 0;
}

bool is_file_exists(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

void getDistance(string fileName, double minIntensity, double minCorrelation, int scanHistorySize, double ppmWindow, int maxNumMzs, string outputFile) {

    //cout << "[mzDeltas]: Starting analysis of " << fileName << "." << endl;

    mzSample* sample = new mzSample();
    sample->loadSample(fileName.c_str());    // load file
    sample->sampleName = mzUtils::cleanFilename(fileName);

	mzPairSpillOverCache.clear();
	mzPairCache.clear();

    map<float, map<float, int> > paircount;
    list<Scan*>fullscans;

    int ms1ScanIndex = 0;
    for (Scan* scan : sample->scans) {
        if (scan->mslevel >= 2 ) continue;
        ms1ScanIndex++;
        int newPairCount = 0;
        int extendedPairCount = 0;

        fullscans.push_back(scan);
        if (fullscans.size() > scanHistorySize) fullscans.pop_front();
        if (fullscans.size() < scanHistorySize) continue;

        // get hight peaks
        vector<double>highpeaks;
        vector<int> intensityOrder = scan->intensityOrderDesc();
        for (int p = 0; p < min(intensityOrder.size(), (size_t) maxNumMzs); p++)
        {
            int j = intensityOrder[p];
            if (scan->intensity[j] > minIntensity) {
                highpeaks.push_back(scan->mz[j]);
            }
        }

    //sort high intensity by mz
    std::sort(highpeaks.begin(), highpeaks.end());

        for (int i = 0; i < highpeaks.size(); i++) {
            float mz1 = highpeaks[i];
            for (int j = i + 1; j < highpeaks.size(); j++) {
                float mz2 = highpeaks[j];

                vector<float>X;
                vector<float>Y;

                for (Scan* f : fullscans) {

                    int pos1 = f->findHighestIntensityPos(mz1, ppmWindow);
                    int pos2 = f->findHighestIntensityPos(mz2, ppmWindow);

                    float x_intensity = (pos1 != -1) ? f->intensity[pos1] : 0;
                    float y_intensity = (pos2 != -1) ? f->intensity[pos2] : 0;

                    X.push_back(x_intensity);
                    Y.push_back(y_intensity);

                }

                double corr = mzUtils::correlation(X, Y);
                if (corr > minCorrelation) {
                    // create a pair
                    mzPair a(mz1, mz2, mzUtils::correlation(X, Y), scan->scannum - scanHistorySize ,  scan->scannum,  ms1ScanIndex );
                    // check if pair exists in our cache
                    int inserted = updateCache(a, ppmWindow);
                    if (inserted == 1 ) {
                        newPairCount += 1;
                    } else {
                        extendedPairCount += 1;
                    }

                }
            }
        }



        // cerr << "Scan: " << scan->scannum << "\trt:" << scan->rt << "\tcache: " << mzPairSpillOverCache.size() << "\t" << mzPairCache.size() << "\tnew=" << newPairCount << "\text: " << extendedPairCount << "\tinputMzs:" << highpeaks.size() << endl;

        // spill
        if (ms1ScanIndex % 10 == 0) {
            spillOver(ms1ScanIndex - 10);
        }
    }

    // dump
    writeMZOutput(outputFile, mzPairSpillOverCache, mzPairCache, sample->sampleName);

    delete(sample);
    mzPairCache.clear();
    mzPairSpillOverCache.clear();

    //cout << "[mzDeltas]: Finished analysis of " << fileName << "." << endl << endl;
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "Usage: ./mzdeltas filename1.mzXML filename2.mzXML ... ";
    }

    string outputFileName = "";

    for (int i = 1; i < argc; i++) {
        string optionString(argv[i]);
        if (ends_with(optionString, "--minintensity")) {
            minIntensity = atof(argv[i + 1]);
        }
        if (ends_with(optionString, "--mincor")) minCorrelation = atof(argv[i + 1]);
        if (ends_with(optionString, "--historylen")) scanHistorySize = atoi(argv[i + 1]);
        if (ends_with(optionString, "--ppm")) ppmWindow = atoi(argv[i + 1]);
        if (ends_with(optionString, "--max_mzs")) maxNumMzs = atoi(argv[i + 1]);
        if (ends_with(optionString, "--output")) outputFileName = argv[i + 1];
    }

    const char* outputFileNameChars = outputFileName.c_str();
    if (is_file_exists(outputFileNameChars)) {
        cerr << "file \"" << outputFileNameChars << "\" already exists! It will be overwritten." << endl;
        remove(outputFileNameChars);
    }

    if (minCorrelation > 1) {
        minCorrelation = 1;
    }

    for (int i = 1; i < argc; i++) {
        string optionString(argv[i]);
        if (ends_with(optionString, "mzXML") || ends_with(optionString, "mzML")) {
            string mzXMLFileName(argv[i]);
            //cerr << mzXMLFileName << endl;
            getDistance(mzXMLFileName, minIntensity, minCorrelation, scanHistorySize, ppmWindow, maxNumMzs, outputFileName);
        }

        if (mzUtils::ends_with(optionString, "/")) {
            vector<string> dirFileNames = mzUtils::getMzSampleFilesFromDirectory(optionString.c_str());

            for (auto mzXMLFileName : dirFileNames) {
                //cerr << mzXMLFileName << endl;
                getDistance(mzXMLFileName, minIntensity, minCorrelation, scanHistorySize, ppmWindow, maxNumMzs, outputFileName);
            }
        }
    }

    cout << "mzDeltas successfully completed without error" << endl;

    return (0);
}

