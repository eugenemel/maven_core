#ifndef FRAGMENT_H
#define FRAGMENT_H

#include<iostream>
#include<vector>
#include<map>
#include<math.h>
#include<float.h>
#include <memory>
#include <unordered_set>

class Compound;
class PeakGroup;
class Scan;
class Adduct;
class TMT;
class DirectInfusionSearchParameters;
class mzSample;

using namespace std;

class TMT {
    public: 
	int scannum=0;
	vector<double>tmtIons;
	int tmtTags=0;
	double tmtTotalIntensity;
	double noise;
	double precursorPurity;

	void copy(const TMT& b) { 
        scannum = b.scannum;
        tmtIons = b.tmtIons;
		tmtTags = b.tmtTags;
		tmtTotalIntensity = b.tmtTotalIntensity;
		noise = b.noise;
		precursorPurity = b.precursorPurity;
    }
};



struct FragmentationMatchScore {

    double fractionMatched;
    double spearmanRankCorrelation;
    double ticMatched;

    double numMatches;
    double numDiagnosticMatches;
    double numSn1Matches; //Issue 312
    double numSn2Matches; //Issue 312
    double numSn3Matches; //Issue 455
    double numSn4Matches; //Issue 455
    double numOxidations; //Issue 455
    bool isHasPrecursorMatch;

    double ppmError;
    double mzFragError;
    double mergedScore;
    double dotProduct;
    double weightedDotProduct;
	double dotProductShuffle;
    double hypergeomScore;
    double mvhScore;
    double ms2purity;
    vector<double> matchedQuantiles;
    vector<int> ranks;

    //Issue 547: Additional scoring metrics
    double dotProductNorm;
    double dotProductNormNL;
    double modifiedDotProduct;
    double spectralEntropyScore;

    static vector<string> getScoringAlgorithmNames() {
        vector<string> names;
        names.push_back("HyperGeomScore");
        names.push_back("MVH");
        names.push_back("DotProduct");
        names.push_back("WeightedDotProduct");
        names.push_back("SpearmanRank");
        names.push_back("TICMatched");
        names.push_back("NumMatches");
        names.push_back("FractionRefMatched");
        names.push_back("NormDotProduct");
        names.push_back("NormNLDotProduct");
        names.push_back("ModifiedDotProduct");
        names.push_back("SpectralEntropy");
        return names;
    }


    double getScoreByName(string scoringAlgorithm) {
        if (scoringAlgorithm == "HyperGeomScore" )         return hypergeomScore;
        else if (scoringAlgorithm == "MVH")                return  mvhScore;
        else if (scoringAlgorithm == "DotProduct")         return dotProduct;
        else if (scoringAlgorithm == "SpearmanRank")       return spearmanRankCorrelation;
        else if (scoringAlgorithm == "TICMatched")         return ticMatched;
        else if (scoringAlgorithm == "WeightedDotProduct") return weightedDotProduct;
        else if (scoringAlgorithm == "NumMatches")        return  numMatches;
        else if (scoringAlgorithm == "FractionRefMatched")  return fractionMatched;
        else if (scoringAlgorithm == "NormDotProduct") return dotProductNorm;
        else if (scoringAlgorithm == "NormNLDotProduct") return dotProductNormNL;
        else if (scoringAlgorithm == "ModifiedDotProduct") return modifiedDotProduct;
        else if (scoringAlgorithm == "SpectralEntropy") return spectralEntropyScore;
        else return hypergeomScore;
    }

    FragmentationMatchScore() {
        fractionMatched=0;
        spearmanRankCorrelation=0;
        ticMatched=0;

        numMatches=0;
        numDiagnosticMatches=0;
        numSn1Matches=0;
        numSn2Matches=0;
        isHasPrecursorMatch = false;

        ppmError=DBL_MAX;
        mzFragError=DBL_MAX;
        mergedScore=0;
        dotProduct=0;
        weightedDotProduct=0;
        hypergeomScore=0;
        mvhScore=0;
        ms2purity=0;
        dotProductShuffle=0;
        dotProductNorm=0;
        dotProductNormNL=0;
        modifiedDotProduct=0;
        spectralEntropyScore=0;
    }

    FragmentationMatchScore& operator=(const FragmentationMatchScore& b) {
        fractionMatched= b.fractionMatched;
        spearmanRankCorrelation=b.spearmanRankCorrelation;
        ticMatched=b.ticMatched;

        numMatches=b.numMatches;
        numDiagnosticMatches=b.numDiagnosticMatches;
        numSn1Matches=b.numSn1Matches;
        numSn2Matches=b.numSn2Matches;

        ppmError=b.ppmError;
        mzFragError=b.mzFragError;
        mergedScore=b.mergedScore;
        dotProduct=b.dotProduct;
        weightedDotProduct=b.weightedDotProduct;
        hypergeomScore=b.hypergeomScore;
        mvhScore=b.mvhScore;
        ms2purity=b.ms2purity;
		matchedQuantiles=b.matchedQuantiles;
		dotProductShuffle = b.dotProductShuffle;
        fractionMatched = b.fractionMatched;
        ranks=b.ranks;

        dotProductNorm = b.dotProductNorm;
        dotProductNormNL = b.dotProductNormNL;
        modifiedDotProduct = b.modifiedDotProduct;
        spectralEntropyScore = b.spectralEntropyScore;

        return *this;
    }
};

class Fragment {
    public: 
        enum SortType {None=0,  Mz=1, Intensity=2 };
        enum ConsensusIntensityAgglomerationType {Mean=0, Median=1, Sum=2};

        double precursorMz;				//parent
        int polarity;					//scan polarity 	+1 or -1
        vector<float> mzs;				//mz values
        vector<float> intensity_array;		//intensity_array
        vector<string> fragment_labels; //Added for Issue #159
        vector<Fragment*> brothers;		//pointers to similar fragments 
        string sampleName;				//name of Sample
        int scanNum;					//scan Number
        float rt;						//retention time of parent scan
        float collisionEnergy;
        int precursorCharge;
        bool isDecoy;
        SortType sortedBy;
        int mergeCount;				    //number of merged events
        float purity;
		void truncateTopN(int n);
		int clusterId;
        float mergedScore;
        map<mzSample*, unordered_set<int>> scanNumMap;

        //this is only populated if isRetainOriginalScanIntensities == true in buildConsensus()
        map<int, vector<float>> consensusPositionToScanIntensities{};

        //Issue 752: NL spectra
        bool isNLSpectrum;

	TMT tmtQuant;
        PeakGroup* group;

        //consensus pattern buuld from brothers generated on buildConsensus call.
        Fragment* consensus; //consensus pattern build on brothers
        vector<int>obscount; // vector size =  mzs vector size, with counts of number of times mz was observed
        map<int,string> annotations; //mz value annotations.. assume that values are sorted by mz

        inline unsigned int nobs() { return mzs.size(); }
        inline void addFragment(Fragment* b) { brothers.push_back(b); }

        //empty constructor
        Fragment();

        //copy constructor
        Fragment(Fragment* other);

        //fragment is just a direct copy of everything in scan, with no filtering or adjustment.
        Fragment(Scan* scan);

        //build fragment based on MS2 scan, possibility of filtering and adjustment.
        Fragment(Scan* scan,
                 float minFractionalIntensity,
                 float minSigNoiseRatio,
                 int maxNumberOfFragments,
                 int baseLineLevel=5,
                 bool isRetainFragmentsAbovePrecursorMz=false,
                 float precursorPurityPpm=10,
                 float minIntensity=0
                 );

        //delete
        ~Fragment();

        void appendBrothers(Fragment* other);
        void printMzList();
        int  findClosestHighestIntensityPos(float _mz, float tolr);
        void printFragment(float productPpmToll,unsigned int limitTopX);
        void printInclusionList(bool printHeader, ostream& outstream, string COMPOUNDNAME);
        void printConsensusMS2(ostream& outstream, double minConsensusFraction);
        void printConsensusMGF(ostream& outstream, double minConsensusFraction);
        void printConsensusNIST(ostream& outstream, double minConsensusFraction, float productPpmToll, Compound* compound, Adduct* adduct);
        FragmentationMatchScore scoreMatch(Fragment* other, float productPpmToll);
        float consensusRt();
        float consensusPurity();

        map<string, int> getDiagnosticMatches(vector<string>& labels, vector<int>& ranks);

	vector<float> asDenseVector(float mzmin, float mzmax, int nbins);
        double compareToFragment(Fragment* other, float productPPMToll);
        static vector<int> compareRanks(Fragment* a, Fragment* b, float productPpmTolr);
        static vector<int> locatePositions( Fragment* a, Fragment* b, float productPpmToll);
        static vector<int> findFragPairsGreedyMz(Fragment* a, Fragment* b, float maxMzDiff, bool isMatchToNLMassShift=false);
        static double normCosineScore(Fragment* a, Fragment* b, vector<int> ranks);
        static double normNLCosineScore(Fragment* library, Fragment* observed, float productPpmTolr);

        string encodeSpectrum(int numDigits=5, string type="encoded");

        void buildConsensus(float productPpmTolr,
                            ConsensusIntensityAgglomerationType consensusIntensityAgglomerationType=Mean,
                            bool isIntensityAvgByObserved=false,
                            bool isNormalizeIntensityArray=true,
                            int minNumMs2ScansForConsensus=0,
                            float minFractionMs2ScansForConsensus=0.0f,
                            bool isRetainOriginalScanIntensities=false);

        //Issue 530: irreversible transformations
        void normalizeIntensityArray(float normValue = 10000);
        void agglomerateMzs(float minMzDelta, bool isMinMzDeltaPpm=false);
        void filterByMinIntensity(float minIntensity);
        void removeCoIsolations(float precursorMz, float ppmTolr);

        vector<unsigned int> intensityOrderDesc();
        vector<int> mzOrderInc();

        void sortByIntensity();
        void sortByMz();
		void mergeFragment(Fragment* brother, float productPpmTolr);


        double spearmanRankCorrelation(const vector<int>& X);
        double fractionMatched(const vector<int>& X);
        double mzErr(const vector<int>& X, Fragment* other);

        double totalIntensity();
        double dotProduct(Fragment* other);
	double dotProductShuffle(Fragment* other,int nbins);
        double ticMatched(const vector<int>& X);
        double mzWeightedDotProduct(const vector<int>& X, Fragment* other);
        bool hasMz(float mzValue, float ppmTolr);
        bool hasNLS(float NLS, float ppmTolr);

        void addNeutralLosses();
        void convertToNLSpectrum();

	void normalizeIntensity(vector<float>&x, int binSize);
    static int getNumDiagnosticFragmentsMatched(string fragLblStartsWith, vector<string> labels, vector<int> ranks);

        static double logNchooseK(int N,int k);
        static double SHP(int matched, int len1, int len2, int N);
        double MVH(const vector<int>& X, Fragment* other);
	vector<double> matchedRankVector(const vector<int>& X, Fragment* other);
        static bool compPrecursorMz(const Fragment* a, const Fragment* b) { return a->precursorMz<b->precursorMz; }
        bool operator<(const Fragment* b) const{ return this->precursorMz < b->precursorMz; }
        bool operator==(const Fragment* b) const{ return fabs(this->precursorMz-b->precursorMz)<0.001; }
};


#endif
