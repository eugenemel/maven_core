#ifndef MZSAMPLE_H
#define MZSAMPLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <sstream>
#include <cstring>
#include  <limits.h>
#include  <float.h>
#include <iomanip>
#include "assert.h"

#include "pugixml.hpp"
#include "base64.h"
#include "statistics.h"
#include "mzUtils.h"
#include "mzPatterns.h"
#include "mzFit.h"
#include "mzMassCalculator.h"
#include "Matrix.h"
#include "Fragment.h"

#ifdef ZLIB
#include <zlib.h>
#endif


#ifdef CDFPARSER
#include "../libcdfread/ms10.h"
#endif

#include "MSReader.h"

#if defined(WIN32) || defined(WIN64)
//define strncasecmp strnicmp
//define isnan(x) ((x) = (x))
#endif /* Def WIN32 or Def WIN64 */

class mzSample;
class Scan;
class Peak;
class PeakGroup;
class mzSlice;
class EIC;
class Compound;
class Adduct;
class mzLink;
class Reaction;
class MassCalculator;
class ChargedSpecies;
class Fragment;
class Isotope;
struct FragmentationMatchScore;

class PeaksSearchParameters;

using namespace pugi;
using namespace mzUtils;
using namespace std;

class mzPoint {
	public:
		mzPoint() {x=y=z=0; }
		mzPoint(double mz,double intensity) { x=mz; y=intensity; z=0;}
		mzPoint(double ix,double iy,double iz) { x=ix; y=iy; z=iz; }
		mzPoint& operator=(const mzPoint& b) { x=b.x; y=b.y; z=b.z; return *this;}
		inline double mz() { return x; }
		inline double intensity() { return y; }
		static bool compX(const mzPoint& a, const mzPoint& b ) { return a.x < b.x; }
		static bool compY(const mzPoint& a, const mzPoint& b ) { return a.y < b.y; }
		static bool compZ(const mzPoint& a, const mzPoint& b ) { return a.z < b.z; }
        double x,y,z;
};


class Scan { 
    public:

    Scan(mzSample* sample, int scannum, int mslevel, float rt, float precursorMz, int polarity);
    void deepcopy(Scan* b); //copy constructor

    inline unsigned int nobs() { return mz.size(); }
    inline mzSample* getSample() { return sample; }
    inline float getPrecMzMin() {return precursorMz - isolationWindowLowerOffset;}
    inline float getPrecMzMax() {return precursorMz + isolationWindowUpperOffset;}

    //Returns instrument-set scan min/max mzs if available.  If not, fall back to observed data.
    float getMinMz();
    float getMaxMz();

    vector<int> findMatchingMzs(float mzmin, float mzmax);

    int findHighestIntensityPos(float mz, float ppm);		//higest intensity pos
    int findHighestIntensityPosAMU(float mz, float amu_tolr); //higehst intensity pos, tol in amu

    int findHighestIntensityPos(float _mz, float ppmMz, float ppm); //highest intensity pos, careful calculation of delta mz

    int findClosestHighestIntensityPos(float mz, float amu_tolr);	//highest intensity pos nearest to the cente mz
    bool isMonoisotopicPrecursor(float monoIsotopeMz, float ppm, int charge=1);

    //returns -1 if not found
    float findNormalizedIntensity(float queryMz, float standardMz, float ppm, float minScanIntensity=0.0f);
    float findClosestMzIntensity(float queryMz, float ppm);

    bool hasMz(float mz, float ppm);
    bool isCentroided() { return centroided; }
    bool isProfile()    { return !centroided; }
    inline int getPolarity() { return polarity; }
    void  setPolarity(int x) { polarity = x; }

    double totalIntensity(){ double sum=0; for(unsigned int i=0;i<intensity.size();i++) sum += intensity[i]; return sum; }
    float maxIntensity()  { float max=0; for(unsigned int i=0;i<intensity.size();i++) if(intensity[i] > max) max=intensity[i]; return max; }
    float minMz()  { if(nobs() > 0) return mz[0]; return 0; }
    float maxMz()  { if(nobs() > 0) return mz[nobs()-1]; return 0; }
    float baseMz();

    vector<int> intensityOrderDesc(); //return postion in a scan from higest to lowerst intensity
    vector<pair<float,float> > getTopPeaks(float minFracIntensity,
                                           float minSNRatio=1,
                                           int baseLinePercentile=5,
                                           float minIntensity=0);

    vector<int>assignCharges(float ppmTolr);

    vector<float> chargeSeries(float Mx, unsigned int Zx);
    ChargedSpecies* deconvolute(float mzfocus, float noiseLevel, float ppmMerge, float minSigNoiseRatio, int minDeconvolutionCharge, int maxDeconvolutionCharge, int minDeconvolutionMass, int maxDeconvolutionMass, int minChargedStates );
	string toMGF();

    void  simpleCentroid();
    void  intensityFilter( int minIntensity);
    void  quantileFilter(int minQuantile);
    void  summary();
    void  log10Transform();

	Scan* getLastFullScan(int historySize);
	vector<mzPoint> getIsolatedRegion(float isolationWindowAmu);
	double getPrecursorPurity(float ppm);

    int mslevel;
    bool centroided;
    float rt;
    int   scannum;

    float precursorMz;
	float isolationWindow;
    float precursorIntensity;
    int precursorCharge;
    string activationMethod;
    int precursorScanNum;
    float productMz;
    float collisionEnergy;
    float injectionTime;

    vector <float> intensity;
    vector <float> mz;
    string scanType;
    string filterLine;
    string filterString = "";
    mzSample* sample;

    float lowerLimitMz = -1.0f;
    float upperLimitMz = -1.0f;

    //matters most for direct infusion data
    //assume 1Da window centered around precursorMz unless otherwise specified
    float isolationWindowLowerOffset = 0.5f;
    float isolationWindowUpperOffset = 0.5f;

    void addChildScan(Scan* s) { children.push_back(s); }
    vector<Scan*> getAllChildren() { return children; }
    Scan* getFirstChild() { if(children.size() == 0) return 0; else return children[0]; }
    TMT tmtQuant(); 

    static bool compRt(Scan* a, Scan* b ) { return a->rt < b->rt; }
    static bool compPrecursor(Scan* a, Scan* b ) { return a->precursorMz < b->precursorMz; }
    static bool compIntensity(Scan* a, Scan* b ) { return a->totalIntensity() > b->totalIntensity(); }
    bool operator< (const Scan& b) { return rt < b.rt; }  //default comparision operation

    vector<Isotope> getIsotopicPattern(float centerMz, float ppm, int maxZ, int maxIsotopes);

    string getSignature(int limitSize=200);

    //Issue 256
    float ms1PrecursorForMs3 = 0.0f;

private:
    	vector<Scan*> children;
        int polarity;


};

//Issue 347
class SRMTransition{
public:

    float precursorMz = 0.0f;
    float productMz = 0.0f;
    float rt = 0.0f;
    Compound *compound = nullptr;
    Adduct *adduct = nullptr;

    vector<pair<mzSample*, mzSlice*>> mzSlices{};

    set<mzSample*> getSamples();
};

//Issue 368
//TODO: this is currently unused, but leaving here in the event may be used in the future.
//Currenly, pair<float, float> of precursorMz, productMz are used instead of a dedicated key class.
//class SRMKey {
//public:
//    int precursorMz = 0;
//    int productMz = 0;
//    int rtStart = -1;
//    int rtEnd = -1;

//    static const int MZ_MULT_FACTOR = 4;
//    static const int RT_MULT_FACTOR = 2;

//    static SRMKey mzKey(float precursorMz, float productMz) {

//        SRMKey key;

//        int precMzInt = mzUtils::mzToIntKey(precursorMz, MZ_MULT_FACTOR);
//        int prodMzInt = mzUtils::mzToIntKey(productMz, MZ_MULT_FACTOR);

//        key.precursorMz = precMzInt;
//        key.productMz = prodMzInt;

//        return key;
//    }

//    static SRMKey mzRtKey(float precursorMz, float productMz, float rtStart, float rtEnd) {

//        SRMKey key = SRMKey::mzKey(precursorMz, productMz);

//        int rtStartKey = mzUtils::mzToIntKey(rtStart, RT_MULT_FACTOR);
//        int rtEndKey = mzUtils::mzToIntKey(rtEnd, RT_MULT_FACTOR);

//        key.rtStart = rtStartKey;
//        key.rtEnd = rtEndKey;

//        return key;
//    }
//};

class SRMTransition;

class mzSlice { 
    public:

    mzSlice(float mzmin, float mzmax, float rtmin, float rtmax) {
        this->mzmin=mzmin;
        this->mzmax=mzmax;
        this->rtmin=rtmin;
        this->rtmax=rtmax;
        this->mz=mzmin+(mzmax-mzmin)/2;
        this->rt=rtmin+(rtmax-rtmin)/2;
        compound=nullptr;
        adduct=nullptr;
        ionCount=0;
        compoundVector=vector<Compound*>();
    }

    mzSlice(float mz, float rt, float ions) {
        this->mz=this->mzmin=this->mzmax=mz;
        this->rt=this->rtmin=this->rtmax=rt;
        this->ionCount=ions;
        compound=nullptr;
        adduct=nullptr;
        ionCount=0;
        compoundVector=vector<Compound*>();
    }

    mzSlice(string filterLine) {
        mzmin=mzmax=rtmin=rtmax=mz=rt=ionCount=0;
        compound=nullptr;
        adduct=nullptr;
        srmId=filterLine;
        compoundVector=vector<Compound*>();
    }

    mzSlice() {
        mzmin=mzmax=rtmin=rtmax=mz=rt=ionCount=0;
        compound=nullptr;
        adduct=nullptr;
        compoundVector=vector<Compound*>();
    }

    mzSlice(SRMTransition srmTransition) {

        mzmin=mzmax=rtmin=rtmax=mz=rt=ionCount=0;

        compound=srmTransition.compound;
        adduct=srmTransition.adduct;

        compoundVector=vector<Compound*>();
        if (srmTransition.compound){
            compoundVector.push_back(srmTransition.compound);
        }

        srmPrecursorMz = srmTransition.precursorMz;
        srmProductMz = srmTransition.productMz;
    }

    mzSlice(const mzSlice& b) {
        mzmin=b.mzmin;
        mzmax=b.mzmax;
        rtmin=b.rtmin;
        rtmax=b.rtmax;
        ionCount=b.ionCount;
        mz=b.mz;
        rt=b.rt;
        compound=b.compound;
        adduct=b.adduct;
        srmId=b.srmId;
        compoundVector=b.compoundVector;
        srmPrecursorMz=b.srmPrecursorMz;
        srmProductMz=b.srmProductMz;
    }

    mzSlice& operator= (const mzSlice& b) {
        mzmin=b.mzmin;
        mzmax=b.mzmax;
        rtmin=b.rtmin;
        rtmax=b.rtmax;
        ionCount=b.ionCount;
        compound=b.compound;
        adduct=b.adduct;
        srmId=b.srmId;
        mz=b.mz;
        rt=b.rt;
        compoundVector=b.compoundVector;
        srmPrecursorMz=b.srmPrecursorMz;
        srmProductMz=b.srmProductMz;
        return *this;
    }

        float mzmin;
        float mzmax;
        float rtmin;
        float rtmax;
        float mz;
        float rt;
        float ionCount;
        Compound* compound;
        Adduct*   adduct;
        vector<Compound*> compoundVector;
		bool deleteFlag=0;

    //SRM-associated
	string srmId;
    float srmPrecursorMz = 0.0f;
    float srmProductMz = 0.0f;

    inline bool isSrmTransitionSlice(){return (srmPrecursorMz > 0.0f && srmProductMz > 0.0f);}
    inline pair<float, float> getSrmMzKey(){return make_pair(srmPrecursorMz, srmProductMz);}

	static bool compIntensity(const mzSlice* a, const mzSlice* b ) { return b->ionCount < a->ionCount; }
	static bool compMz(const mzSlice* a, const mzSlice* b ) { return a->mz < b->mz; }
	static bool compRt(const mzSlice* a, const mzSlice* b ) { return a->rt < b->rt; }
	bool operator< (const mzSlice* b) { return mz < b->mz; }
	bool operator< (const mzSlice& b) { return mz < b.mz; }
};

class mzLink { 
		public: 
			mzLink(){ mz1=mz2=0; value1=value2=0.0; data1=data2=NULL; correlation=0;}
			mzLink( int a, int b, string n ) { value1=a; value2=b; note=n; correlation=0; }
			mzLink( float a,float b, string n ) { mz1=a; mz2=b; note=n; correlation=0; }
			~mzLink() {}
            //link between two mzs
			float mz1;          //from 
			float mz2;          //to
			string note;        //note about this link

            //generic placeholders to attach values to links
			float value1;		
			float value2;	

            //generic  placeholders to attach objects to the link
            void* data1;   
            void* data2;
			//

            //Issue 337
            //for use with Isotope widget
            float isotopeFrac = 0.0f;
            float percentExpected = 0.0f;
            float percentRelative = 0.0f;

           float correlation;
           static bool compMz(const mzLink& a, const mzLink& b ) { return a.mz1 < b.mz1; }
           static bool compCorrelation(const mzLink& a, const mzLink& b) { return a.correlation > b.correlation; }
       };

class mzSample {
public:
    mzSample();                         			// constructor
    ~mzSample();
    void openStream(const char* filename);	   	// constructor : load from file
    void loadSample(const char* filename);      // constructor: load from file
    void loadSample(const char* filename, bool isCorrectPrecursor);	   	// constructor : load from file
    void loadMsToolsSample(const char* filename);	// constructor : load using MSToolkit library
    void parseMzData(const char*);			// load data from mzData file
    void parseMzCSV(const char*);			// load data from mzCSV file
    void parseMzXML(const char*);			// load data from mzXML file
    void parseMzML(const char*);			// load data from mzML file
    int  parseCDF (char *filename, int is_verbose);     // load netcdf files
    Scan* parseMzXMLScan(const xml_node& scan, int scannum);		// parse individual scan
    Scan* randomAccessMzXMLScan(int seek_pos_start, int seek_pos_end);
    void writeMzCSV(const char*);
    string cleanSampleName(string fileName);

    map<string,string> mzML_cvParams(xml_node node);
    void parseMzMLChromatogromList(xml_node);
    void parseMzMLSpectrumList(xml_node);
    void summary();						   	// print info about sample
    void calculateMzRtRange();                      //compute min and max values for mz and rt
    float getAverageFullScanTime();			// average time difference between scans
    void enumerateSRMScans();			//srm->scan mapping for QQQ


    float correlation(float mz1,  float mz2, float ppm, float rt1, float rt2, bool debug=false); //correlation in EIC space
    float getNormalizationConstant() { return _normalizationConstant; }
    void  setNormalizationConstant(float x) { _normalizationConstant = x; }

    Scan* getScan(unsigned int scanNum);				//get indexes for a given scan
    Scan* getAverageScan(float rtmin, float rtmax, int mslevel, int polarity, float resolution);

    EIC* getEIC(float minMz, float maxMz, float minRt, float maxRt, int mslevel, string scanFilterString="");	//get eic based on minMz, maxMz, minRt, maxRt, mslevel
    EIC* getEIC(string srmId);	//get eic based on srmId
    EIC* getEIC(float precursorMz, float collisionEnergy, float productMz, float amuQ1, float amuQ2 );
    EIC* getEIC(pair<float, float> mzKey); //get eic based on SRM precursor, product ion mzs
    EIC* getTIC(float,float,int);		//get Total Ion Chromatogram
    EIC* getBIC(float,float,int);		//get Base Peak Chromatogram
    double getMS1PrecursorMass(Scan* ms2scan,float ppm);
    vector<Scan*> getFragmentationEvents(mzSlice* slice);

    deque <Scan*> scans;
    int sampleId;
    string sampleName;
    string fileName;
    bool isSelected;
    bool isBlank;
    float color[4];	                                //sample display color, [r,g,b, alpha]
    float minMz;
    float maxMz;
    float minRt;
    float maxRt;
    float maxIntensity;
    float minIntensity;
    float totalIntensity;
    int   _sampleOrder;                              //sample display order
    bool  _C13Labeled;
    bool  _N15Labeled;
    bool  _S34Labeled; //Feng note: added to track S34 labeling state
    bool  _D2Labeled;  //Feng note: added to track D2 labeling state
    float _normalizationConstant;
    string  _setName;
    map<string,vector<int> >srmScans;		//srm to scan mapping
    map<string,string> instrumentInfo;		//tags associated with this sample
    map<pair<float, float>,vector<int>> srmScansByMzs{};    //srm to scan mapping, key = <precursorMz, productMz>

    float _averageFullScanTime;

    //saving and restoring retention times
    vector<float>originalRetentionTimes;       			//saved retention times prior to alignment
    vector<double>polynomialAlignmentTransformation;		//parameters for polynomial transform

    //Instrument PPM Calibration
    vector<double> ppmCorrectionCoef;

    void saveOriginalRetentionTimes();
    void restoreOriginalRetentionTimes();       
	void applyPolynomialTransform();

    //class functions
    void addScan(Scan*s);
    inline int getPolarity() { if( scans.size()>0) return scans[0]->getPolarity(); return 0; }
    inline unsigned int   scanCount() { return(scans.size()); }
    inline string getSampleName() { return sampleName; }
    void setSampleOrder(int x) { _sampleOrder=x; }
    inline int	  getSampleOrder() { return _sampleOrder; }
    inline string  getSetName()  { return _setName; }
    void   setSetName(string x) { _setName=x; }
    void   setSampleName(string x) { sampleName=x; }
    void   setSampleId(int id) { sampleId=id; }
    inline int	  getSampleId() { return sampleId; }

    static float getMaxRt(const vector<mzSample*>&samples);
    bool C13Labeled(){ return _C13Labeled; }
    bool N15Labeled(){ return _N15Labeled; }

    static bool compSampleOrder(const mzSample* a, const mzSample* b ) { return a->_sampleOrder < b->_sampleOrder; }
    static bool compSampleName(const mzSample* a, const mzSample* b ) { return a->sampleName < b->sampleName; }
    static mzSlice getMinMaxDimentions(const vector<mzSample*>& samples);


    static void setFilter_minIntensity(int x ) { filter_minIntensity=x; }
    static void setFilter_centroidScans( bool x) { filter_centroidScans=x; }
    static void setFilter_intensityQuantile(int x ) { filter_intensityQuantile=x; }
    static void setFilter_mslevel(int x ) { filter_mslevel=x; }
    static void setFilter_polarity(int x ) { filter_polarity=x; }

    static int getFilter_minIntensity() { return filter_minIntensity; }
    static int getFilter_intensityQuantile() { return filter_intensityQuantile; }
    static int getFilter_centroidScans() { return filter_centroidScans; }
    static int getFilter_mslevel()  { return filter_mslevel; }
    static int getFilter_polarity() { return filter_polarity; }

    vector<float> getIntensityDistribution(int mslevel);

    private:
        static int filter_minIntensity;
        static bool filter_centroidScans;
        static int filter_intensityQuantile;
		static int filter_mslevel;
        static int filter_polarity;
        ifstream   _iostream;

};

struct mzSample_name_less {
    bool operator()(mzSample* lhs, mzSample *rhs) const {
        return lhs->sampleName < rhs->sampleName;
    }
};

class EIC {
	
	public:

    EIC() {
        sample=NULL;
        mzmin=mzmax=rtmin=rtmax=0;
        maxIntensity=totalIntensity=0;
        eic_noNoiseObs=0;
        smootherType=GAUSSIAN;
        baselineSmoothingWindow=5;  //baseline smoothin window
        baselineDropTopX=60;    //fraction of point to remove when computing baseline
        for(unsigned int i=0; i<4;i++) color[i]=0;
    }

		~EIC();


        enum SmootherType { GAUSSIAN=0, AVG=1, SAVGOL=2 };
        vector <int> scannum;
		vector <float> rt;
		vector <float> mz;
		vector <float> intensity;
		vector <Peak>  peaks;
		string sampleName;

        mzSample* sample;       //pointer to originating sample
        float color[4];         //color of the eic line, [r,g,b, alpha]
        vector<float> spline;   //smoothed intensity values
        vector<float> baseline; //baseline

        float maxIntensity;     //maxItensity in eics
        float totalIntensity;   //sum of all intensities in EIC
        int   eic_noNoiseObs;   //number of observatiosn above baseline.

        float mzmin;
        float mzmax;
        float rtmin;
        float rtmax;

        float baselineQCutVal; // computed by EIC::computeBaseline(). Deals with raw intensities (not with spline).

		Peak* addPeak(int peakPos);
		void deletePeak(unsigned int i);
        void getPeakPositions(int smoothWindow);
        void getPeakPositionsB(int smoothWindow, float minSmoothedPeakIntensity);
        void getPeakPositionsC(int smoothWindow, bool debug);
        void getSingleGlobalMaxPeak(int smoothWindow);
		void getPeakDetails(Peak& peak);
		void getPeakWidth(Peak& peak);
        void computeBaseLine(int smoothingWindow, int dropTopX);
        void computeSpline(int smoothWindow);
		void findPeakBounds(Peak& peak);
		void getPeakStatistics();
		void checkGaussianFit(Peak& peak);
        vector<Scan*> getFragmentationEvents();
		void subtractBaseLine();
		void removeOverlapingPeaks();
		EIC* clone(); //make a copy of self


		vector<mzPoint> getIntensityVector(Peak& peak);
		void summary();
        void setSmootherType(EIC::SmootherType x) { smootherType=x; }
        void setBaselineSmoothingWindow(int x) { baselineSmoothingWindow=x;}
        void setBaselineDropTopX(int x) { baselineDropTopX=x; }
        void interpolate();

		inline unsigned int size() { return intensity.size();}
		inline mzSample* getSample() { return sample; } 
        static vector<PeakGroup> groupPeaks(vector<EIC*>&eics, int smoothingWindow, float maxRtDiff);
        static vector<PeakGroup> groupPeaksB(vector<EIC*>&eics, int smoothingWindow, float maxRtDiff, float minSmoothedPeakIntensity);

		static EIC* eicMerge(const vector<EIC*>& eics);
		static void removeLowRankGroups(vector<PeakGroup>&groups, unsigned int rankLimit );
		static bool compMaxIntensity(EIC* a, EIC* b ) { return a->maxIntensity > b->maxIntensity; }

        private:
                SmootherType smootherType;
                int baselineSmoothingWindow;
                int baselineDropTopX;

};

class Peak {
	public:
		Peak();
		Peak(EIC* e, int p);
        Peak(const Peak& p);
		Peak& operator=(const Peak& o);
                void copyObj(const Peak& o);

		//~Peak() { cerr << "~Peak() " << this << endl; }

		unsigned int pos;
		unsigned int minpos;
		unsigned int maxpos;

		float rt;
		float rtmin;
		float rtmax;
		float mzmin;
		float mzmax;

		unsigned int scan;
		unsigned int minscan;
		unsigned int maxscan;

		float peakArea;                    //non corrected sum of all intensities
		float peakAreaCorrected;           //baseline substracted area
		float peakAreaTop;                  //top 3points of the peak
		float peakAreaFractional;          //area of the peak dived by total area in the EIC
		float peakRank;                     //peak rank (sorted by peakAreaCorrected)

		float peakIntensity;               //not corrected intensity at top of the pek
		float peakBaseLineLevel;            //baseline level below the highest point

		float peakMz;	//mz value at the top of the peak
		float medianMz;	//averaged mz value across all points in the peak
		float baseMz;	//mz value across base of the peak

		float quality;	//from 0 to 1. indicator of peak goodness
		unsigned int width;		//width of the peak at the baseline
		float gaussFitSigma;	//fit to gaussian curve
		float gaussFitR2;		//fit to gaussian curve
		int groupNum;

		unsigned int noNoiseObs; 
		float noNoiseFraction;
		float symmetry;
		float signalBaselineRatio;
		float groupOverlap;			// 0 no overlap, 1 perfect overlap
		float groupOverlapFrac;		

		int   	chargeState;		    	
		int   	ms2EventCount;
		float   selectionScore;
		bool    isMonoIsotopic;

		bool localMaxFlag;
		bool fromBlankSample;		//true if peak is from blank sample

		char label;		//classification label
        mzSample *sample;  //pointer to sample

    private:
		EIC* eic; 		//pointer to eic 

	public: 
		void setEIC(EIC* e) { eic=e; }
		inline EIC*	 getEIC() { return eic;    }
		inline bool hasEIC() { return eic != NULL; }
		Scan* getScan() { if(sample) return sample->getScan(scan); else return NULL; }	
        vector<Scan*> getFragmentationEvents(float ppmWindow);
        Fragment getConsensusFragmentation(float ppmWindow, float productPpmTolr);

		int getChargeState();

		void   setSample(mzSample* s ) { sample=s; }
		inline mzSample* getSample() { return sample; }
		inline bool hasSample() 	{  return sample != NULL; }
		
		void setLabel(char label) { this->label=label;}
		inline char getLabel() { return label;}
		
        static bool compRtMin(const Peak& a, const Peak& b ) { return a.rtmin < b.rtmin; }
		static bool compRt(const Peak& a, const Peak& b ) { return a.rt < b.rt; }
		static bool compIntensity(const Peak& a, const Peak& b ) { return b.peakIntensity < a.peakIntensity; }
		static bool compArea(const Peak& a, const Peak& b ) { return b.peakAreaFractional < a.peakAreaFractional; }
		static bool compMz(const Peak& a, const Peak& b ) { return a.peakMz < b.peakMz; }
		static bool compSampleName(const Peak& a, const Peak& b ) { return a.sample->getSampleName() < b.sample->getSampleName(); }
		static bool compSampleOrder(const Peak& a, const Peak& b ) { return a.sample->getSampleOrder() < b.sample->getSampleOrder(); }
		inline static float overlap(const Peak& a, const Peak& b) {	return( checkOverlap(a.rtmin, a.rtmax, b.rtmin, b.rtmax)); }
		vector<mzLink> findCovariants();
};

//Issue 371: moved from maven/classifier.h
class Classifier
{

public:
    Classifier();
    virtual ~Classifier();

    //reimplemented in children
    virtual void classify(PeakGroup* grp)=0;
    virtual void train(vector<PeakGroup*>& groups)=0;
    virtual void refineModel(PeakGroup* grp)=0;
    virtual void saveModel(string filename)=0;
    virtual void loadModel(string filename)=0;
    virtual void loadDefaultModel()=0;
    virtual bool hasModel()=0;

    //common non virtal functions
    void classify(vector<PeakGroup*>& groups);
    void saveFeatures(vector<PeakGroup*>& groups, string filename);
    vector<Peak*> removeRedundacy(vector<Peak*>&peaks);

    inline string getModelFilename() { return _filename; }
    int getLabelsCount () { return labels.size(); }
    int getNumFeatures() { return num_features; }
    string getClassifierType() { return clsf_type; }
    void printLabelDistribution();

    virtual vector<float> getFeatures(Peak& p) = 0;

protected:
    string _filename;

    string clsf_type;
    int num_features;
    vector<string> features_names;
    vector< vector<float> > FEATURES;
    vector<char>labels;
};

enum IsotopeParametersType{INVALID=0,SAVED=1,FROM_GUI=2};

struct IsotopeParameters {

    string searchVersion = "8.1.27.25";
    IsotopeParametersType isotopeParametersType = IsotopeParametersType::INVALID;

    float ppm = 20;
    double maxIsotopeScanDiff = 10;
    double maxNaturalAbundanceErr = 100;
    double minIsotopicCorrelation=0;
    bool   isC13Labeled=false;
    bool   isN15Labeled=false;
    bool   isS34Labeled=false;
    bool   isD2Labeled=false;

    EIC::SmootherType eic_smoothingAlgorithm = EIC::SmootherType::GAUSSIAN;
    float eic_smoothingWindow;

    bool isIgnoreNaturalAbundance = true;
    bool isExtractNIsotopes = false;
    int maxIsotopesToExtract = 5;

    float avgScanTime = 0.1f; //TODO: must set this based on sample info!

    Adduct *adduct = nullptr;
    Classifier *clsf = nullptr;

    string adductName = ""; // only exists to assist in encoding/decoding adduct
    string clsfFile = ""; //only exists to assist in encoding/decoding classifier

    inline bool isIsotopes() {return (isC13Labeled || isN15Labeled || isS34Labeled || isD2Labeled);}

    string encodeParams();
    static IsotopeParameters decode(string encodedIsotopeParameters);

};

class PeakGroup {

	public:
        enum GroupType {None=0, C13Type=1, AdductType=2, CovariantType=3, IsotopeType=4, SRMTransitionType=5 };     //group types
        enum QType	   {AreaTop=0, Area=1, Height=2, AreaNotCorrected=3, RetentionTime=4, Quality=5, SNRatio=6, MS2Count=7 };
		PeakGroup();
		PeakGroup(const PeakGroup& o);
		PeakGroup& operator=(const PeakGroup& o);
		bool operator==(const PeakGroup* o);
        bool operator==(const PeakGroup& o);
        void copyObj(const PeakGroup& o);

		~PeakGroup();

		PeakGroup* parent;
        Compound* compound;
        Adduct*   adduct;

		vector<Peak> peaks;
        vector<PeakGroup> children;

        string srmId;
        string tagString;
        string searchTableName;

        //Issue 271: enables lookup of compound information when loaded compound dbs change
        string compoundId;
        string compoundDb;

        string displayName; //Issue 75: For use with tabledockwidget, other GUI displays

        //Issue 402: Saved isotopeParameters should only apply to bookmarked peaks
        IsotopeParameters isotopeParameters;

        /**
         * @brief importedCompoundName
         *
         * Introduced as a part of MAVEN Issue 146
         *
         * For use when a peak group should be connected to a compound, but should
         * be named something other than the linked compounds' name.
         *
         * The imported compound name should take precedence over the name linked compound.
         *
         * This is useful when an identification algorithm connects to a compound for display
         * of an original spectrum, but suggests summarizing the compound to a higher (less specific)
         * level (as in lipidomics).
         */
        string importedCompoundName;

        string getName();               //compound name + tagString + srmid

        vector<char> labels; //Issue 127: for use with a more intricate tagging system

        bool isFocused;

        int groupId;
        int metaGroupId;

        float maxIntensity;
        float meanRt;
        float meanMz;
        int  ms2EventCount;

        //isotopic information
        float expectedAbundance;
        int   isotopeC13count;
        
        //Flag for deletion
        bool deletedFlag;

        float minRt;
        float maxRt;
        float minMz;
        float maxMz;

        float blankMax;
        float blankMean;
        unsigned int blankSampleCount;

        int sampleCount;
        float sampleMean;
        float sampleMax;

        unsigned int maxNoNoiseObs;
        unsigned int  maxPeakOverlap;
        float maxQuality;
        float maxPeakFracionalArea;
        float maxSignalBaseRatio;
        float maxSignalBaselineRatio;
        int goodPeakCount;
        float expectedRtDiff;
        float groupRank;

        //for sample contrasts  ratio and pvalue
        float changeFoldRatio;
        float changePValue;

        //internal flags
        bool isComputedGroupStatistics;


        //LOSS
        int     chargeState;
        int     isotopicIndex;

        bool  	hasSrmId()  { return srmId.empty(); }
        void  	setSrmId(string id)	  { srmId=id; }
        inline  string getSrmId() { return srmId; }

        /**
         * Issue 368: SRM support
         */
        float srmPrecursorMz;

        /**
         * Issue 368: SRM support
         */
        float srmProductMz;

        bool isPrimaryGroup();
        inline bool hasCompoundLink()  { if(compound != NULL) return true ; return false; }
        inline bool isEmpty()   { if(peaks.size() == 0) return true; return false; }
		inline unsigned int peakCount()  { return peaks.size(); 	  }
		inline unsigned int childCount() { return children.size(); }
		inline Compound* getCompound() { return compound; }

		inline PeakGroup* getParent() { return parent; }

		inline vector<Peak>& getPeaks() { return peaks; } 
        inline vector<PeakGroup>& getChildren()  { return children; }

        vector<Scan*> getRepresentativeFullScans();
        vector<Scan*> getFragmentationEvents();
        void findHighestPurityMS2Pattern(float precPpmTolr);
        Scan* getAverageFragmentationScan(float productPpmTolr);

        void computeFragPattern(float productPpmTolr);
        void computeDIFragPattern(shared_ptr<DirectInfusionSearchParameters> directInfusionSearchParameters);
        void computePeaksSearchFragPattern(shared_ptr<PeaksSearchParameters> params);

        Peak* getHighestIntensityPeak();
        int getChargeStateFromMS1(float ppm);
        bool isMonoisotopic(float ppm);


		inline void setParent(PeakGroup* p) {parent=p;}
        inline float ppmDist(float cmass) { return mzUtils::ppmDist(cmass,meanMz); }

		inline void addPeak(const Peak& peak) { peaks.push_back(peak); peaks.back().groupNum=groupId; }
		inline void addChild(const PeakGroup& child) { children.push_back(child); children.back().parent = this;   }
		Peak* getPeak(mzSample* sample);

		GroupType _type;
		inline GroupType type() { return _type; }
		inline void setType(GroupType t)  { _type = t; }
        inline bool isIsotope() { return _type == IsotopeType; }
        inline bool isAdduct() {  return _type == AdductType; }

		void summary();
        void groupStatistics(bool isForceRecomputation=true);
		void updateQuality();
		float medianRt();
		float meanRtW();
        void pullIsotopes(IsotopeParameters isotopeParameters); //Issue 371: refactor to dedicated method

        bool isGroupGood();
        bool isGroupBad();
        bool isGroupLabeled(char label);
        void markGroupGood();
        void markGroupBad();
        void addLabel(char label);
        void toggleLabel(char label);
        string getPeakGroupLabel();

		void reduce();
		void fillInPeaks(const vector<EIC*>& eics);
		void computeAvgBlankArea(const vector<EIC*>& eics);
		void groupOverlapMatrix();

        Peak* getSamplePeak(mzSample* sample);
        FragmentationMatchScore fragMatchScore;
        Fragment fragmentationPattern;

		void deletePeaks();
        bool deletePeak(unsigned int index);

		void clear();
		void deleteChildren();
        bool deleteChild(unsigned int index);
        bool deleteChild(PeakGroup* child);
        void copyChildren(const PeakGroup& other);

		vector<float> getOrderedIntensityVector(vector<mzSample*>& samples, QType type);
		void reorderSamples();

		static bool compRt(const PeakGroup& a, const PeakGroup& b ) { return(a.meanRt < b.meanRt); }
		static bool compMz(const PeakGroup& a, const PeakGroup& b ) { return(a.meanMz > b.meanMz); }
		static bool compIntensity(const PeakGroup& a, const PeakGroup& b ) { return(a.maxIntensity > b.maxIntensity); }
		static bool compArea(const PeakGroup& a, const PeakGroup& b ) { return(a.maxPeakFracionalArea > b.maxPeakFracionalArea); }
		static bool compQuality(const PeakGroup& a, const PeakGroup& b ) { return(a.maxQuality > b.maxQuality); }
		static bool compRankPtr(const PeakGroup* a, const PeakGroup* b ) { return(a->groupRank < b->groupRank); }
		static bool compRatio(const PeakGroup& a, const PeakGroup& b ) { return(a.changeFoldRatio < b.changeFoldRatio); }
		static bool compPvalue(const PeakGroup* a, const PeakGroup* b ) { return(a->changePValue< b->changePValue); }
		static bool compC13(const PeakGroup* a, const PeakGroup* b) { return(a->isotopeC13count < b->isotopeC13count); }
		static bool compMetaGroup(const PeakGroup& a, const PeakGroup& b) { return(a.metaGroupId < b.metaGroupId); }
		//static bool operator< (const PeakGroup* b) { return this->maxIntensity < b->maxIntensity; }
        static bool compMaxIntensity(const PeakGroup* a, const PeakGroup* b) { return(a->maxIntensity > b->maxIntensity); }
        static void clusterGroups(vector<PeakGroup> &allgroups, vector<mzSample*>samples, double maxRtDiff, double minSampleCorrelation, double minPeakShapeCorrelation, double ppm, vector<double> mzDeltas=vector<double>());

private:
        void processLabel(char label, bool isToggle);
};

class Compound { 
        static MassCalculator* mcalc;

		private:
			PeakGroup _group;			//link to peak group
            bool      _groupUnlinked;
            float exactMass;

		public:
            Compound(string id, string name, string formula, int charge, float exactMass);
			Compound(string id, string name, string formula, int charge );
            virtual ~Compound(){} //empty destructor

            PeakGroup* getPeakGroup() { return &_group; }
			void setPeakGroup(const PeakGroup& group ) { _group = group; _group.compound = this; }
			bool hasGroup()    { if(_group.meanMz != 0 ) return true; return false; } 
			void clearGroup()  { _group.clear(); }
			void unlinkGroup() { _group.clear(); _groupUnlinked = true; }
            bool groupUnlinked() { return _groupUnlinked; }
            FragmentationMatchScore scoreCompoundHit(Fragment* f, float productPpmTolr, bool searchProton);

            int    cid;
            string id;
            string name;
            string formula;
            string smileString;
            string adductString;
            int    ionizationMode;

            string srmId;
            float expectedRt;
            int charge;

            //QQQ mapping
            string method_id;
            float precursorMz;	//QQQ parent ion
            float productMz;	//QQQ child ion
            float collisionEnergy; //QQQ collision energy
            float logP;

            string db;			//name of database for example KEGG, ECOCYC.. etc..
            bool  isDecoy;
            bool  virtualFragmentation;

            vector<float>fragment_mzs;
            vector<float>fragment_intensity;
            vector<string>fragment_labels;  //Added in Issue 159

            map<int,string>fragment_iontype;
            vector<string> category;

            float ajustedMass(int charge);
            static bool compMass(const Compound* a, const Compound* b )      { return(a->exactMass < b->exactMass);       }
            static bool compName(const Compound* a, const Compound* b )    { return(a->name < b->name);       }
            static bool compFormula(const Compound* a, const Compound* b ) { return(a->formula < b->formula); }

            inline float getExactMass() { return exactMass; }
            void setExactMass(float value) { exactMass = value; }

            inline void setFormula(string formulastr) { formula = formulastr; }
            string getFormula() { return formula; }
            map<string, string> metaDataMap = {};

            virtual vector<Compound*> getChildren();
            virtual vector<int> getConstituentMzs();
};

//Issue 416
enum SummarizedCompoundType{GENERAL=0, SUM_COMPOSITION=1, ACYL_CHAIN=2};

class SummarizedCompound : public Compound {

public:
    vector<Compound*> children;
    string summarizedName;
    SummarizedCompoundType type = SummarizedCompoundType::GENERAL;

    //Issue 314
    unordered_set<int> constituentMzsSet;

    SummarizedCompound(string summarizedName, vector<Compound*> childrenCompounds) : Compound("", summarizedName, "", 0){
        children = childrenCompounds;
        this->summarizedName = summarizedName;
        this->db = "summarized";
    }
    vector<Compound*> getChildren();

    //Issue 314
    vector<int> getConstituentMzs();

    static vector<pair<string, string>> parseCompoundId(string generalSummarizedCompoundName, bool debug=false);

    //Relies on children
    void computeSummarizedData();
    virtual ~SummarizedCompound(){}
};

class Ms3Compound : public Compound {
public:

    /**
     * @brief ms3_fragment_mzs, ms3_fragment_intensity, ms3_fragment_labels
     *
     * key is "int" value from mzUtils::mzToIntKey() and mzUtils::intKeyToMz()
     *
     */

    map<int, vector<float>> ms3_fragment_mzs;
    map<int, vector<float>> ms3_fragment_intensity;
    map<int, vector<string>> ms3_fragment_labels;

    Compound* baseCompound = nullptr;

    Ms3Compound(Compound* baseCompound) : Compound(string(baseCompound->id+"-ms3"), baseCompound->name, baseCompound->formula, baseCompound->charge, baseCompound->getExactMass()){
        this->baseCompound = baseCompound;
        computeMs3Spectra();
    }

    virtual ~Ms3Compound(){}

    //Relies on ms2 fragment labels
    void computeMs3Spectra();

    vector<Compound*> getChildren();
};

class Isotope {
public:
    string name;
    double mz;
    double abundance;
    int charge;
    int N15;
    int C13;
    int S34;
    int H2;

    Isotope(string name, float mass, int c=0, int n=0, int s=0, int h=0) {
        this->mz=mass; this->name=name; charge=0;
        C13=c; N15=n; S34=s; H2=h;
    }

    Isotope() {
        mz=0; abundance=0; N15=0; C13=0; S34=0; H2=0; charge=0;
    }

    Isotope(const Isotope& b) {
        name=b.name;
        mz=b.mz;
        abundance=b.abundance;
        N15=b.N15; S34=b.S34; C13=b.C13; H2=b.H2;
        charge = b.charge;
    }

    Isotope& operator=(const Isotope& b) {
        name=b.name; mz=b.mz; abundance=b.abundance;
        N15=b.N15; S34=b.S34; C13=b.C13; H2=b.H2;
        charge =b.charge;
        return *this;
    }

};

class Adduct { 

	public:
        Adduct(){ mass=0; charge=1; nmol=1; }

        Adduct(string name, float mass, int charge, int nmol){
            this->name=name; this->mass=mass; this->charge=charge; this->nmol=nmol;
        }

        string name;
		int	  nmol;
		float mass;
		float charge;
		
        //given adduct mass, compute monoisotopic parent mass
		inline float computeParentMass(float mz)  { return  (mz*abs(charge)-mass)/nmol; }

        //given parent monoisotopic mass, compute adduct mass
		inline float computeAdductMass(float pmz) { return (pmz*nmol+mass)/abs(charge); }

        static vector<Adduct*> loadAdducts(string filename);
};

struct AlignmentSegment {
		string sampleName;
		float  seg_start;
		float  seg_end;
		float  new_start;
		float  new_end;
		float updateRt(float oldRt);

};

/**
 * @brief The AnchorPoint struct
 * Container for storing information about rt value, whether or not the rt value was interpolated,
 * and the sample associated with the rt value.
 */
struct AnchorPoint {

        mzSample* sample;
        float rt;           //observedRt (AlignmentSegment.seg_start or AlignmentSegment.seg_end)
        bool isRtFromEIC;

        explicit AnchorPoint(mzSample* sample) {
            this->sample = sample;
        }

        inline void setInterpolatedRtValue(float rt){this->rt = rt; isRtFromEIC = false;}

        bool setEICRtValue(mzSlice* slice, int eic_smoothingWindow, float minPeakIntensity);
};

/**
 * @brief The AnchorPointSet class
 * a collection of individual AnchorPoints, each of which know if they were interpolated
 * across the RT range of the experiment, or if they were explicitly extracted.
 */
class AnchorPointSet {
public:

    //from constructor
    mzSlice* slice;
    vector<mzSample*> eicSamples;

    //used for last point in file
    AnchorPointSet() {
        slice = nullptr;
    }

    //used for reading from file.
    AnchorPointSet(double mzmin, double mzmax, double rtmin, double rtmax, int eic_smoothingWindowVal, float minPeakIntensityVal) {
        slice = new mzSlice(mzmin, mzmax, rtmin, rtmax);
        eic_smoothingWindow = eic_smoothingWindowVal;
        minPeakIntensity = minPeakIntensityVal;
    }

    AnchorPointSet(const PeakGroup& pg) {
        slice = new mzSlice(pg.minMz, pg.maxMz, pg.minRt, pg.maxRt);

        if (pg.peaks.size() >= minNumObservedSamples) {
            for (auto &p : pg.peaks) {
                eicSamples.push_back(p.sample);
            }
        } else {
            isValid = false;
        }

    }

    //computed by compute()
    map<mzSample*, AnchorPoint*> sampleToPoints{};

    /**
     * @brief compute
     * Method to determine sampleToPoints map.
     */
    void compute(const vector<mzSample*>& allSamples);

    /**
     * @brief minNumObservedSamples
     */
    void setEICSamplesByFilter(const vector<mzSample*>& allSamples, string stringFilter);

    int minNumObservedSamples = 2;
    bool isValid = true;

    //keep these attached as fields
    int eic_smoothingWindow = 5;
    float minPeakIntensity = 0.0f;

    static AnchorPointSet lastRt(vector<mzSample*>& allSamples);

    string toString();

};

class Aligner {

	public:
		Aligner();
		void doAlignment(vector<PeakGroup*>& peakgroups);
		vector<double> groupMeanRt();
		double checkFit();
        void Fit(int ideg);
		void PolyFit(int maxdegree);
		void saveFit();
		void restoreFit();
		void setMaxItterations(int x) { maxItterations=x; }
		void setPolymialDegree(int x) { polynomialDegree=x; }
		void setSamples(vector<mzSample*>set) { samples=set; }

		void loadAlignmentFile(string alignmentFile); //load alignment information from a file
		void doSegmentedAligment();	 //ralign scans using guided alignment

        inline void addSegment(string sampleName, AlignmentSegment* s) {
			alignmentSegments[sampleName].push_back(s); 
		}

        //AnchorPoint related updates
        vector<AnchorPointSet> groupsToAnchorPoints(vector<mzSample*>& samples, vector<PeakGroup*>& peakGroups, int eic_smoothingWindow, float minPeakIntensity);
        map<mzSample*, vector<pair<float, float>>> anchorPointSetToUpdatedRtMap(vector<AnchorPointSet>& anchorPoints, mzSample* refSample);

        void exportAlignmentFile(vector<AnchorPointSet>& anchorPoints, mzSample* refSample, string outputFile);


	private:
		vector< vector<float> > fit;
		vector<mzSample*> samples;
		vector<PeakGroup*> allgroups;
        map<string,vector<AlignmentSegment*> > alignmentSegments;
		int maxItterations;
		int polynomialDegree;
		
};


class ChargedSpecies{
public:
    ChargedSpecies(){
        deconvolutedMass=0; minZ=0; maxZ=0; countMatches=0; error=0; upCount=0; downCount=0; scan=NULL; totalIntensity=0; isotopicClusterId=0;
        minRt=maxRt=meanRt=0; istotopicParentFlag=false; minIsotopeMass=0; maxIsotopeMass=0; massDiffectCluster=-100; filterOutFlag=false;
        qscore=0;
    }

    float deconvolutedMass;     //parent mass guess
    float error;
    int minZ;
    int maxZ;
    int countMatches;
    float totalIntensity;
    float qscore;
    int upCount;
    int downCount;
    Scan* scan;
    vector<float> observedMzs;
    vector<float> observedIntensities;
    vector<float> observedCharges;

    int isotopicClusterId;
    int massDiffectCluster;
    bool istotopicParentFlag;
    bool filterOutFlag;

    float meanRt;
    float minRt;
    float maxRt;

    vector<ChargedSpecies*> isotopes;
    int minIsotopeMass;
    int maxIsotopeMass;

    static bool compIntensity(ChargedSpecies* a, ChargedSpecies* b ) { return a->totalIntensity > b->totalIntensity; }
    static bool compMass(ChargedSpecies* a, ChargedSpecies* b ) { return a->deconvolutedMass < b->deconvolutedMass; }
    static bool compMatches(ChargedSpecies* a, ChargedSpecies* b ) { return a->countMatches > b->countMatches; }
    static bool compRt(ChargedSpecies* a, ChargedSpecies* b ) { return a->meanRt < b->meanRt; }

    static bool compMetaGroup(ChargedSpecies* a, ChargedSpecies* b ) {
        return (a->isotopicClusterId * 100000 + a->deconvolutedMass  < b->isotopicClusterId* 100000 + b->deconvolutedMass);
    }

    void updateIsotopeStatistics() {
        minIsotopeMass = maxIsotopeMass =deconvolutedMass;
        for(unsigned int i=0; i < isotopes.size(); i++) {
            if (isotopes[i]->deconvolutedMass < minIsotopeMass) minIsotopeMass=isotopes[i]->deconvolutedMass;
            if (isotopes[i]->deconvolutedMass > maxIsotopeMass) maxIsotopeMass=isotopes[i]->deconvolutedMass;
        }
    }
};

/**
 * @brief The SearchParameters class
 *
 * base class for all parameters associated with spectral library / full dataset scans
 */
class SearchParameters {

    public:

    /** =======================
     * PROGRAM LEVEL
     * searchVersion: version of search protocol used to generate results.
     * ========================*/
     string searchVersion = "UNKNOWN";  //Issue 335

    /** =======================
     * SCAN FILTER ASSOCIATED
     * All parameters are arguments of Fragment::Fragment() constructor.
     * ========================*/

    float scanFilterMinFracIntensity = 0;
    float scanFilterMinSNRatio = 0;
    int scanFilterMaxNumberOfFragments = -1;
    int scanFilterBaseLinePercentile = 0;
    bool scanFilterIsRetainFragmentsAbovePrecursorMz = true;
    float scanFilterPrecursorPurityPpm = 0;
    float scanFilterMinIntensity = 0;

    //scan filter for MS1 scans
    float scanFilterMs1MinRt = -1.0f;
    float scanFilterMs1MaxRt = -1.0f;

    //scan filter for MS2 scans
    float scanFilterMs2MinRt = -1.0f;
    float scanFilterMs2MaxRt = -1.0f;

    /** =======================
     * CONSENSUS SPECTRUM ASSOCIATED
     * All parameters are arguments Fragment::buildConsensus() method.
     * ========================*/
    Fragment::ConsensusIntensityAgglomerationType consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
    bool consensusIsIntensityAvgByObserved = true;
    bool consensusIsNormalizeTo10K = false;

    //consensus spectrum formation of MS1 scans
    float consensusMs1PpmTolr = 10;
    int consensusMinNumMs1Scans = 0;
    float consensusMinFractionMs1Scans = 0;

    //consensus spectrum formation of MS2 scans
    float consensusPpmTolr = 10;
    int consensusMinNumMs2Scans = 0;
    float consensusMinFractionMs2Scans = 0;

    /** ===================
     * MS1 SEARCH RELATED
     * @param ms1PpmTolr: tolerance value used for matching theoretical ion m/z to an m/z peak in an MS1 scan
     * ==================== */
    float ms1PpmTolr = 5;

    /** ===================
     * MS2 SEARCH RELATED
     * @param ms2MinNumMatches: min number of reference peaks found in observed spectrum
     * @param ms2MinNumDiagnosticMatches: min number of reference peaks labeled as diagnostic found in observed spectrum
     * @param ms2MinNumUniqueMatches: min num of reference peaks that only match to a single compound.
     *          in a DIMS context, the library compounds may be pre-summarized.
     * @param ms2PpmTolr: m/z tolerance value used for matching reference <--> observed spectra in an MS2 scan
     *
     * Should only be used for comparing MS2 m/zs
     *
     * Contexts:
     *
     * [1] DirectInfusionMatchAssessment::computeMs2MatchAssessment():
     * > float maxDeltaMz = (params->ms2PpmTolr * static_cast<float>(t.precursorMz))/ 1000000;
     * > fragmentationMatchScore.ranks = Fragment::findFragPairsGreedyMz(&t, observedSpectrum, maxDeltaMz);
     *
     * [2] DirectInfusionMatchInformation::computeMs1PartitionFractions():
     * > float queryIntensity = scan->findClosestMzIntensity(queryMz, params->ms2PpmTolr); # scan is an MS2 scan, queryMz is an MS2 m/z
     *
     * [3] Fragment::scoreMatch():
     * >float maxDeltaMz = (productPpmTolr * static_cast<float>(precursorMz))/ 1000000; #variable is renamed to productPpmTolr
     *
     * @param ms2MinIntensity: minimum intensity value for an MS2 spectral peak to be considered real
     * ==================== */

    int ms2MinNumMatches = 5;
    int ms2MinNumDiagnosticMatches = 0; //Previously, diagnostic fragments always started with '*'
    int ms2MinNumUniqueMatches = 0;
    float ms2PpmTolr = 20;
    float ms2MinIntensity = 0;

    /** New Param 2020-05-04 **/
    map<string, int> ms2MinNumDiagnosticMatchesMap {};

    virtual string encodeParams() = 0;

    virtual ~SearchParameters() = default; //c++11 way
};

enum PeakGroupCompoundMatchingPolicy {
    ALL_MATCHES,
    SINGLE_TOP_HIT,
    TOP_SCORE_HITS
};

/**
 * @brief The PeaksSearchParameters class
 *
 * Used by Peaks Dialog in maven GUI
 */
class PeaksSearchParameters : public SearchParameters {

public:

    //baseline
    float baselineSmoothingWindow = 0;
    int baselineDropTopX = 80;

    //eic
    float eicSmoothingWindow = 0;
    string eicEicSmoothingAlgorithm = "GAUSSIAN";
    float eicMaxPeakGroupRtDiff = 0.5;
    int eicNoNoiseObs = 0;

    //quality
    float qualitySignalBaselineRatio = 1.00;
    float qualitySignalBlankRatio = 1.00;
    float qualityMinPeakGroupIntensity = 10000;
    int qualityMinPeakWidth = 5; //scans, no noise observations
    float qualityMinGoodPeakPerGroup = 0;
    float qualityMinPeakQuality = 0;
    string qualityClassifierModelName = "default.model";

    //isotopes
    bool isotopesIsRequireMonoisotopicPeaks = true;
    bool isotopesExtractIsotopicPeaks = false;
    float isotopesMzTolerance = 20.0f;

    //ms1 matching
    // ms1PpmTolr (SearchParameters)
    bool ms1IsMatchRtFlag = false;
    float ms1RtTolr = 2.0f;
    float ms1MassSliceMergePpm = 20.0f;
    float ms1MassSliceMergeNumScans = 10.0f;

    //ms2 matching
    bool ms2IsMatchMs2 = false;
    //ms2MinNumMatches (SearchParameters)
    //ms2MinNumDiagnosticMatches (SearchParameters)
    //ms2PpmTolr (SearchParameters)
    string ms2ScoringAlgorithm = "HyperGeomScore";
    float ms2MinScore = 0.0;

    //matching options
    bool matchingIsRequireAdductPrecursorMatch = true;
    bool matchingIsRetainUnknowns = false;
    string matchingLibraries = "";
    PeakGroupCompoundMatchingPolicy matchingPolicy = PeakGroupCompoundMatchingPolicy::SINGLE_TOP_HIT;

    //default constructor
    PeaksSearchParameters() {

        //ms1 matching
        ms1PpmTolr = 20;

        //scan filter parameters
        scanFilterMinFracIntensity = 0.01;
        scanFilterMinSNRatio = 1;
        scanFilterMaxNumberOfFragments = 1024;
        scanFilterBaseLinePercentile = 5;
        scanFilterIsRetainFragmentsAbovePrecursorMz = false;
        scanFilterPrecursorPurityPpm = 10;
        scanFilterMinIntensity = 0;

        //consensus spectrum formation parameters
        consensusPpmTolr = 10;
        consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
        consensusIsIntensityAvgByObserved = false;
        consensusIsNormalizeTo10K = true;
        consensusMinNumMs2Scans = 0;
        consensusMinFractionMs2Scans = 0.0f;
    }

public:
    string encodeParams();
    static shared_ptr<PeaksSearchParameters> decode(string encodedParams);
};


#endif
