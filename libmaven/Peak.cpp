#include "mzSample.h"

Peak::Peak() {
	pos = 0;
    eic = nullptr;
    sample = nullptr;
}

Peak::Peak(EIC* e, int p) { 
        pos = static_cast<unsigned int>(p);
		eic = e;
		baseMz=0;
		fromBlankSample=false;
		groupOverlap=0;
        groupNum=0;
		groupOverlapFrac=0;
		label=0;
		localMaxFlag=0;
		maxpos=pos; 
		maxscan=0;
		medianMz=0;
		minpos=pos; 
		minscan=0;
		mzmax=0;
		mzmin=0;
		noNoiseObs=0;  
		noNoiseFraction=0;
		peakArea=0; 
		peakAreaCorrected=0;
		peakAreaFractional=0;
		peakAreaTop=0;
		peakIntensity=0; 
		peakMz=0;
		peakBaseLineLevel=0;
		quality=0;
		rt=0;
		rtmax=0;
		rtmin=0;
        sample = nullptr; //possibly reassigned near end of constructor
		scan=0;
		signalBaselineRatio=0;
		symmetry=0;
		width=0;
		gaussFitSigma=10;
		gaussFitR2=FLT_MAX;
		peakRank=INT_MAX;

		chargeState=0;
		isMonoIsotopic=false;
		selectionScore=0;
        ms2EventCount=0;

        if (eic and eic->sample) sample = eic->sample;

        //Issue 549: new fields
        smoothedIntensity = 0.0f;
        smoothedPeakArea = 0.0f;
        smoothedPeakAreaCorrected = 0.0f;
        smoothedPeakAreaTop = 0.0f;
        smoothedSignalBaselineRatio = 0.0f;

        minPosFWHM = 0;
        maxPosFWHM = 0;
        minScanFWHM = 0;
        maxScanFWHM = 0;
        rtminFWHM = 0.0f;
        rtmaxFWHM = 0.0f;
        peakAreaFWHM = 0.0f;
        smoothedPeakAreaFWHM = 0.0f;
}


void Peak::copyObj(const Peak& o ) {
		pos= o.pos;
		eic = o.eic;
		sample = o.sample;
		baseMz= o.baseMz;
		fromBlankSample= o.fromBlankSample;
		groupOverlap= o.groupOverlap;
		groupOverlapFrac= o.groupOverlapFrac;
        groupNum=o.groupNum;
		label= o.label;
		localMaxFlag= o.localMaxFlag;
		maxpos= o.maxpos;
		maxscan= o.maxscan;
		medianMz= o.medianMz;
		minpos= o.minpos;
		minscan= o.minscan;
		mzmax= o.mzmax;
		mzmin= o.mzmin;
		noNoiseObs= o.noNoiseObs;
		noNoiseFraction= o.noNoiseFraction;
		peakArea= o.peakArea;
		peakAreaCorrected=o.peakAreaCorrected;
		peakAreaFractional=o.peakAreaFractional;
		peakAreaTop=o.peakAreaTop;
		peakIntensity=o.peakIntensity;
		peakMz= o.peakMz;
		quality= o.quality;
		peakBaseLineLevel=o.peakBaseLineLevel;
		rt=o.rt;
		rtmax= o.rtmax;
		rtmin= o.rtmin;
		scan= o.scan;
		signalBaselineRatio= o.signalBaselineRatio;
		symmetry= o.symmetry;
		width= o.width;
		gaussFitSigma=o.gaussFitSigma;
		gaussFitR2= o.gaussFitR2;
		peakRank=o.peakRank;

		chargeState=o.chargeState;
		isMonoIsotopic=o.isMonoIsotopic;
		ms2EventCount=o.ms2EventCount;
		selectionScore=o.selectionScore;

        //Issue 549: new fields
        smoothedIntensity = o.smoothedIntensity;
        smoothedPeakArea = o.smoothedPeakArea;
        smoothedPeakAreaCorrected = o.smoothedPeakAreaCorrected;
        smoothedPeakAreaTop = o.smoothedPeakAreaTop;
        smoothedSignalBaselineRatio = o.smoothedSignalBaselineRatio;

        minPosFWHM = o.minPosFWHM;
        maxPosFWHM = o.maxPosFWHM;
        minScanFWHM = o.minScanFWHM;
        maxScanFWHM = o.maxScanFWHM;
        rtminFWHM = o.rtminFWHM;
        rtmaxFWHM = o.rtmaxFWHM;
        peakAreaFWHM = o.peakAreaFWHM;
        smoothedPeakAreaFWHM = o.smoothedPeakAreaFWHM;
}

Peak& Peak::operator=(const Peak& o)  {
	copyObj(o);
	return *this;
}

Peak::Peak(const Peak& o ) { copyObj(o); }

//test
int Peak::getChargeState() {
	Scan* s = this->getScan();
	int peakPos = -1;
	if (s) peakPos = s->findHighestIntensityPos(peakMz,0.5);

	if (s and peakPos > 0 and peakPos < s->nobs() ) {
        float ppm = peakMz/1e6*10;
		vector<int> parentPeaks = s->assignCharges(ppm);
		cerr << "pos=" << peakPos << endl; //peakPos
		cerr << "peakMz=" << peakMz << endl;
		cerr << "peakMzscan=" << s->mz[peakPos] << endl;
		cerr << "chargeState=" << parentPeaks[peakPos] << endl;
        this->chargeState=parentPeaks[peakPos];

        return parentPeaks[peakPos];
	} else {
		return 0;
	}
}

vector<Scan*> Peak::getFragmentationEvents(float ppmWindow) {

    vector<Scan*> matchedscans;

    if (!sample) return matchedscans;

	float amuTol = peakMz/1e6*ppmWindow;

    double minMz = peakMz - amuTol;
    double maxMz = peakMz + amuTol;

    if (mzmax - mzmin > 0.5f) { //direct infusion peak
        minMz = mzmin;
        maxMz = mzmax;
    }

    for( unsigned int j=0; j < sample->scans.size(); j++ ) {
         	Scan* scan = sample->scans[j];
            if (scan->mslevel != 2) continue; //ms2 scans only
       		if (scan->rt < rtmin) continue;
            if (scan->rt > rtmax) break;
            if( scan->precursorMz >= minMz and scan->precursorMz <= maxMz) {
                matchedscans.push_back(scan);
            }
    }
    return matchedscans;
}


Fragment Peak::getConsensusFragmentation(float ppmWindow,float productPpmTolr=20) {
            vector<Scan*> ms2events =  this->getFragmentationEvents(ppmWindow);
            if (ms2events.size()==0) return Fragment();
            Fragment* f = new Fragment(ms2events[0],0.0001,0,1024);
            for(Scan* s: ms2events) f->addFragment(new Fragment(s,0.0001,0,1024));

            f->buildConsensus(productPpmTolr);
            Fragment fcons(f->consensus);
            delete(f);
            return fcons;
}



vector<mzLink> Peak::findCovariants() {

		vector<mzLink>covariants;
		if (sample == NULL ) return covariants;

		//find scan range in which we will be checking for covariance
		vector<Scan*> scans;
		for(unsigned int i=0; i<sample->scans.size(); i++ ) {
				Scan* _scan = sample->scans[i];
				if(_scan == NULL ) continue;
				if(_scan->mslevel != 1) continue;
				if(_scan->rt < rt-0.1) continue;
				if(_scan->rt > rt+0.1) break;
				scans.push_back(_scan);
		}
		if ( scans.size() == 0 ) return covariants;
		int scanCount= scans.size();
		cerr << "Scans=" << scanCount << endl;

		//constuct Map
		map<int, vector<float> >M;
		map<int, vector<float> >::iterator itr;
		for(unsigned int i=0; i<scans.size(); i++ ) {
			Scan* _scan = scans[i];
			for(unsigned int j=0; j<_scan->nobs(); j++ ) {
					int rmz = int(_scan->mz[j]*1000);
					if (M[rmz].size()==0 )  M[rmz].resize(scanCount);
					M[rmz][i]= _scan->intensity[j];
			}
		} 

		//merge ajacent slices
		int lastMz=0;
		for(itr = M.begin(); itr != M.end(); ++itr) {
				float rmz = (*itr).first;
				if ( lastMz != 0 && rmz-lastMz < 2 ) {//merge
						for (int j=0; j < scanCount; j++ ) { 
								M[rmz][j] = max( M[lastMz][j], M[rmz][j]);
								M[lastMz][j]=0;
						}
				}
				lastMz = (*itr).first;
		}

		//refvector
		int refbinguess = int(peakMz*1000);
		int maxobs=0;
		vector<float> yref;

		for (int i=refbinguess-5; i < refbinguess+5; i++ ) {
				if ( M.count(i) ) {
						vector<float> y =  M[i];
						int nonzerCount=0; 
						for (int j=0; j<scanCount; j++ ) if (y[j]!=0) nonzerCount++;
						if ( nonzerCount > maxobs ) { yref=y; maxobs=nonzerCount; }
				}
		}
		if (yref.size()==0) return covariants;

		cerr << "Reference" << endl;
        for (int j=0; j<scanCount; j++ ) {
            cerr << yref[j];
        }
        cerr << endl;

		for(itr = M.begin(); itr != M.end(); ++itr) { 
				int rmz = (*itr).first;
				vector<float> y =  (*itr).second;
				//float score=matchScore(yref, y );
				float score = mzUtils::correlation(yref,y);
				if ( (float) score < 0.5) continue;
				mzLink link;
				link.mz1 = peakMz;
				link.mz2 = rmz/1000.0;
				link.value1 = pos;
				link.correlation = score;
				link.note = "Covariant";
				covariants.push_back(link);
		}
		return covariants;
}
