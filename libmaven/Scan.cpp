#include <numeric>
#include "mzSample.h"

Scan::Scan(mzSample* sample, int scannum, int mslevel, float rt, float precursorMz, int polarity) {
    this->sample = sample;
    this->rt = rt;
    this->scannum = scannum;
    this->precursorMz = precursorMz;
    this->mslevel = mslevel;
    this->polarity = polarity;
    this->productMz=0;
    this->collisionEnergy=0;
    this->centroided=0;
    this->precursorCharge=0;
    this->precursorIntensity=0;
    this->isolationWindow=0;
    this->injectionTime=0;
    this->ms1PrecursorForMs3=0;

    /*if ( polarity != 1 && polarity != -1 ) {
        cerr << "Warning: polarity of scan is not 1 or -1 " << polarity << endl;
    }*/
}

void Scan::deepcopy(Scan* b) {
    this->sample = b->sample;
    this->rt = b->rt;
    this->scannum = b->scannum;
    this->precursorMz = b->precursorMz;
    this->precursorIntensity= b->precursorIntensity;
    this->precursorCharge= b->precursorCharge;
    this->mslevel = b->mslevel;
    this->polarity = b->polarity;
    this->productMz= b->productMz;
    this->collisionEnergy= b->collisionEnergy;
    this->centroided= b->centroided;
    this->intensity = b->intensity;
    this->mz    = b->mz;    
    this->scanType = b->scanType;
    this->filterLine = b->filterLine;
    this->filterString = b->filterString;
    this->setPolarity( b->getPolarity() );
    this->isolationWindow = b->isolationWindow;
    this->injectionTime = b->injectionTime;
    this->ms1PrecursorForMs3 = b->ms1PrecursorForMs3;
    this->lowerLimitMz = b->lowerLimitMz;
    this->upperLimitMz = b->upperLimitMz;
}

int Scan::findHighestIntensityPos(float _mz, float ppm) {
    return findHighestIntensityPos(_mz, _mz, ppm);
}

/**
 * @brief Scan::findHighestIntenstyPosDeltaMz
 * @param _mz
 * @param deltaMz
 * @return
 *
 * If the ppm must depend on a different value than _mz, include this ppmMz value.
 */
int Scan::findHighestIntensityPos(float _mz, float ppmMz, float ppm){
    float mzmin = _mz - ppmMz/1e6*ppm;
    float mzmax = _mz + ppmMz/1e6*ppm;

    vector<float>::iterator itr = lower_bound(mz.begin(), mz.end(), mzmin-0.1);
    int lb = itr-mz.begin();
    int bestPos=-1;  float highestIntensity=0;
    for(unsigned int k=lb; k < nobs(); k++ ) {
            if (mz[k] < mzmin) continue;
            if (mz[k] > mzmax) break;
            if (intensity[k] > highestIntensity ) {
                    highestIntensity=intensity[k];
                    bestPos=k;
            }
    }
    return bestPos;
}

//returns -1 if not found
//find highest intensity peak if multiple hits
int Scan::findHighestIntensityPosAMU(float queryMz, float tolr) {

    int pos = -1;

    //m/zs cannot possibly be found in the range
    if (queryMz < getMinMz() || queryMz > getMaxMz()) return pos;

    float minQueryMz = queryMz - tolr;
    float maxQueryMz = queryMz + tolr;

    auto lbQuery = lower_bound(mz.begin(), mz.end(), minQueryMz);
    auto lbQueryPos = static_cast<unsigned int>(lbQuery - mz.begin());

    float highestIntensity = -1;
    for (unsigned int i = lbQueryPos; i < mz.size(); i++) {
        if (mz[i] <= maxQueryMz) {
            if (intensity[i] > highestIntensity) {
                highestIntensity = intensity[i];
                pos = static_cast<int>(i);
            }
        } else {
            break;
        }
    }

    return pos;
}

//AMU Matching
int Scan::findClosestHighestIntensityPos(float _mz, float tolr) {
			float mzmin = _mz - tolr-0.001;
			float mzmax = _mz + tolr+0.001;

			vector<float>::iterator itr = lower_bound(mz.begin(), mz.end(), mzmin-0.1);
			int lb = itr-mz.begin();
			float highestIntensity=0; 
			for(unsigned int k=lb; k < mz.size(); k++ ) {
				if (mz[k] < mzmin) continue; 
				if (mz[k] > mzmax) break;
				if (intensity[k] > highestIntensity) highestIntensity=intensity[k];
			}
				
			int bestPos=-1; float bestScore=0;
			for(unsigned int k=lb; k < mz.size(); k++ ) {
				if (mz[k] < mzmin) continue; 
				if (mz[k] > mzmax) break;
				float deltaMz = (mz[k]-_mz); 
				float alignScore = sqrt(intensity[k] / highestIntensity)-(deltaMz*deltaMz);
			//	cerr << _mz << "\t" << k << "\t" << deltaMz << " " << alignScore << endl;
				if (bestScore < alignScore) { bestScore=alignScore; bestPos=k; }
			}
			//if(bestPos>=0) cerr << "best=" << bestPos << endl;
			return bestPos;
}


vector<int> Scan::findMatchingMzs(float mzmin, float mzmax) {

    vector<int> matches;
    auto itr = lower_bound(mz.begin(), mz.end(), mzmin);

    auto lb = static_cast<unsigned int>(itr - mz.begin());

    for(unsigned int k = lb; k < nobs(); k++ ) {
        if (mz[k] <= mzmax) {
            matches.push_back(static_cast<int>(k));
        } else {
            break;
        }
	}

//	cerr << "matches:" << mzmin << " " << mzmax << " " << matches.size() << endl;
	return matches;
}

//returns -1 if not found
//find closest m/z to queryMz if multiple hits
float Scan::findClosestMzIntensity(float queryMz, float ppm) {

    float queryIntensity = -1.0f;

    //m/zs cannot possibly be found in the range
    if (queryMz < getMinMz() || queryMz > getMaxMz()) return queryIntensity;

    float minQueryMz = queryMz - queryMz*ppm/1e6f;
    float maxQueryMz = queryMz + queryMz*ppm/1e6f;

    auto lbQuery = lower_bound(mz.begin(), mz.end(), minQueryMz);
    auto lbQueryPos = static_cast<unsigned int>(lbQuery - mz.begin());

    float mzDist = 99999999;
    for (unsigned int i = lbQueryPos; i < mz.size(); i++) {
        if (mz[i] <= maxQueryMz) {
            if (abs(mz[i]-queryMz) < mzDist) {
                mzDist = mz[i]-queryMz;
                queryIntensity = intensity[i];
            }
        } else {
            break;
        }
    }

    return queryIntensity;
}

//returns -1 if not found
//find closest m/z to queryMz if multiple hits
float Scan::findNormalizedIntensity(float queryMz, float standardMz, float ppm, float minScanIntensity){

    float normalizedIntensity = -1.0f;

    float queryIntensity = -1.0f;
    float standardIntensity = -1.0f;

    //m/zs cannot possibly be found in the range
    if (queryMz < getMinMz() || queryMz > getMaxMz() || standardMz < getMinMz() || standardMz > getMaxMz()) return normalizedIntensity;

    float minQueryMz = queryMz - queryMz*ppm/1e6f;
    float maxQueryMz = queryMz + queryMz*ppm/1e6f;

    auto lbQuery = lower_bound(mz.begin(), mz.end(), minQueryMz);
    auto lbQueryPos = static_cast<unsigned int>(lbQuery - mz.begin());

    float mzDist = 99999999;
    for (unsigned int i = lbQueryPos; i < mz.size(); i++) {
        if (mz[i] <= maxQueryMz) {
            if (abs(mz[i]-queryMz) < mzDist) {
                mzDist = mz[i]-queryMz;
                queryIntensity = intensity[i];
            }
        } else {
            break;
        }
    }

    // could not find query m/z with sufficient intensity in scan
    if (queryIntensity < 0 || queryIntensity < minScanIntensity) return normalizedIntensity;

    float minStandardMz = standardMz - standardMz*ppm/1e6f;
    float maxStandardMz = standardMz + standardMz*ppm/1e6f;

    auto lbStandard = lower_bound(mz.begin(), mz.end(), minStandardMz);
    auto lbStandardPos = static_cast<unsigned int>(lbStandard-mz.begin());

    mzDist = 99999999;
    for (unsigned int i = lbStandardPos; i < mz.size(); i++) {
        if (mz[i] <= maxStandardMz) {
            if (abs(mz[i]-standardMz) < mzDist) {
                mzDist = mz[i]-standardMz;
                standardIntensity = intensity[i];
            }
        } else {
            break;
        }
    }

    // could not find standard m/z with sufficient intensity in scan
    if (standardIntensity < 0 || standardIntensity < minScanIntensity) return normalizedIntensity;

    normalizedIntensity = queryIntensity/standardIntensity;

    return normalizedIntensity;
}

//removes intensities from scan that lower than X
void Scan::quantileFilter(int minQuantile) {
        if (intensity.size() == 0 ) return;
        if( minQuantile <= 0 || minQuantile >= 100 ) return;

        int vsize=intensity.size();
        vector<float>dist = quantileDistribution(this->intensity);
        vector<float>cMz;
        vector<float>cIntensity;
        for(int i=0; i<vsize; i++ ) {
            if ( intensity[i] > dist[ minQuantile ]) {
                cMz.push_back(mz[i]);
                cIntensity.push_back(intensity[i]);
            }
        }
        vector<float>(cMz).swap(cMz);
        vector<float>(cIntensity).swap(cIntensity);
        mz.swap(cMz);
        intensity.swap(cIntensity);
}


//removes intensities from scan that lower than X
void Scan::intensityFilter(int minIntensity) {
        if (intensity.size() == 0 ) return;

        //first pass.. find local maxima in intensity space
        int vsize=intensity.size();
        vector<float>cMz;
        vector<float>cIntensity;
        for(int i=0; i<vsize; i++ ) {
           if ( intensity[i] > minIntensity) { //local maxima
                cMz.push_back(mz[i]);
                cIntensity.push_back(intensity[i]);
            }
        }
        vector<float>(cMz).swap(cMz);
        vector<float>(cIntensity).swap(cIntensity);
        mz.swap(cMz);
        intensity.swap(cIntensity);
}

void Scan::simpleCentroid() { //centroid data

        if (intensity.size() < 5 ) return;

        //pass zero smooth..
		int smoothWindow = intensity.size() / 20;
		int order=2;

		if (smoothWindow < 1 )  { smoothWindow = 2; }
		if (smoothWindow > 10 ) { smoothWindow = 10; }

        mzUtils::SavGolSmoother smoother(smoothWindow,smoothWindow,order);
		//smooth once
        vector<float>spline = smoother.Smooth(intensity);
		//smooth twice
        spline = smoother.Smooth(spline);

        //find local maxima in intensity space
        int vsize=spline.size();
        vector<float>cMz;
        vector<float>cIntensity;
        for(int i=1; i<vsize-2; i++ ) {
            if ( spline[i] > spline[i-1] &&  spline[i] > spline[i+1] ) { //local maxima in spline space
					//local maxima in real intensity space
					float maxMz=mz[i]; float maxIntensity=intensity[i];
					for(int j=i-1; j<i+1; j++) {
							if (intensity[i] > maxIntensity) { maxIntensity=intensity[i]; maxMz=mz[i]; }
					}
                	cMz.push_back(maxMz);
                	cIntensity.push_back(maxIntensity);
            }
        }

        vector<float>(cMz).swap(cMz);
        vector<float>(cIntensity).swap(cIntensity);
        mz.swap(cMz);
        intensity.swap(cIntensity);

        centroided = true;
} 

bool Scan::hasMz(float _mz, float ppm) {
    float mzmin = _mz - _mz/1e6*ppm;
    float mzmax = _mz + _mz/1e6*ppm;
	vector<float>::iterator itr = lower_bound(mz.begin(), mz.end(), mzmin);
	//cerr << _mz  << " k=" << lb << "/" << mz.size() << " mzk=" << mz[lb] << endl;
	for(unsigned int k=itr-mz.begin(); k < nobs(); k++ ) {
        if (mz[k] >= mzmin && mz[k] <= mzmax )  return true;
		if (mz[k] > mzmax ) return false;
    } 
    return false;
}

float Scan::getTIC() {
    if (tic <= 0.0f) {
        tic = std::accumulate(intensity.begin(), intensity.end(), 0.0f);
    }
    return tic;
}

void Scan::summary() {

    cerr << "Polarity=" << getPolarity()
         << " msLevel="  << mslevel
         << " rt=" << rt
         << " m/z size=" << mz.size()
         << " ints size=" << intensity.size()
         << " precursorMz=" << precursorMz
         << " productMz=" << productMz
         << " srmID=" << filterLine
         << " totalIntensty=" << this->totalIntensity()
         << endl;

    for(int i=0; i<mz.size(); i++) {
        cerr << "\t" << mz[i] << "\t" << intensity[i] << endl;
    }

}

//generate multi charges series..endingin in change Zx,Mx
vector<float> Scan::chargeSeries(float Mx, unsigned int Zx) {
    //Mx  = observed m/z
    //Zz  = charge of Mx
    //n =  number of charge states to g
    vector<float>chargeStates(Zx+20,0);
    double M = (Mx*Zx)-Zx;
    for(unsigned int z=1; z<Zx+20; z++) chargeStates[z]=(M+z)/z;  
    return(chargeStates);
}


ChargedSpecies* Scan::deconvolute(float mzfocus, float noiseLevel, float ppmMerge, float minSigNoiseRatio, int minDeconvolutionCharge, int maxDeconvolutionCharge, int minDeconvolutionMass, int maxDeconvolutionMass, int minChargedStates ) {

    int mzfocus_pos = this->findHighestIntensityPos(mzfocus,ppmMerge);
    if (mzfocus_pos < 0 ) { cout << "ERROR: Can't find parent " << mzfocus << endl; return NULL; }
    float parentPeakIntensity=this->intensity[mzfocus_pos];
    float parentPeakSN=parentPeakIntensity/noiseLevel;
    if(parentPeakSN <=minSigNoiseRatio) return NULL;

    //cout << "Deconvolution of " << mzfocus << " pSN=" << parentPeakSN << endl;

    int scanTotalIntensity=0;
    for(unsigned int i=0; i<this->nobs();i++) scanTotalIntensity+=this->intensity[i];

    ChargedSpecies* x = new ChargedSpecies();

    for(int z=minDeconvolutionCharge; z <= maxDeconvolutionCharge; z++ ) {
        float expectedMass = (mzfocus*z)-z;     //predict what M ought to be
        int countMatches=0;
        float totalIntensity=0;
        int upCount=0;
        int downCount=0;
        int minZ=z;
        int maxZ=z;

        if (expectedMass >= maxDeconvolutionMass || expectedMass <= minDeconvolutionMass ) continue;

        bool lastMatched=false;
        float lastIntensity=parentPeakIntensity;
        for(int ii=z; ii < z+50 && ii<maxDeconvolutionCharge; ii++ ) {
            float brotherMz = (expectedMass+ii)/ii;
            int pos = this->findHighestIntensityPos(brotherMz,ppmMerge);
            float brotherIntensity = pos>=0?this->intensity[pos]:0;
            float snRatio = brotherIntensity/noiseLevel;
            if (brotherIntensity < 1.1*lastIntensity && snRatio > 2 && withinXppm(this->mz[pos]*ii-ii,expectedMass,ppmMerge)) {
                maxZ = ii;
                countMatches++;
                upCount++;
                totalIntensity += brotherIntensity;
                lastMatched=true;
                lastIntensity=brotherIntensity;
                //cout << "up.." << ii << " pos=" << pos << " snRa=" << snRatio << "\t"  << " T=" << totalIntensity <<  endl;
            } else if (lastMatched == true) {   //last charge matched ..but this one didn't..
                break;
            }
        }

        lastMatched = false;
              lastIntensity=parentPeakIntensity;
              for(int ii=z-1; ii > z-50 && ii>minDeconvolutionCharge; ii--) {
                  float brotherMz = (expectedMass+ii)/ii;
                  int pos = this->findHighestIntensityPos(brotherMz,ppmMerge);
                  float brotherIntensity = pos>=0?this->intensity[pos]:0;
                  float snRatio = brotherIntensity/noiseLevel;
                  if (brotherIntensity < 1.1*lastIntensity && snRatio > 2 && withinXppm(this->mz[pos]*ii-ii,expectedMass,ppmMerge)) {
                      minZ = ii;
                      countMatches++;
                      downCount++;
                      totalIntensity += brotherIntensity;
                      lastMatched=true;
                      lastIntensity=brotherIntensity;
                      //cout << "up.." << ii << " pos=" << pos << " snRa=" << snRatio << "\t"  << " T=" << totalIntensity <<  endl;
                  } else if (lastMatched == true) {   //last charge matched ..but this one didn't..
                      break;
                  }
              }

              if (x->totalIntensity < totalIntensity && countMatches>minChargedStates && upCount >= 2 && downCount >= 2 ) {
                      x->totalIntensity = totalIntensity;
                      x->countMatches=countMatches;
                      x->deconvolutedMass = (mzfocus*z)-z;
                      x->minZ = minZ;
                      x->maxZ = maxZ;
                      x->scan = this;
                      x->observedCharges.clear();
                      x->observedMzs.clear();
                      x->observedIntensities.clear();
                      x->upCount = upCount;
                      x->downCount = downCount;

                      float qscore=0;
                      for(int ii=minZ; ii <= maxZ; ii++ ) {
                              int pos = this->findHighestIntensityPos( (expectedMass+ii)/ii, ppmMerge );
                              if (pos > 0 ) {
                                      x->observedCharges.push_back(ii);
                                      x->observedMzs.push_back( this->mz[pos] );
                                      x->observedIntensities.push_back( this->intensity[pos] );
                                      float snRatio = this->intensity[pos]/noiseLevel;
                                      qscore += log(pow(0.97,(int)snRatio));
                              //      if(ii == z) cout << '*';
                              //      cout << setprecision(2) << snRatio << ",";
                              }
                      }
                      x->qscore = -20*qscore;
                      //cout << " upC=" << x->upCount << " downC=" << x->downCount << " qscore=" << -qscore <<  " M=" << x->deconvolutedMass << endl;
              }
          }
    // done..
    if ( x->countMatches > minChargedStates ) {
            float totalError=0; float totalIntensity=0;
            for(unsigned int i=0; i < x->observedCharges.size(); i++ ) {
                    float My = (x->observedMzs[i]*x->observedCharges[i]) - x->observedCharges[i];
                    float deltaM = abs(x->deconvolutedMass - My);
                    totalError += deltaM*deltaM;
                    totalIntensity += x->observedIntensities[i];
            }
            //cout << "\t" << mzfocus << " matches=" << x->countMatches << " totalInts=" << x->totalIntensity << " Score=" << x->qscore << endl;
            x->error = sqrt(totalError/x->countMatches);
            //cout << "-------- total Error= " << sqrt(totalError/x->countMatches) << " total Intensity=" << totalIntensity << endl;

            return x;
    } else {
            delete(x);
            x=NULL;
            return(x);





    }
}

/**
 * @brief Scan::getTopPeaks
 * @param minFracIntensity: retain peaks with intensity at or above this proportion of max intensity peak
 * @param minSNRatio: retain peaks with S:N at or above this value (use @param baseLineLevel for noise level)
 * @param baseLinePercentile: (expressed as a percentage) intensity percentile corresponding to S:N ratio of 1
 * @param minIntensity: minimum intensity value for a peak to be retained
 * @return
 */
vector <pair<float,float> > Scan::getTopPeaks(float minFracIntensity, float minSNRatio, int baseLinePercentile, float minIntensity) {
   unsigned int N = nobs();
   float baseline = -1;

    vector< pair<float,float> > selected;
	if(N == 0) return selected;

    vector<int> positions = this->intensityOrderDesc();	// [0]=highest-> [999]..lowest
	float maxI = intensity[positions[0]];				// intensity at highest position

   //compute baseline intensity --> used as proxy for noise level, SN ratio calculation
   //note that positions vector is sorted by descending intensity, hence 100-X
   if(baseLinePercentile>0 && baseLinePercentile<100) {
        float cutvalueF = (100.0f-static_cast<float>(baseLinePercentile))/100.0f;	// baseLinePercentile=5 -> cutValue 0.95 --> baseline=1
        unsigned int mid = static_cast<unsigned int>(N * cutvalueF);                // baseline position
        if(mid < N) baseline = intensity[positions[mid]];                           // intensity at baseline
   }

   for(unsigned int i=0; i<N; i++) {
           unsigned int pos = static_cast<unsigned int>(positions[i]);
           if (
                   intensity[pos] >= minIntensity &&
                   intensity[pos]/maxI >= minFracIntensity &&
                   (baseline < 0 || intensity[pos]/baseline >= minSNRatio)
               ) {
                   selected.push_back(make_pair(intensity[pos], mz[pos]));
           } else {
                   break;
           }
   }

   return selected;
}


vector<int> Scan::intensityOrderDesc() {
    vector<pair<float,int> > mzarray(nobs());
    vector<int>position(nobs());
    for(unsigned int pos=0; pos < nobs(); pos++ ) {
        mzarray[pos] = make_pair(intensity[pos],pos);
    }

    //Issue 195: Break ties based on position in intensity array
    sort(mzarray.begin(), mzarray.end(), [](const pair<float, int>& lhs, const pair<float, int>& rhs){
        if (lhs.first == rhs.first) {
            return lhs.second < rhs.second;
        } else {
            return lhs.first > rhs.first;
        }
    });

   //return positions in order from highest to lowest intensity
   for(unsigned int i=0; i < mzarray.size(); i++) { position[i] = mzarray[i].second; }
   return position;
}

string Scan::toMGF() { 
    std::stringstream buffer;
    buffer << "BEGIN IONS" << endl;
    if (sample) { buffer << "TITLE=" <<  sample->sampleName << "." << scannum << "." << scannum << "." << precursorCharge << endl; }
    buffer << "PEPMASS=" << setprecision(8) << precursorMz << " " << setprecision(3) << precursorIntensity << endl;
    buffer << "RTINSECONDS=" << setprecision(9) << rt*60 << "\n";
    buffer << "CHARGE=" << precursorCharge; if(polarity < 0) buffer << "-"; else buffer << "+"; buffer << endl;
    for(unsigned int i=0; i < mz.size(); i++) {
        buffer << setprecision(8) << mz[i] << " " << setprecision(3) << intensity[i] << endl;
    }
    buffer << "END IONS" << endl;
    //cout << buffer;
    return buffer.str();
}


vector<int> Scan::assignCharges(float ppmTolr) {
    if ( nobs() == 0) {
        vector<int>empty;
        return empty;
    }

    int N = nobs(); //current scan size
    vector<int>chargeStates (N,0);
    vector<int>peakClusters = vector<int>(N,0);
    vector<int>parentPeaks = vector<int>(N,0);
    int clusterNumber=0;

    //order intensities from high to low
    vector<int>intensityOrder = intensityOrderDesc();
    double NMASS=C13_MASS-12.00;

    //a little silly, required number of peaks in a series in already to call a charge
                          //z=0,   z=1,    z=2,   z=3,    z=4,   z=5,    z=6,     z=7,   z=8,
    int minSeriesSize[9] = { 1,     2,     3,      3,      3,     4,      4,       4,     5  } ;

    //for every position in a scan
    for(int i=0; i < N; i++ ) {
        int pos=intensityOrder[i];
        float centerMz = mz[pos];
        float centerInts = intensity[pos];
       // float ppm = (0.125/centerMz)*1e6;
        float ppm = ppmTolr*2;
       // cerr << pos << " " <<  centerMz << " " << centerInts << " " << clusterNumber << endl;
        if (chargeStates[pos] != 0) continue;  //charge already assigned

        //check for charged peak groups
        int bestZ=0; int maxSeriesIntenisty=0;
        vector<int>bestSeries;


        //determine most likely charge state
        for(int z=5; z>=1; z--) {
            float delta = NMASS/z;
            int zSeriesIntensity=centerInts;
            vector<int>series;


            for(int j=1; j<6; j++) { //forward
                float mz=centerMz+(j*delta);

                int matchedPos = findHighestIntensityPos(mz,ppm);
                if (matchedPos>0 && intensity[matchedPos]<centerInts) {
                    series.push_back(matchedPos);
                    zSeriesIntensity += intensity[matchedPos];
                } else break;
            }

            for(int j=1; j<3; j++) {  //back
                float mz=centerMz-(j*delta);
                int matchedPos = findHighestIntensityPos(mz,ppm);
                if (matchedPos>0 && intensity[matchedPos]<centerInts) {
                    series.push_back(matchedPos);
                    zSeriesIntensity += intensity[matchedPos];
                } else break;
            } 
            //cerr << endl;
            if (zSeriesIntensity>maxSeriesIntenisty) { bestZ=z; maxSeriesIntenisty=zSeriesIntensity; bestSeries=series; }
        }

        //if ( i < 50) cerr << centerMz << " " << bestZ << " " << bestSeries.size() << " " << minSeriesSize[bestZ] << endl;

        //series with highest intensity is taken to be be the right one
        if(bestZ > 0 and bestSeries.size() >= minSeriesSize[bestZ] ) {
            clusterNumber++;
            int parentPeakPos=pos;
            for(unsigned int j=0; j<bestSeries.size();j++) {
                int brother_pos =bestSeries[j];
                if(bestZ > 1 and mz[brother_pos] < mz[parentPeakPos]
                        and intensity[brother_pos] < intensity[parentPeakPos]
                        and intensity[brother_pos] > intensity[parentPeakPos]*0.25)
                        parentPeakPos=brother_pos;
                chargeStates[brother_pos]=bestZ;
                peakClusters[brother_pos]=clusterNumber;
             }

           //if ( i < 50 ) cerr << "c12parent: " << mz[parentPeakPos] << endl;
            peakClusters[parentPeakPos]=clusterNumber;
            parentPeaks[parentPeakPos]=bestZ;


            //cerr << "z-series: " <<  mz[pos] << endl;
            //cerr << "parentPeak=" << mz[parentPeakPos] << " " << bestZ << endl;
        }
    }
    return parentPeaks;
}

float Scan::baseMz() {
    float maxIntensity = 0;
    float baseMz = 0;
    for(unsigned int i=0; i < nobs(); i++ ) {
        if ( intensity[i] > maxIntensity) {
            maxIntensity = intensity[i];
            baseMz = mz[i];
        }
    }
    return baseMz;
}


bool Scan::isMonoisotopicPrecursor(float queryMz, float ppm, int charge) {
    const float C13_DELTA_MASS = 13.0033548378-12.0;

    //intensity o monoisotopic peak must at least 5% of the C13 peak (
    //Fix: this needs to be tunable, could be an issue for isotopic labeling experiments
    //if C12 completly disappers
    const float minParentIntensityFrac = 0.05;


    if (charge <= 0) return true;  //senity check default to true

    //look back, to see if there is a parent peak.
    int parentPos=this->findHighestIntensityPos(queryMz-C13_DELTA_MASS/charge,ppm);

    //found a potential C12 peak.
    if (parentPos != -1) {
        // compare intensities of presumbed parent and this peak, parent nees to be above some threshold
        int peakPos =  this->findHighestIntensityPos(queryMz,ppm);
        if (peakPos != -1 and intensity[parentPos] > intensity[peakPos]* minParentIntensityFrac) {
            return false;
        }
    }

    return true;
}


vector<Isotope> Scan::getIsotopicPattern(float centerMz, float ppm, int maxZ=6, int maxIsotopes=10) {

    //determine most likely charge state
    int bestZ=0;
    int maxSeriesIntenisty=0;
    vector<Isotope> isotopes;

    const double NMASS= 13.0033548378-12.0;

    int pos=this->findHighestIntensityPos(centerMz,ppm);
    if (pos<=0) return isotopes;
    float focusMz = this->mz[pos];
    float focusIntensity = this->intensity[pos];

    for(int z=1; z<maxZ; z++) {
        float delta = NMASS/z;
        double zSeriesIntensity=0;
        vector<int>pattern;
        int isoCount=0;

        for(int j=0; j<maxIsotopes; j++) { //forward
            float mz=focusMz+(j*delta);
            int matchedPos =  this->findHighestIntensityPos(mz,ppm);
            if (matchedPos>0) {
                float matchedInt =  this->intensity[matchedPos];
                if (matchedInt/focusIntensity<0.1) break;

                zSeriesIntensity += log(this->intensity[matchedPos]);
                pattern.push_back(matchedPos);
                isoCount++;
            } else break;
        }

        for(int j=1; j<maxIsotopes; j++) {  //back
            float mz=focusMz-(j*delta);
            int matchedPos = this->findHighestIntensityPos(mz,ppm);
            if (matchedPos>0) {
                float matchedInt =  this->intensity[matchedPos];
                if (matchedInt/focusIntensity>0.5) break;
                if (matchedInt/focusIntensity<0.1) break;

                zSeriesIntensity += log(this->intensity[matchedPos]);
                pattern.push_back(matchedPos);
                isoCount++;
            } else break;
        }

        if (zSeriesIntensity>maxSeriesIntenisty and isoCount>1) {
            bestZ=z;
            sort(pattern.begin(),pattern.end());
            isotopes.clear();;

            for(int i=0; i<pattern.size();i++) {
                int pos = pattern[i];
                Isotope iso("I",this->mz[pos],i,0,0,0);
                iso.charge=bestZ;
                iso.abundance=this->intensity[pos];
                isotopes.push_back(iso);;
            }
        }
    }
    return isotopes;
}

void Scan::log10Transform() {
    for(int i=0; i < this->intensity.size(); i++) {
        if (this->intensity[i] > 0) {
            this->intensity[i] = log10(this->intensity[i]);
        }
    }
}


Scan* Scan::getLastFullScan(int historySize=50) {
	if (!this->sample) return 0;
    int scanNum = this->scannum;
    for(int i=scanNum; i>scanNum-historySize; i--) {
        Scan* lscan = this->sample->getScan(i);
        if (!lscan or lscan->mslevel > 1) continue;
		return lscan; // found ms1 scan, all is good
	}
	return 0;
}

string Scan::getSignature(int limitSize) {

    stringstream SIG;
    map<int,bool>seen;

    int mz_count=0;
    for(int pos: intensityOrderDesc() ) {
        float mzround = (int) mz[pos];
        int peakIntensity = (int) intensity[pos];

        if(! seen.count(mzround)) {
            SIG << "[" << setprecision(9) << mz[pos] << "," << peakIntensity << "]";
            seen[mzround]=true;
        }

        if (mz_count++ >= limitSize) break;
    }

    return SIG.str();
}

vector<mzPoint> Scan::getIsolatedRegion(float isolationWindowAmu=1.0) {

	vector<mzPoint>isolatedSegment;
	if(! this->sample) return isolatedSegment;

	//find last ms1 scan or get out
	Scan* lastFullScan = this->getLastFullScan(50);
	if (!lastFullScan) return isolatedSegment;

	//no precursor information
	if (this->precursorMz <= 0) return isolatedSegment;

	//extract isolated region 
	float minMz = this->precursorMz-isolationWindowAmu/2.0;
	float maxMz = this->precursorMz+isolationWindowAmu/2.0;
	
	for(int i=0; i< lastFullScan->nobs(); i++ ) {
		if (lastFullScan->mz[i] < minMz) continue;
		if (lastFullScan->mz[i] > maxMz) break;
		isolatedSegment.push_back(mzPoint(lastFullScan->mz[i], lastFullScan->intensity[i]));
	}
	return isolatedSegment;
}


double Scan::getPrecursorPurity(float ppm=10.0) {
	//did not locate precursor mass
    if (this->precursorMz <= 0 ) return 0;
    if (this->sample == 0 ) return 0;

	//extract isolated window
    vector<mzPoint>isolatedSegment = this->getIsolatedRegion(this->isolationWindow);
	if (isolatedSegment.size() == 0) return 0;

	//get last full scan
	Scan* lastFullScan = this->getLastFullScan();
	if (!lastFullScan) return 0;

	//locate intensity of isoloated mass
    int pos = lastFullScan->findHighestIntensityPos(this->precursorMz,ppm);
	if (pos < 0) return 0;
	double targetInt = lastFullScan->intensity[pos];

	//calculate total intensity in isolted segment
	double totalInt=0;
	for (mzPoint& x: isolatedSegment) totalInt += x.intensity();

	if (totalInt>0) {
		return targetInt/totalInt;
	} else {
		return 0;
	}
}

TMT Scan::tmtQuant() {
    float TMT_ION_MASS[] = { 126.127726,127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145,131.138180 };

    std::vector<double>tmtReporterIons(10,0); 

    //signal
    double totalSignal=0;
    int    tmtTagCount=0;
    for(int i=0; i<10; i++ ) {
        int pos = this->findHighestIntensityPos(TMT_ION_MASS[i],20.00);

        if (pos>0)  {
            tmtReporterIons[i] = this->intensity[pos]; 
            totalSignal+= this->intensity[pos]; 
            tmtTagCount++;
        }
    }

    TMT tmtquant;
    tmtquant.scannum = this->scannum;
    tmtquant.tmtTags = tmtTagCount;
    tmtquant.tmtTotalIntensity = totalSignal;
    tmtquant.tmtIons = tmtReporterIons;

    //calculate noise
    double totalI=0; int n=0;
    for(int i=0; i< nobs(); i++ ) {
	if (mz[i] < 126) continue;
	if (mz[i] > 131.2) break;
	totalI += intensity[i];
	n++; 
    }

    //remove signal contribution
    totalI = totalI-totalSignal;
    n = n - tmtTagCount;

    tmtquant.noise =  0;
    if(n>0 and totalI >0) tmtquant.noise =  totalI/n;
    //if(tmtquant.noise < MIN_TMT_ION) tmtquant.noise=MIN_TMT_ION;

    /*
    cerr << "TMT SCAN:" << endl;
    for(int i=0; i<this->nobs(); i++ ) {
	if (this->mz[i] < 126) continue;
	if (this->mz[i] > 131.2) break;
        cerr << setprecision(6) << this->mz[i] << "\t" << setprecision(3) << this->intensity[i] << endl;
    }
    */
    return tmtquant;
}

float Scan::getMinMz(){
    if (lowerLimitMz > 0) return lowerLimitMz;
    if (!mz.empty()) return mz[0];
    return -1.0f;
}

float Scan::getMaxMz(){
    if (upperLimitMz > 0) return upperLimitMz;
    if (!mz.empty()) return mz[mz.size()-1];
    return -1.0f;
}

/**
 * @brief Scan::snapToGrid
 * Update the m/z and intensity values based on a grid that starts at a minimum m/z value of
 * Scan::getMinMz() (usually metadata collected in the mzML file), and continues with spacing
 * provided by params->binMzWidth.
 *
 * Note that this function mutates the m/z values in the Scan object, so should be used with caution.
 *
 * @param params
 * @param debug
 */
void Scan::snapToGrid(shared_ptr<ScanParameters> params, bool debug) {
    if (params->binMzWidth < 0) return;

    //snapping to grid is illegal if the Scan has already been snapped to grid.
    if (this->snappedToGridSize > 0) return;

    double minMz = getMinMz(); // instrument setting, if available.

    double currentGrid = minMz;
    double nextGrid = minMz + params->binMzWidth;

    //  <mz,   intensity>
    map<double, vector<float>> snappedMzAndIntensity{};

    //organize intensities into mz bins, where the mz bins correspond to grid points.
    for (unsigned int i = 0; i < mz.size(); i++) {
        double mzVal = this->mz[i];

        // Adjust currentGrid and nextGrid such that currentGrid is less than (or equal to) mz,
        // and nextGrid is greater than mz.
        // Because the m/z list is sorted, this only manifests as mz being greater than nextGrid.
        if (mzVal > nextGrid) {
            double delta = mzVal - currentGrid;
            int numChunks = floor(delta/params->binMzWidth);

            currentGrid = currentGrid + (numChunks * params->binMzWidth);
            nextGrid = nextGrid + (numChunks * params->binMzWidth);

            if (debug) {
                cout << "i=" << i
                     << ": mzVal=" << mzVal
                     << ", delta=" << delta
                     << "--> currentGrid=" << currentGrid
                     << " nextGrid=" << nextGrid
                     << endl;
            }
        }

        //little adjustments, in case previous step wasn't quite enough
        while(mzVal > nextGrid) {
            nextGrid = nextGrid + params->binMzWidth;
            currentGrid = currentGrid + params->binMzWidth;
        }

        double deltaToCurrent = mzVal - currentGrid;
        double deltaToNext = nextGrid - mzVal;

        double snapMz = 0.0f;
        if (deltaToCurrent <= deltaToNext) {
            snapMz = currentGrid;
        } else {
            snapMz = nextGrid;
        }

        if (debug) {
            cout << currentGrid << "  <  " << mzVal << "   <  " << nextGrid
                 << " --> snapToMz=" << snapMz
                 << endl;
        }

        if (snappedMzAndIntensity.find(snapMz) == snappedMzAndIntensity.end()) {
            snappedMzAndIntensity.insert(make_pair(snapMz, vector<float>{}));
            if (debug) {
                cout << "Adding mz=" << snapMz << endl;
            }
        }
        snappedMzAndIntensity[snapMz].push_back(this->intensity[i]);
    }

    // reduction step
    vector<float> snappedMz(snappedMzAndIntensity.size());
    vector<float> snappedIntensity(snappedMzAndIntensity.size());

    unsigned int counter = 0;

    for (auto it = snappedMzAndIntensity.begin(); it != snappedMzAndIntensity.end(); ++it) {

        float intensity = 0.0f;
        if (params->binIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
            intensity = median(it->second);
        } else if (params->binIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
            intensity = accumulate(it->second.begin(), it->second.end(), 0.0f) / it->second.size();
        } else if (params->binIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Sum) {
            intensity = accumulate(it->second.begin(), it->second.end(), 0.0f);
        } else if (params->binIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Max) {
            intensity = *max_element(it->second.begin(), it->second.end());
        }

        snappedMz[counter] = it->first;
        snappedIntensity[counter] = intensity;

        if (debug) {
            cout << "(" << it->first << ", " << intensity << ")" << endl;
        }

        counter++;
    }

    this->mz = snappedMz;
    this->intensity = snappedIntensity;

    this->snappedToGridSize = params->binMzWidth;
}

string Scan::getBinaryEncodedData(bool includeTags) {
    //base64::

    //TODO

    return "";
}
