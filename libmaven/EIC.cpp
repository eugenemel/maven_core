#include "mzSample.h"

EIC::~EIC() {  
	peaks.clear();
}

EIC* EIC::clone() {
	EIC* clonedEIC = new EIC();

	clonedEIC->rtmin = rtmin;
	clonedEIC->rtmax = rtmax;
	clonedEIC->mzmin = mzmin;
	clonedEIC->mzmax = mzmax;

	clonedEIC->sampleName = sampleName;
	clonedEIC->sample = sample;
	
    for(unsigned int i=0; i<4;i++) clonedEIC->color[i]=color[i];

	clonedEIC->maxIntensity = maxIntensity;
	clonedEIC->totalIntensity = totalIntensity;
	clonedEIC->eic_noNoiseObs = eic_noNoiseObs;

    clonedEIC->smootherType= smootherType;
	clonedEIC->baselineSmoothingWindow=baselineSmoothingWindow;  
    clonedEIC->baselineDropTopX=baselineDropTopX;
	clonedEIC->peaks = peaks;


	clonedEIC->rt = rt;
	clonedEIC->intensity = intensity;
	clonedEIC->mz = mz;
	clonedEIC->scannum = scannum;
    clonedEIC->spline = spline;
    clonedEIC->baseline = baseline;

	return clonedEIC;
}

EIC* EIC::eicMerge(const vector<EIC*>& eics) {

	EIC* meic = new EIC();

	unsigned int maxlen = 0;	
	float minRt = DBL_MAX;
	float maxRt = DBL_MIN;
	for (unsigned int i=0; i < eics.size(); i++ )  {
		if ( eics[i]->size() > maxlen )  maxlen = eics[i]->size() + 1; 
		if ( eics[i]->rtmin  < minRt  )  minRt = eics[i]->rtmin;
		if ( eics[i]->rtmax  > maxRt  )  maxRt = eics[i]->rtmax;
	}
	
	if (maxlen == 0 ) return meic;

	//create new EIC
    meic->sample = nullptr;
	vector<float>intensity(maxlen,0);
	vector<float>rt(maxlen,0);
	vector<int>scans(maxlen,0);
	vector<float>mz(maxlen,0);
    vector<int>mzcount(maxlen,0);

	//smoothing 	//initalize time array
	for (unsigned int i=0; i < maxlen; i++ ) { 
        rt[i] = minRt + i*((maxRt-minRt)/maxlen); 
		scans[i]=i; 
    } 

	//combine intensity data from all pulled eics
	for (unsigned int i=0; i< eics.size(); i++ ) {
            EIC* e = eics[i];
            for (unsigned int j=0; j< e->size(); j++ ) {
                unsigned int bin = ((e->rt[j] - minRt ) / (maxRt-minRt) * maxlen);
                if (bin >= maxlen) bin=maxlen-1; 


                if(e->spline.size() and e->spline[j] > 0) {
                    intensity[bin] += e->spline[j];
                } else {
                    intensity[bin] += e->intensity[j];
				}

                if(e->mz[j] > 0) {
                    mz[bin] += e->mz[j];
                    mzcount[bin]++;
                }
            }
	}

   unsigned int eicCount=eics.size();
   for(unsigned int i=0; i<maxlen;i++ ) {
        intensity[i] /= eicCount;
         if( intensity[i] > meic->maxIntensity) meic->maxIntensity=intensity[i];
         if (mzcount[i]) mz[i] /=  mzcount[i];
         meic->totalIntensity += intensity[i];
	}

	//copy to new EIC
    meic->rtmin =  minRt;
	meic->rtmax =  maxRt;
	meic->intensity = intensity;
	meic->rt = rt;
	meic->scannum = scans;
	meic->mz = mz;
	meic->sampleName = eics[0]->sampleName;
	meic->sample 	 = eics[0]->sample;
	return meic;
}

void EIC::computeBaseLine(int smoothing_window, int dropTopX) {

    if (baseline.size()) {  //delete previous baseline if exists
        baseline.clear();
        eic_noNoiseObs=0;
        baselineQCutVal=0.0f;
    }

    unsigned int n = static_cast<unsigned int>(intensity.size());
	if (n == 0)  return;

	try { 
		baseline = vector<float>(n,0);
	}  catch(...) { 
        cout << "Exception caught while allocating memory " << n << "floats " << endl;
	}

    //Issue 572: Dropping no values amounts to assuming no measured intensities
    //are below the baseline.
    if (dropTopX < 0 || dropTopX > 100) {
        baselineQCutVal = -1.0f; // all intensity values are greater, everything matches
        baseline = intensity;
    } else {
        //sort intensity vector
        vector<float> tmpv = intensity;
        std::sort(tmpv.begin(),tmpv.end());

        //compute maximum intensity of baseline, any point above this value will be dropped
        //user specifies quantile of points to keep, for example
        //drop 60% of highest intensities = cut at 40% value;
        float cutvalueF = (100.0-(float) dropTopX)/101;
        unsigned int pos = tmpv.size() * cutvalueF;
        float qcut=0;
        pos < tmpv.size() ? qcut = tmpv[pos] : qcut = tmpv.back();

        baselineQCutVal = qcut;

        //drop all points above maximum baseline value
        for(unsigned int i=0; i<n; i++ ) {
            if (intensity[i] > qcut) {
                baseline[i]=qcut;
            } else {
                baseline[i] = intensity[i];
            }
        }
    }

    //smooth baseline
    mzUtils::GaussianSmoother smoother(static_cast<unsigned int>(smoothing_window));
    vector<float> smoothed = smoother.smooth(baseline);
    for (unsigned int i = 0; i < smoothed.size(); i++) baseline[i] = smoothed.at(i);

    //count number of observation in EIC above baseline
    for(unsigned int i=0; i < n; i++) {
        if(intensity[i]>baseline[i]) eic_noNoiseObs++;
    }
}

//Issue 549
void EIC::computeBaselineByNonPeakIntensity(shared_ptr<PeakPickingAndGroupingParameters> params, bool debug) {

    if (params->eicBaselineEstimationType != EICBaselineEstimationType::EIC_NON_PEAK_MAX_SMOOTHED_INTENSITY &&
            params->eicBaselineEstimationType != EICBaselineEstimationType::EIC_NON_PEAK_MEDIAN_SMOOTHED_INTENSITY) {
        return;
    }

    if (peaks.empty()) return;

    unsigned long currentPos = 0;

    vector<float> nonPeakIntensities{};

    for (auto peak : peaks) {

        unsigned long peakMinPos = peak.minpos;
        unsigned long peakMaxPos = peak.maxpos;

        if (debug) cout << "peakMinPos: " << peakMinPos  << endl;

        while (currentPos < peakMinPos) {
            if (debug) cout << currentPos << ": " << spline[currentPos] << endl;
            nonPeakIntensities.push_back(spline[currentPos]);
            currentPos++;
        }

        if (debug) cout << "peakMaxPos: " << peakMaxPos  << endl;

        //skip over any positions associated with a peak.
        currentPos = peakMaxPos + 1;
    }

    while(currentPos < spline.size()) {
        if (debug) cout << currentPos << ": " << spline[currentPos] << endl;
        nonPeakIntensities.push_back(spline[currentPos]);
        currentPos++;
    }

    if (params->eicBaselineEstimationType == EICBaselineEstimationType::EIC_NON_PEAK_MAX_SMOOTHED_INTENSITY) {
        baselineQCutVal =*std::max_element(nonPeakIntensities.begin(), nonPeakIntensities.end());
    } else if (params->eicBaselineEstimationType == EICBaselineEstimationType::EIC_NON_PEAK_MEDIAN_SMOOTHED_INTENSITY) {
        baselineQCutVal = mzUtils::median(nonPeakIntensities);
    }

    //recompute baseline if there are non-peak scans to use for baseline estimation.
    if (!nonPeakIntensities.empty()) {

        //drop all points above updated baselineQCutVal
        for(unsigned int i=0; i< intensity.size(); i++ ) {
            if (intensity[i] > baselineQCutVal) {
                baseline[i] = baselineQCutVal;
            } else {
                baseline[i] = intensity[i];
            }
        }

        //smooth baseline
        mzUtils::GaussianSmoother smoother(static_cast<unsigned int>(params->peakBaselineSmoothingWindow));
        vector<float> smoothed = smoother.smooth(baseline);
        for (unsigned int i = 0; i < smoothed.size(); i++) baseline[i] = smoothed.at(i);

        //count number of observation in EIC above baseline
        eic_noNoiseObs = 0;
        for(unsigned int i=0; i < intensity.size(); i++) {
            if(intensity[i] > baseline[i]) eic_noNoiseObs++;
        }

        if (debug) {
            cout << "\nBaseline with baselineQCutVal=" << baselineQCutVal << ":" << endl;
            for (unsigned int i=0; i < intensity.size(); i++) {
                cout << "i=" << i << ": " << baseline[i] << endl;
            }
        }
    } else if (debug) {
        cout << "EIC::computeBaselineByNonPeakIntensity(): nonPeakIntensities was empty, baseline will not be recomputed." << endl;
    }
}

void  EIC::subtractBaseLine() {

	if (baseline.size() == 0 ) {
        cout << "subtractBaseLine() failed. empty baseline vector\n";
		return;
	}

	eic_noNoiseObs=0;
	for(unsigned int i=0; i<intensity.size(); i++ ) {
		intensity[i] -= baseline[i];
		if(intensity[i]<0) intensity[i]=0;
        if(intensity[i]>0) eic_noNoiseObs++;
	}
}

void EIC::computeSpline(int smoothWindow) {
	int n = intensity.size();

    if (n == 0) return;
        if ( this->spline.size() ) { spline.clear(); }

        try {
                this->spline = vector<float>(n,0);
        }  catch(...) {
                cout << "Exception caught while allocating memory " << n << "floats " << endl;
        }

        //initalize spline, set to intensity vector
        for(int i=0; i<n; i++) spline[i] = intensity[i];

        if (smoothWindow > n/3 ) smoothWindow=n/3; //smoothing window is too large
        if (smoothWindow <= 1) return; 	//nothing to smooth get out

        if( smootherType == SAVGOL) { //SAVGOL SMOOTHER

            mzUtils::SavGolSmoother smoother(smoothWindow,smoothWindow,4);
            vector<float>smoothed = smoother.Smooth(intensity);
            for(int i=0; i<n; i++) spline[i] = smoothed[i];

        } else if (smootherType == GAUSSIAN) { //GAUSSIAN SMOOTHER

            mzUtils::GaussianSmoother smoother(smoothWindow, 3, 1);
            vector<float> smoothed = smoother.smooth(intensity);
            for (int i = 0; i<n; i++) spline[i] = smoothed[i];

        } else if ( smootherType == AVG) { //MOVING AVERAGE SMOOTHER

            mzUtils::MovingAverageSmoother smoother(smoothWindow);
            vector<float> smoothed = smoother.smooth(intensity);
            for (int i =0; i<n; i++) spline[i] = smoothed[i];

        }
}


Peak* EIC::addPeak(int peakPos) { 
    peaks.push_back(Peak(this,peakPos));
	return &peaks[peaks.size()-1];
}
		
void  EIC::getPeakPositions(int smoothWindow) {

    //Need to explicitly clear out peaks, else an extra peak will be added with each getPeaks() call.
    peaks.clear();

    //cout << "getPeakPositions() " << " sWindow=" << smoothWindow << " sType=" << smootherType << endl;

    unsigned int N = intensity.size();
    if ( N == 0 ) return;

    computeSpline(smoothWindow);
    if (spline.size() == 0) return;

    for (unsigned int i=1; i < N-1; i++ ) {
        if ( spline[i] > spline[i-1] && spline[i] > spline[i+1]) {
            addPeak(i);
        } else if ( spline[i] > spline[i-1] && spline[i] == spline[i+1] ) {
            float highpoint = spline[i];
            while(i<N-1) { 
                i++; 
                if ( spline[i+1] == highpoint) continue;
                if ( spline[i+1] > highpoint) break; 
                if ( spline[i+1] < highpoint) { addPeak(i); break; }
            }
        }
    }	

    //baseline always uses Gaussian smoothing.
    computeBaseLine(baselineSmoothingWindow, baselineDropTopX);

    getPeakStatistics();
}

/**
 * Description:
 * get the 5-point global maximum for this EIC.
 * This should only be used when the EIC is thought to contain exactly one valid peak,
 * which is also the most intense peak in the window.
 *
 * @brief EIC::getSingleGlobalMaxPeak
 * @param smoothWindow
 */
void EIC::getSingleGlobalMaxPeak(int smoothWindow) {

    //Need to explicitly clear out peaks, else an extra peak will be added with each getPeaks() call.
    peaks.clear();

    //cout << "getPeakPositions() " << " sWindow=" << smoothWindow << " sType=" << smootherType << endl;

    unsigned int N = intensity.size();
    if ( N == 0 ) return;

    computeSpline(smoothWindow);
    if (spline.size() == 0) return;

    //Isuse 698: Switch from 5-point max to 3-point max
    unsigned int globalMax = 0;
    for (unsigned int i=1; i < N-1; i++ ) {
        if ( spline[i] > spline[i-1] &&
             spline[i] > spline[i+1] &&
             spline[i] > spline[globalMax]) {
            globalMax = i;
        }
    }

    if (globalMax > 0) {
        addPeak(globalMax);
    }

    //baseline always uses Gaussian smoothing.
    computeBaseLine(baselineSmoothingWindow, baselineDropTopX);

    getPeakStatistics();
}

/**
 * @brief EIC::getPeakPositionsB
 * @param smoothWindow
 * @param intensityThreshold
 * Use intensity threshold, and write rt information to peak
 */
void EIC::getPeakPositionsB(int smoothWindow, float minSmoothedPeakIntensity) {

    //Need to explicitly clear out peaks, else an extra peak will be added with each getPeaks() call.
    peaks.clear();

    unsigned int N = intensity.size();
    if ( N == 0 ) return;

    computeSpline(smoothWindow);
    if (spline.size() == 0) return;

    for (unsigned int i=1; i < N-1; i++ ) {
        if (spline[i] < minSmoothedPeakIntensity ) continue;
        if (spline[i] > spline[i-1] && spline[i] > spline[i+1]) {
            addPeak(i);
        } else if ( spline[i] > spline[i-1] && spline[i] == spline[i+1] ) {
            float highpoint = spline[i];
            while(i<N-1) {
                i++;
                if ( spline[i+1] == highpoint) continue;
                if ( spline[i+1] > highpoint) break;
                if ( spline[i+1] < highpoint) {
                    addPeak(i);
                    break;
                }
            }
        }
    }

    //baseline always uses Gaussian smoothing.
    computeBaseLine(baselineSmoothingWindow, baselineDropTopX);

    getPeakStatistics();
}

/**
 * @brief EIC
 *
 * @param smoothWindow
 *
 * Approach that guarantees that peaks in the same EIC do not overlap. This approach is an alternative to
 * the EIC::findPeakBounds() method (in EIC::getPeakStatistics()).
 *
 * In this approach, peaks are identified using 3-point local maxima, just as in EIC::getPeakPositions()
 * and EIC::getPeakPositionsB().
 *
 * In EIC::getPeakPositionsB(), only peaks that satisfy a minimum intensity threshold are retained.
 * In EIC::getPeakPositions(), all peaks are retained.
 *
 * Here, only peaks that are above the baseline are retained.
 *
 * After that, minima are determined by finding the minimum smoothed intensity point between two maxima,
 *  or between a maximum and either end of the EIC.
 *
 * Peaks always range from one minima to another. If an EIC contains multiple maxima, adjacent maxima will share
 * a minimum between them.
 *
 * Issue 569: Alternative approach to determine RT bounds of a peak
 * @param rtBoundsMaxIntensityFraction
 * When computing Rt bounds, stop when the intensity crosses this threshold
 *
 * Issue 572: Slope-based peak boundary detection.
 * @param rtBoundsSlopeThreshold: once the slope between adjacent points in a peak is lower than this value,
 * a peak boundary has been reached.
 *
 * @author phillipseitzer
 * @since 20191104
 */
void EIC::getPeakPositionsC(int smoothWindow, bool debug, bool isComputePeakBounds, float rtBoundsMaxIntensityFraction, float rtBoundsSlopeThreshold) {

    peaks.clear();

    unsigned int N = static_cast<unsigned int>(intensity.size());
    if (N == 0) return;

    computeSpline(smoothWindow);
    if (spline.size() == 0) return;

    //baseline always uses Gaussian smoothing.
    computeBaseLine(baselineSmoothingWindow, baselineDropTopX);
    if (baseline.size() == 0) return;

    enum SplineAnnotation {
        MAX,
        MIN,
        NONE
    };

    vector<SplineAnnotation> splineAnnotation(N, SplineAnnotation::NONE);

    int lastMax = -1;
    int firstMax = -1;

    for (unsigned int i=1; i < N-1; i++ ) {

        //only consider peaks above baseline
        if (spline[i] > spline[i-1] && spline[i] > spline[i+1] && intensity[i] > baselineQCutVal) {

            addPeak(static_cast<int>(i));
            if (!isComputePeakBounds) peaks[peaks.size()-1].rt = rt[i];

            splineAnnotation[i] = SplineAnnotation::MAX;

            if (firstMax == -1){ //only assigned first time MAX appears
                firstMax = static_cast<int>(i);
            }
            lastMax = static_cast<int>(i); //reassigned every time MAX appears

        }
    }

    if (firstMax == -1) return; //no peaks determined based on 3-point max rule

    if (debug && !isComputePeakBounds) {
        cout << "WARNING: peak bounds were not determined.  Returning peaks with no peak boundaries defined." << endl;
    }

    if (!isComputePeakBounds) return;

    if (debug) cout << "baselineQCutVal=" << baselineQCutVal << endl;

    if (debug) {
        cout << "===================================" << endl;
        cout << "BEFORE ASSIGNING MINIMA:" << endl;
        for (unsigned int i = 0; i < N; i++) {
            string type = "";
            if (splineAnnotation[i] == SplineAnnotation::MAX) {
                type = "MAX";
            } else if (splineAnnotation[i] == SplineAnnotation::MIN) {
                type = "MIN";
            } else if (splineAnnotation[i] == SplineAnnotation::NONE) {
                type = "NONE";
            } else {
                cout << "ILLEGAL TYPE." << endl;
                abort();
            }
            cout << "i=" <<  i << " rt=" << rt[i] << " " << intensity[i] << " (spline=" << spline[i] << ") " << type << endl;
        }
        cout << "===================================" << endl;
    }

   //Issue 482: re-work peak boundary assignments (maxima and minima)
    peaks.clear();

    for (unsigned int i = 1; i < N-1; i++) {

        //indicates a new peak
        if (splineAnnotation[i] == SplineAnnotation::MAX) {

            if (debug) {
                cout << "Peak @ i=" << i << ":" << endl;
            }

            float smoothedPeakIntensity = spline[i];
            float intensityThreshold = baselineQCutVal;

            //Issue 569: Use rtBoundsMaxIntensityFraction instead of intensity value
            if (rtBoundsMaxIntensityFraction > 0.0f && rtBoundsMaxIntensityFraction < 1.0f) {
                intensityThreshold = max(intensityThreshold, rtBoundsMaxIntensityFraction * intensity[i]);
            }

            if (debug) {
                cout << "intensityThreshold: " << intensityThreshold << endl;
                cout << "rtBoundsMaxIntensityFraction: " << rtBoundsMaxIntensityFraction << endl;
                cout << "rtBoundsSlopeThreshold: " << rtBoundsSlopeThreshold << endl;
            }

            Peak* peak = addPeak(static_cast<int>(i));

            //descend to the left
            unsigned int leftMinimumIntensityIndex = i-1;
            unsigned int leftIndex = i-1;
            unsigned int leftNextIndex = i-2;

            while(true) {

                //if this point is below the baseline, it is invalid, stop immediately
                if (intensity[leftIndex] < intensityThreshold) {
                    break;
                }

                //if this point is another maximum, it is invalid, stop immediately
                if (splineAnnotation[leftIndex] == SplineAnnotation::MAX) {
                    break;
                }

                //this point is valid, compare for minimum intensity
                if (spline[leftIndex] < spline[leftMinimumIntensityIndex]) {
                    leftMinimumIntensityIndex = leftIndex;
                }

                //stop at the end of the EIC - this is the last valid point to the left.
                if (leftIndex == 0) {
                    break;
                }

                //Issue 572: Use slope based peak boundary detection
                if (rtBoundsSlopeThreshold > 0) {

                    float diff = spline[leftIndex] - spline[leftNextIndex];

                    //avoid possible divide-by-zero error
                    if (diff < 1e-10f) {
                        leftMinimumIntensityIndex = leftIndex;
                        if (debug) {
                            cout << "diff L: (" << leftIndex << ", " << leftNextIndex << "): ("
                                 << spline[leftIndex] << ", " << spline[leftNextIndex] << ")"
                                 << " diff = " << diff
                                 << endl;
                            cout << "slope < rtBoundsSlopeThreshold --> L=" <<  leftMinimumIntensityIndex << endl;
                        }
                        break;
                    }

                    float slope = ( diff / smoothedPeakIntensity) / (rt[leftIndex]-rt[leftNextIndex]);
                    if (debug) {
                        cout << "slope L: (" << leftIndex << ", " << leftNextIndex << "): ("
                             << spline[leftIndex] << ", " << spline[leftNextIndex] << ")"
                             << " diff = " << diff
                             << " frac = " << ( -1* diff / smoothedPeakIntensity)
                             << " slope = " << slope
                             << endl;
                    }

                    if (slope < rtBoundsSlopeThreshold) {
                        leftMinimumIntensityIndex = leftIndex;
                        if (debug) {
                            cout << "slope < rtBoundsSlopeThreshold --> L=" << leftMinimumIntensityIndex << endl;
                        }
                        break;
                    }
                }

                //continue moving to the left, check for more valid points
                leftIndex--;
                leftNextIndex--;
            }

            if (debug) {
                cout << "leftMinimumIntensityIndex: " << leftMinimumIntensityIndex << endl;
            }

            //descend to the right
            unsigned int rightMinimumIntensityIndex = i+1;
            unsigned int rightIndex = i+1;
            unsigned int rightNextIndex = i+2;

            while(true) {

                //if this point is below the baseline, it is invalid, stop immediately
                if (intensity[rightIndex] < intensityThreshold) {
                    break;
                }

                //if this point is another maximum, it is invalid, stop immediately
                if (splineAnnotation[rightIndex] == SplineAnnotation::MAX) {
                    break;
                }

                //this point is valid, compare for minimum intensity
                if (spline[rightIndex] < spline[rightMinimumIntensityIndex]) {
                    rightMinimumIntensityIndex = rightIndex;
                }

                //stop at the end of the EIC - this is the last valid point to the right
                if (rightIndex == N-1) {
                    break;
                }

                //Issue 572: Use slope based peak boundary detection
                if (rtBoundsSlopeThreshold > 0 && rightNextIndex <= N-1) {

                    float diff = spline[rightIndex] - spline[rightNextIndex];

                    //avoid possible divide-by-zero error
                    if (diff < 1e-10f) {
                        rightMinimumIntensityIndex = rightIndex;
                        if (debug) {
                            cout << "diff R: (" << rightIndex << ", " << rightNextIndex << "): ("
                                 << spline[rightIndex] << ", " << spline[rightNextIndex] << ")"
                                 << " diff = " << diff
                                 << endl;
                            cout << "slope < rtBoundsSlopeThreshold --> R=" << rightMinimumIntensityIndex << endl;
                        }
                        break;
                    }

                    float slope = ( -1* diff / smoothedPeakIntensity)/(rt[rightIndex]-rt[rightNextIndex]);
                    if (debug) {
                        cout << "slope R: (" << rightIndex << ", " << rightNextIndex << "): ("
                             << spline[rightIndex] << ", " << spline[rightNextIndex] << ")"
                             << " diff = " << diff
                             << " frac = " << ( -1* diff / smoothedPeakIntensity)
                             << " slope = " << slope
                             << endl;
                    }

                    if (slope < rtBoundsSlopeThreshold) {
                        rightMinimumIntensityIndex = rightIndex;
                        if (debug) {
                            cout << "slope < rtBoundsSlopeThreshold --> R=" << rightMinimumIntensityIndex << endl;
                        }
                        break;
                    }
                }

                //continue moving to the right, check for more valid points
                rightIndex++;
                rightNextIndex++;
            }

            if (debug) {
                cout << "rightMinimumIntensityIndex: " << rightMinimumIntensityIndex << endl;
            }

            peak->minpos = leftMinimumIntensityIndex;
            peak->rtmin = rt[leftMinimumIntensityIndex];
            peak->mzmin = mz[leftMinimumIntensityIndex];
            peak->minscan = static_cast<unsigned int>(scannum[leftMinimumIntensityIndex]);

            peak->maxpos = rightMinimumIntensityIndex;
            peak->rtmax = rt[rightMinimumIntensityIndex];
            peak->mzmax = mz[rightMinimumIntensityIndex];
            peak->maxscan = static_cast<unsigned int>(scannum[rightMinimumIntensityIndex]);
        }
    }

    // TODO: recompute baseline based on peak boundaries, if necessary

    //Issue 549: get peak details after the fact
    for (auto& peak : peaks) {
        getPeakDetails(peak, false);
    }

    //assign peak ranks based on total area of the peak
    sort(peaks.begin(), peaks.end(), Peak::compArea);
    for(unsigned int i=0; i<peaks.size(); i++) peaks[i].peakRank = i;

    //sort peaks by RT
    sort(peaks.begin(), peaks.end(), Peak::compRt);

    if (debug) {
        for (auto peak : peaks) {
            cout << "PEAK: pos=" << peak.pos << ", mz=" << peak.peakMz << ", rt=" << peak.rt << endl;
            cout << "\t min=" << peak.minpos << ", minmz=" << peak.mzmin << ", rtmin=" << peak.rtmin << endl;
            cout << "\t max=" << peak.maxpos << ", maxmz=" << peak.mzmax << ", rtmax=" << peak.rtmax << endl;
            cout << "\t SN=" << peak.signalBaselineRatio << ", (" << peak.peakIntensity << "/" << peak.peakBaseLineLevel << ")" << endl;
        }
    }

}

void EIC::getPeakPositionsD(shared_ptr<PeakPickingAndGroupingParameters> params, bool debug){

    peaks.clear();

    unsigned int N = static_cast<unsigned int>(intensity.size());
    if (N == 0) return;

    computeSpline(params->peakSmoothingWindow);
    if (spline.size() == 0) return;

    //baseline always uses Gaussian smoothing.
    computeBaseLine(params->peakBaselineSmoothingWindow, params->peakBaselineDropTopX);
    if (baseline.size() == 0) return;

    enum SplineAnnotation {
        MAX,
        MIN,
        NONE
    };

    vector<SplineAnnotation> splineAnnotation(N, SplineAnnotation::NONE);

    int lastMax = -1;
    int firstMax = -1;

    for (unsigned int i=1; i < N-1; i++ ) {

        //only consider peaks above baseline
        if (spline[i] > spline[i-1] && spline[i] > spline[i+1] && intensity[i] > baselineQCutVal) {

            addPeak(static_cast<int>(i));
            if (!params->peakIsComputeBounds) peaks[peaks.size()-1].rt = rt[i];

            splineAnnotation[i] = SplineAnnotation::MAX;

            if (firstMax == -1){ //only assigned first time MAX appears
                firstMax = static_cast<int>(i);
            }
            lastMax = static_cast<int>(i); //reassigned every time MAX appears

        }
    }

    //Issue 740: Support option to pick maxima on the edge of a trace as peaks, if
    //they are higher than the point directly adjacent to them.
    //This option requires at least 3 points; a trace of 2 points would always have one of the edges picked
    //a peak, which doesn't make a lot of sense.
    if (params->peakIsPickEdgePeaks && N >= 3) {
        if (spline[0] > spline[1]) {
            splineAnnotation[0] = SplineAnnotation::MAX;
            if (firstMax == -1) {
                firstMax = 0;
            }
        }
        if (spline[N-1] > spline[N-2]){
            splineAnnotation[N-1] = SplineAnnotation::MAX;
            if (firstMax == -1) {
                firstMax = N-1;
            }
        }
    }

    if (firstMax == -1) return; //no peaks determined based on 3-point max rule

    if (debug && !params->peakIsComputeBounds) {
        cout << "WARNING: peak bounds were not determined.  Returning peaks with no peak boundaries defined." << endl;
    }

    if (!params->peakIsComputeBounds) return;

    if (debug) cout << "baselineQCutVal=" << baselineQCutVal << endl;

    if (debug) {
        cout << "===================================" << endl;
        cout << "BEFORE ASSIGNING MINIMA:" << endl;
        for (unsigned int i = 0; i < N; i++) {
            string type = "";
            if (splineAnnotation[i] == SplineAnnotation::MAX) {
                type = "MAX";
            } else if (splineAnnotation[i] == SplineAnnotation::MIN) {
                type = "MIN";
            } else if (splineAnnotation[i] == SplineAnnotation::NONE) {
                type = "NONE";
            } else {
                cout << "ILLEGAL TYPE." << endl;
                abort();
            }
            cout << "i=" <<  i << " rt=" << rt[i] << " " << intensity[i] << " (spline=" << spline[i] << ") " << type << endl;
        }
        cout << "===================================" << endl;
    }

   //Issue 482: re-work peak boundary assignments (maxima and minima)
    peaks.clear();

    for (unsigned int i = 0; i < N; i++) {

        //indicates a new peak
        if (splineAnnotation[i] == SplineAnnotation::MAX) {

            if (debug) {
                cout << "Peak @ i=" << i << ":" << endl;
            }

            float smoothedPeakIntensity = spline[i];
            float halfMaxIntensity = 0.5f * smoothedPeakIntensity;
            float intensityThreshold = baselineQCutVal;

            //Issue 569: Use rtBoundsMaxIntensityFraction instead of intensity value
            if (params->peakRtBoundsMaxIntensityFraction > 0.0f && params->peakRtBoundsMaxIntensityFraction < 1.0f) {
                intensityThreshold = max(intensityThreshold, params->peakRtBoundsMaxIntensityFraction * intensity[i]);
            }

            if (debug) {
                cout << "intensityThreshold: " << intensityThreshold << endl;
                cout << "rtBoundsMaxIntensityFraction: " << params->peakRtBoundsMaxIntensityFraction << endl;
                cout << "rtBoundsSlopeThreshold: " << params->peakRtBoundsSlopeThreshold << endl;
                cout << "halfMaxIntensity: " << halfMaxIntensity << endl;
            }

            Peak* peak = addPeak(static_cast<int>(i));

            //descend to the left
            unsigned int leftMinimumIntensityIndex = i-1;
            unsigned int leftIndex = i-1;
            unsigned int leftNextIndex = i-2;

            bool reachedLeftHalfMax = false;
            unsigned int leftHalfMaxIntensityIndex = leftIndex; //initialization

            if (debug) {
                cout << "\nSTART LEFT DESCENT" << endl;
            }

            while(true) {

                //left-edge peak - special case from params->peakIsPickEdgePeaks
                if (i == 0) {
                    leftMinimumIntensityIndex = 0;
                    break;
                }

                if (debug) {
                    cout << "Evaluating " << leftIndex << ": "
                         << "rt=" << rt[leftIndex]
                         << ", intensity=" << intensity[leftIndex]
                         << ", intensityThreshold=" << intensityThreshold
                         << ", reachedLeftHalfMax? " << (reachedLeftHalfMax ? "yes" : "no")
                         << endl;
                }

                //if this point is below the baseline, it is invalid, stop immediately
                if (intensity[leftIndex] < intensityThreshold) {
                    break;
                }

                //if this point is another maximum, it is invalid, stop immediately
                if (splineAnnotation[leftIndex] == SplineAnnotation::MAX) {
                    break;
                }

                //Issue 603: half-max quant metrics
                if (!reachedLeftHalfMax && spline[leftIndex] < halfMaxIntensity) {
                    if (debug) {
                        cout << "FWHM Left: " << leftIndex
                             << ", spline=" << spline[leftIndex]
                             << endl;
                    }
                    reachedLeftHalfMax = true;
                }

                if (!reachedLeftHalfMax) {
                    leftHalfMaxIntensityIndex = leftIndex;
                }

                //this point is valid, compare for minimum intensity
                if (spline[leftIndex] < spline[leftMinimumIntensityIndex]) {
                    leftMinimumIntensityIndex = leftIndex;
                }

                //stop at the end of the EIC - this is the last valid point to the left.
                if (leftIndex == 0) {
                    break;
                }

                //Issue 572: Use slope based peak boundary detection
                if (params->peakRtBoundsSlopeThreshold > 0) {

                    float diff = spline[leftIndex] - spline[leftNextIndex];

                    //avoid possible divide-by-zero error
                    if (diff < 1e-10f) {
                        leftMinimumIntensityIndex = leftIndex;
                        if (debug) {
                            cout << "diff L: (" << leftIndex << ", " << leftNextIndex << "): ("
                                 << spline[leftIndex] << ", " << spline[leftNextIndex] << ")"
                                 << " diff = " << diff
                                 << endl;
                            cout << "slope < rtBoundsSlopeThreshold --> L=" <<  leftMinimumIntensityIndex << endl;
                        }
                        break;
                    }

                    float slope = ( diff / smoothedPeakIntensity) / (rt[leftIndex]-rt[leftNextIndex]);
                    if (debug) {
                        cout << "slope L: (" << leftIndex << ", " << leftNextIndex << "): ("
                             << spline[leftIndex] << ", " << spline[leftNextIndex] << ")"
                             << " diff = " << diff
                             << " frac = " << ( -1* diff / smoothedPeakIntensity)
                             << " slope = " << slope
                             << endl;
                    }

                    if (slope < params->peakRtBoundsSlopeThreshold) {
                        leftMinimumIntensityIndex = leftIndex;
                        if (debug) {
                            cout << "slope < rtBoundsSlopeThreshold --> L=" << leftMinimumIntensityIndex << endl;
                        }
                        break;
                    }
                }

                //continue moving to the left, check for more valid points
                leftIndex--;
                leftNextIndex--;
            }

            if (debug) {
                cout << "leftMinimumIntensityIndex: " << leftMinimumIntensityIndex << "\n" << endl;
            }

            //descend to the right
            unsigned int rightMinimumIntensityIndex = i+1;
            unsigned int rightIndex = i+1;
            unsigned int rightNextIndex = i+2;

            bool reachedRightHalfMax = false;
            unsigned int rightHalfMaxIntensityIndex = rightIndex; //initialization

            if (debug) {
                cout << "\nSTART RIGHT DESCENT" << endl;
            }

            while(true) {

                //right-edge peak - special case from params->peakIsPickEdgePeaks
                if (i == (N-1)) {
                    rightMinimumIntensityIndex = N-1;
                    break;
                }

                if (debug) {
                    cout << "Evaluating " << rightIndex << ": "
                         << "rt=" << rt[rightIndex]
                         << ", intensity=" << intensity[rightIndex]
                         << ", intensityThreshold=" << intensityThreshold
                         << ", reachedRightHalfMax? " << (reachedRightHalfMax ? "yes" : "no")
                         << endl;
                }

                //if this point is below the baseline, it is invalid, stop immediately
                if (intensity[rightIndex] < intensityThreshold) {
                    break;
                }

                //if this point is another maximum, it is invalid, stop immediately
                if (splineAnnotation[rightIndex] == SplineAnnotation::MAX) {
                    break;
                }

                //Issue 603: half-max quant metrics
                if (!reachedRightHalfMax && spline[rightHalfMaxIntensityIndex] < halfMaxIntensity) {
                    if (debug) {
                        cout << "FWHM Right: " << rightIndex
                             << ", spline=" << spline[rightIndex]
                             << endl;
                    }
                    reachedRightHalfMax = true;
                }

                if (!reachedRightHalfMax) {
                    rightHalfMaxIntensityIndex = rightIndex;
                }

                //this point is valid, compare for minimum intensity
                if (spline[rightIndex] < spline[rightMinimumIntensityIndex]) {
                    rightMinimumIntensityIndex = rightIndex;
                }

                //stop at the end of the EIC - this is the last valid point to the right
                if (rightIndex == N-1) {
                    break;
                }

                //Issue 572: Use slope based peak boundary detection
                if (params->peakRtBoundsSlopeThreshold > 0 && rightNextIndex <= N-1) {

                    float diff = spline[rightIndex] - spline[rightNextIndex];

                    //avoid possible divide-by-zero error
                    if (diff < 1e-10f) {
                        rightMinimumIntensityIndex = rightIndex;
                        if (debug) {
                            cout << "diff R: (" << rightIndex << ", " << rightNextIndex << "): ("
                                 << spline[rightIndex] << ", " << spline[rightNextIndex] << ")"
                                 << " diff = " << diff
                                 << endl;
                            cout << "slope < rtBoundsSlopeThreshold --> R=" << rightMinimumIntensityIndex << endl;
                        }
                        break;
                    }

                    float slope = ( -1* diff / smoothedPeakIntensity)/(rt[rightIndex]-rt[rightNextIndex]);
                    if (debug) {
                        cout << "slope R: (" << rightIndex << ", " << rightNextIndex << "): ("
                             << spline[rightIndex] << ", " << spline[rightNextIndex] << ")"
                             << " diff = " << diff
                             << " frac = " << ( -1* diff / smoothedPeakIntensity)
                             << " slope = " << slope
                             << endl;
                    }

                    if (slope < params->peakRtBoundsSlopeThreshold) {
                        rightMinimumIntensityIndex = rightIndex;
                        if (debug) {
                            cout << "slope < rtBoundsSlopeThreshold --> R=" << rightMinimumIntensityIndex << endl;
                        }
                        break;
                    }
                }

                //continue moving to the right, check for more valid points
                rightIndex++;
                rightNextIndex++;
            }

            if (debug) {
                cout << "rightMinimumIntensityIndex: " << rightMinimumIntensityIndex  << "\n" << endl;
            }

            peak->minpos = leftMinimumIntensityIndex;
            peak->rtmin = rt[leftMinimumIntensityIndex];
            peak->mzmin = mz[leftMinimumIntensityIndex];
            peak->minscan = static_cast<unsigned int>(scannum[leftMinimumIntensityIndex]);

            peak->maxpos = rightMinimumIntensityIndex;
            peak->rtmax = rt[rightMinimumIntensityIndex];
            peak->mzmax = mz[rightMinimumIntensityIndex];
            peak->maxscan = static_cast<unsigned int>(scannum[rightMinimumIntensityIndex]);

            if (reachedLeftHalfMax && reachedRightHalfMax) {

                peak->minPosFWHM = leftHalfMaxIntensityIndex;
                peak->rtminFWHM = rt[leftHalfMaxIntensityIndex];
                peak->minScanFWHM = static_cast<unsigned int>(scannum[leftHalfMaxIntensityIndex]);

                peak->maxPosFWHM = rightHalfMaxIntensityIndex;
                peak->rtmaxFWHM = rt[rightHalfMaxIntensityIndex];
                peak->maxScanFWHM = static_cast<unsigned int>(scannum[rightHalfMaxIntensityIndex]);
            }
        }
    }

    //Issue 549: recompute baseline before computing peak details, if appropriate.
    if (params->eicBaselineEstimationType != EICBaselineEstimationType::DROP_TOP_X) {
        computeBaselineByNonPeakIntensity(params, debug);
    }

    for (auto& peak : peaks) {
        getPeakDetails(peak, params->peakIsReassignPosToUnsmoothedMax);
    }

    //assign peak ranks based on total area of the peak
    sort(peaks.begin(), peaks.end(), Peak::compArea);
    for(unsigned int i=0; i<peaks.size(); i++) peaks[i].peakRank = i;

    //sort peaks by RT
    sort(peaks.begin(), peaks.end(), Peak::compRt);

    if (debug) {
        cout << "===================================\n";
        for (auto peak : peaks) {
            cout << "PEAK: pos=" << peak.pos << ", mz=" << peak.peakMz << ", rt=" << peak.rt << "\n";
            cout << "\t min=" << peak.minpos << ", minmz=" << peak.mzmin << ", rtmin=" << peak.rtmin << "\n";
            cout << "\t max=" << peak.maxpos << ", maxmz=" << peak.mzmax << ", rtmax=" << peak.rtmax << "\n";

            cout << "FWHM stats:\n";
            cout << "\tFWHM min=" << peak.minPosFWHM << ", FWHM rtmin=" << peak.rtminFWHM << "\n";
            cout << "\tFWHM max=" << peak.maxPosFWHM << ", FWHM rtmax=" << peak.rtmaxFWHM << "\n";

            cout << "Quant Stats:\n";
            cout << "\t SN=" << peak.signalBaselineRatio << ", (" << peak.peakIntensity << "/" << peak.peakBaseLineLevel << ")" << "\n";
            cout << "\t Smoothed SN=" << peak.smoothedSignalBaselineRatio << ", (" << peak.smoothedIntensity << "/" << peak.peakBaseLineLevel << ")" << "\n";
            cout << "\t peakArea=" << peak.peakArea << ", smoothedPeakArea=" << peak.smoothedPeakArea << "\n";
            cout << "\t FWHM peakArea=" << peak.peakAreaFWHM << ", FWHM smoothedPeakArea=" << peak.smoothedPeakAreaFWHM << "\n" << endl;
        }
    }
}

void EIC::findPeakBounds(Peak& peak) {
	int apex = peak.pos;

	int ii = apex-1;
	int jj = apex+1;
	int lb = ii;
	int rb = jj; 

	unsigned int N = intensity.size();
	if (N==0) return;
    if (spline.empty())   return;
    if (baseline.empty()) return;

    //cout << "findPeakBounds:" << apex << " " << rt[apex] << endl;

	int directionality = 0;
	float lastValue = spline[apex];
	while(ii > 0 && ii < (int) N) { //walk left
		float relSlope = (spline[ii]-lastValue)/lastValue;
		relSlope > 0.01 ? directionality++ : directionality=0;
		//if (spline[ii]<=spline[lb] ) lb=ii; 
		if (intensity[ii]<=intensity[lb] ) lb=ii; 
		if (spline[ii] == 0 ) break;
		if (spline[ii]<=baseline[ii]) break;
		if (spline[ii]<=spline[apex]*0.01) break;

		if (directionality >= 2) break;
		lastValue=spline[ii];
		ii=ii-1;
	}

	directionality = 0;
	lastValue = spline[apex];

	while(jj >0 && jj < (int) N ){ //walk right
		float relSlope = (spline[jj]-lastValue)/lastValue;
		relSlope > 0.01 ? directionality++ : directionality=0;
		//if (spline[jj]<=spline[rb] ) rb=jj;
		if (intensity[jj]<=intensity[rb] ) rb=jj;
		if (spline[jj] == 0 ) break;
		if (spline[jj]<=baseline[ii]) break;
		if (spline[jj]<=spline[apex]*0.01) break;

		if (directionality >= 2) break;
		lastValue=spline[jj];
		jj=jj+1;
	}

	//find maximum point in the span from min to max position
	for(unsigned int k=lb; k<rb && k < N; k++ ) {
		if(intensity[k]> intensity[peak.pos] && mz[k] > 0 ) peak.pos = k;
	}

	//remove zero intensity points on the left
	for(unsigned int k=lb; k<peak.pos && k < N; k++ ) {
		if (intensity[k] > 0 ) break;
		lb = k;
	}

	//remove zero intensity points on the right
	for(unsigned int k=rb; k>peak.pos && k < N; k-- ) {
		if (intensity[k] > 0 ) break; 
		rb = k;
	}

	//for rare cases where peak is a single observation
	if (lb == apex && lb-1 >= 0) lb=apex-1;
	if (rb == apex && rb+1 < N) rb=apex+1;

	peak.minpos = lb;
	peak.maxpos = rb;
    //cout << "\tfindPeakBounds:" << lb << " " << rb << " " << rb-lb+1 << endl;
}

/**
 * @brief EIC::getPeakDetails
 * @param peak
 * @param isCorrectPeakByMaxIntensity
 *
 * Issue 381: "isCorrectPeakByMaxIntensity" should be false when using spline-based maxima finding,
 * as in EIC::getPeakPositionsC().
 *
 * Sometimes the intensity maxima will disagree with the spline maxima, which can lead to peak boundaries
 * being reassigned, and can cause pathological behavior.
 *
 */
void EIC::getPeakDetails(Peak& peak, bool isCorrectPeakByMaxIntensity) {

    unsigned long N = intensity.size();

    if (N == 0) return;
    if (peak.pos >= N) return;

    //intensity and mz at the apex of the peaks
    peak.peakIntensity = intensity[peak.pos];
    peak.smoothedIntensity = spline[peak.pos];

    peak.noNoiseObs = 0;
    peak.peakAreaCorrected  = 0;
    peak.peakArea=0;
    float baselineArea=0;

    if (sample && sample->isBlank) {
        peak.fromBlankSample = true;
    }

    StatisticsVector<float>allmzs;
    string bitstring;
    if(peak.maxpos >= N) peak.maxpos= static_cast<unsigned int>(N)-1;
    if(peak.minpos >= N) peak.minpos = peak.pos; //unsigned number weirdness.

    float lastValue = intensity[peak.minpos];

    for(unsigned int j=peak.minpos; j<= peak.maxpos;j++ ){

        peak.peakArea += intensity[j];
        peak.smoothedPeakArea += spline[j];
        baselineArea +=   baseline[j];

        if (intensity[j] > baseline[j]) peak.noNoiseObs++;

        if (isCorrectPeakByMaxIntensity) {
            if(peak.peakIntensity < intensity[j]) {
                peak.peakIntensity = intensity[j];
                peak.smoothedIntensity = spline[j];
                peak.pos = j;
            }
        }

        if (mz.size() > 0 && mz[j] > 0) allmzs.push_back(mz[j]);

        if (intensity[j] <= baseline[j]) {
                bitstring += "0";
        } else if (intensity[j] > lastValue) {
                  bitstring += "+";
        } else if( intensity[j] < lastValue) {
                  bitstring += "-";
        } else if( intensity[j] == lastValue) {
                if (bitstring.length()>1) bitstring += bitstring[bitstring.length()-1]; else bitstring += "0";
        }

        lastValue = intensity[j];
    }

    if (peak.minPosFWHM > 0 && peak.maxPosFWHM > 0) {
        for (unsigned int j = peak.minPosFWHM; j<= peak.maxPosFWHM; j++ ){
            peak.peakAreaFWHM += intensity[j];
            peak.smoothedPeakAreaFWHM += spline[j];
        }
    }

    getPeakWidth(peak);

    if (rt.size() > 0 && rt.size() == N ) {
        peak.rt =    rt[ peak.pos ];
        peak.rtmin = rt[ peak.minpos ];
        peak.rtmax = rt[ peak.maxpos ];
    }

    if (scannum.size() && scannum.size() == N ) {
        peak.scan    = scannum[ peak.pos ];		//scan number at the apex of the peak
        peak.minscan = scannum[ peak.minpos ];	//scan number at left most bound
        peak.maxscan = scannum[ peak.maxpos ];	//scan number at the right most bound
    }

    int n =1;
    peak.peakAreaTop = intensity[peak.pos];
    peak.smoothedPeakAreaTop = spline[peak.pos];

    if (peak.pos-1 < N){
        peak.peakAreaTop += intensity[peak.pos-1];
        peak.smoothedPeakAreaTop += spline[peak.pos-1];
        n++;
    }

    if (peak.pos+1 < N){
        peak.peakAreaTop += intensity[peak.pos+1];
        peak.smoothedPeakAreaTop += spline[peak.pos+1];
        n++;
    }

    //Issue 603, 632: peakBaselineLevel cannot be 0, causes problems with computing S:N ratios
    float peakEdgeMaxIntensity = MAX(intensity[peak.minpos], intensity[peak.maxpos]);
    float peakEdgeMinIntensity = MIN(intensity[peak.minpos], intensity[peak.maxpos]);
    float peakBaselineLevel = 0;

    if (baseline[peak.pos] > 0) {
        peakBaselineLevel = baseline[peak.pos];
    } else if (peakEdgeMinIntensity > 0) {
        peakBaselineLevel = peakEdgeMinIntensity;
    } else if (peakEdgeMaxIntensity > 0) {
        peakBaselineLevel = peakEdgeMaxIntensity;
    } else {
        peakBaselineLevel = 10; // last resort
    }


    peak.peakMz = mz[ peak.pos ];
    peak.peakAreaTop /= n;
    peak.peakBaseLineLevel = peakBaselineLevel;
    peak.noNoiseFraction = (float) peak.noNoiseObs/(this->eic_noNoiseObs+1);
    peak.peakAreaCorrected = peak.peakArea-baselineArea;
    peak.peakAreaFractional = peak.peakAreaCorrected/(totalIntensity+1);
    peak.signalBaselineRatio = peak.peakIntensity/peak.peakBaseLineLevel; //Issue 603: was maxBaseLine

    peak.smoothedPeakAreaCorrected = peak.smoothedPeakArea-baselineArea;
    peak.smoothedSignalBaselineRatio = peak.smoothedIntensity/peak.peakBaseLineLevel; //Issue 603: was maxBaseline

    if (allmzs.size()> 0 ) {
        peak.medianMz = allmzs.median();
        peak.baseMz =   allmzs.mean();
        peak.mzmin  =   allmzs.minimum();
        peak.mzmax  =   allmzs.maximum();
    }

    if ( peak.medianMz == 0) { peak.medianMz = peak.peakMz; }
    //cout << peak.peakMz << " " << peak.medianMz << " " << bitstring << endl;

	mzPattern p(bitstring);
    if (peak.width >= 5) peak.symmetry =  p.longestSymmetry('+','-');
    checkGaussianFit(peak);
}



void EIC::getPeakWidth(Peak& peak) { 

	int width=1;
	int left=0;
	int right=0;
	unsigned int N = intensity.size();

	for(unsigned int i=peak.pos-1; i>peak.minpos && i < N; i--) {
           if(intensity[i] > baseline[i]) left++;  else break;
        }

        for(unsigned int j=peak.pos+1; j<peak.maxpos && j < N; j++ ) {
           if(intensity[j] > baseline[j]) right++; else break;
        }

	peak.width = width + left + right;
}

vector<mzPoint> EIC::getIntensityVector(Peak& peak) { 
	vector<mzPoint>y;

	if (intensity.size()>0 ) {
        unsigned int  maxi=peak.maxpos;
        unsigned int  mini=peak.minpos;
        if(maxi >= intensity.size()) maxi=intensity.size()-1;

        for(unsigned int i=mini; i <= maxi; i++ ) {
            if(baseline.size() and intensity[i] > baseline[i] )  {
                y.push_back(mzPoint(rt[i],intensity[i],mz[i]));
            } else {
                y.push_back(mzPoint(rt[i],intensity[i],mz[i]));
            }
        }
	}
	return y;
}


void EIC::checkGaussianFit(Peak& peak) { 

        //set default values for fit (in case fit cannot be assessed).
        peak.gaussFitSigma = 0;
        peak.gaussFitR2 = 0.03;

		int left =  peak.pos - peak.minpos;
		int right = peak.maxpos - peak.pos;
		if (left <= 0 || right <= 0 ) return;
		int moves = min(left, right);
		if (moves < 3 ) return;

		//copy intensities into seperate vector
		//dim
		vector<float>pints(moves*2+1);

		int j=peak.pos+moves; if (j >= intensity.size() ) j=intensity.size()-1; if(j<1) j=1;
		int i=peak.pos-moves; if (i < 1 ) i=1;

		int k=0;
		for(; i<=j; i++) { pints[k]=intensity[i]; k++; }

        //   <sigma, minR>
        pair<float, float> gaussFitParams = mzUtils::gaussFit(pints, peak.gaussFitSigma, peak.gaussFitR2);

        peak.gaussFitSigma = gaussFitParams.first;
        peak.gaussFitR2 = gaussFitParams.second;

        //cout << "\tcheckGaussianFit(): Best Sigma=" << peak.gaussFitSigma <<  " minRsqr=" << peak.gaussFitR2 << endl;
}



void  EIC::getPeakStatistics() {

    for(unsigned int i=0; i<peaks.size(); i++) {
        findPeakBounds(peaks[i]);
        getPeakDetails(peaks[i]);
    }

	removeOverlapingPeaks();

	//assign peak ranks based on total area of the peak
    sort(peaks.begin(),peaks.end(),Peak::compArea);
    for(unsigned int i=0; i<peaks.size(); i++) peaks[i].peakRank = i; 
}

void  EIC::deletePeak(unsigned int i) {
    if (i<peaks.size()) {
     	peaks.erase( peaks.begin()+i );
    }
}

void EIC::summary() { 
        cout << "EIC: mz=" << mzmin <<  "-" << mzmax << " rt=" << rtmin << "-" << rtmax << endl;
        cout << "   : maxIntensity=" << maxIntensity << endl;
        cout << "   : peaks=" << peaks.size() << endl;
}		

void EIC::removeLowRankGroups( vector<PeakGroup>& groups, unsigned int rankLimit ) {
	if (groups.size() < rankLimit ) return;
	std::sort(groups.begin(), groups.end(), PeakGroup::compIntensity);
	for(unsigned int i=0; i < groups.size(); i++ ) {
		if( i > rankLimit ) { 
			groups.erase(groups.begin()+i); i--; 
        }
	}
}

/**
 * @brief EIC::groupPeaksB
 * @param eics
 * @param smoothingWindow
 * @param maxRtDiff
 * @param noiseThreshold
 * this noise threshold applies to the actual raw signal in the files
 *
 * @return
 * peaks are picked in individual EICs, and grouped together based on similarity in RT.
 *
 * peaks are not retained if they are below the noiseThreshold.
 * Note that the smoothed value is compared to the noiseThreshold, not the raw intensity.
 *
 * @deprecated in favor of EIC::groupPeaksC()
 */
vector<PeakGroup> EIC::groupPeaksB(vector<EIC*>& eics, int smoothingWindow, float maxRtDiff, float minSmoothedPeakIntensity) {

        //TODO: check that this flag is not too slow
        bool debug = false;

        if (debug){
            cout <<"smoothingWindow=" << smoothingWindow << ", maxRtDiff=" << maxRtDiff << endl;
        }

        //list filled and return by this function
        vector<PeakGroup> pgroups;

        //case with empty eics
        if (eics.empty()) return pgroups;

        //case there is only a single EIC, there is nothing to group
        if (eics.size() == 1 && eics[0]) {
            EIC* m=eics[0];
            for(unsigned int i=0; i< m->peaks.size(); i++ ) {
                PeakGroup grp;
                grp.groupId = static_cast<int>(i);
                grp.addPeak(m->peaks[i]);
                grp.groupStatistics();
                pgroups.push_back(grp);
            }
            if (debug) {
                cout << "Case with 1 eic produces " << pgroups.size() <<  "peak groups." << endl;
            }
            return pgroups;
        }

        int numTotalPeaks = 0;
        for (auto eic : eics){
            eic->getPeakPositionsB(smoothingWindow, minSmoothedPeakIntensity);
            numTotalPeaks += eic->peaks.size();
        }

        if (debug) {
            cout << "Discovered " << numTotalPeaks << " peaks in " << eics.size() << " samples." << endl;
        }

                  //<sample id,   peak>
        vector<pair<unsigned int, Peak>> peakSamplePairs = vector<pair<unsigned int,Peak>>(numTotalPeaks);

        int k = 0;
        for (unsigned int i = 0; i < eics.size(); i++) {
            EIC *eic = eics.at(i);
            for (auto peak : eic->peaks) {
                peakSamplePairs.at(k) = make_pair(i, peak);
                k++;
            }
        }

        sort(peakSamplePairs.begin(), peakSamplePairs.end(),
             [](const pair<int, Peak>& lhs, const pair<int, Peak>& rhs){
                return lhs.second.rt - rhs.second.rt < 0;
            });

        vector<pair<double, pair<unsigned int, unsigned int>>> dissimilarities;

        for (unsigned int i = 0; i < peakSamplePairs.size(); i++){

            pair<unsigned int, Peak> peakPairI = peakSamplePairs.at(i);

            for (unsigned int j = i+1; j < peakSamplePairs.size(); j++) {

                pair<unsigned int, Peak> peakPairJ = peakSamplePairs.at(j);

                //skip peaks from the same sample.
                if (peakPairI.first == peakPairJ.first) continue;

                float deltaRt = peakPairJ.second.rt - peakPairI.second.rt;

                //out of tolerance condition
                if (deltaRt > maxRtDiff) {
                    continue;
                }

                //else, create a pair
                dissimilarities.push_back(make_pair(deltaRt, make_pair(i, j)));
            }
        }

        if (debug) {
            cout << "Computed " << dissimilarities.size() << " dissimilarities." << endl;
        }

        sort(dissimilarities.begin(), dissimilarities.end(),
             [](const pair<double, pair<unsigned int, unsigned int>>& lhs,
                const pair<double, pair<unsigned int, unsigned int>>& rhs){
            if (abs(lhs.first - rhs.first) < 1e-6) {
              if (lhs.second.first == rhs.second.first) {
                return lhs.second.second < rhs.second.second;
              } else {
                return lhs.second.first < rhs.second.first;
              }
            } else {
              return lhs.first < rhs.first;
            }
        });

        if (debug) {
            cout << "Dissimilarities: " << endl;
            for (auto diss : dissimilarities) {
                cout << "(" << diss.second.first << ", " << diss.second.second << "): " << diss.first << endl;
            }
        }

        // <unsigned int> --> peakSamplePairs index
        //Initially, all peaks start in their own clusters
        vector<vector<unsigned int>> peakGroups = vector<vector<unsigned int>> (peakSamplePairs.size());
        for (unsigned int i = 0; i < peakSamplePairs.size(); i++){
            peakGroups.at(i) = {i};
        }

        if (debug) {
           cout << endl;
           cout << "peakGroups initial status:" << endl;
           for (auto peakGroup : peakGroups) {
               if (peakGroup.empty()) continue;
               cout << "peakGroup Indexes: ";
               for (auto index : peakGroup){
                   cout << index << " ";
               }
               cout << endl;
           }
        }

        int counter = 0;
        for (auto dissimilarity : dissimilarities) {

            if (debug) {
                cout << "***********************" << endl;
                cout << "ITERATION " << counter << endl;
                counter++;
            }

            //refers to index in peakSamplePair (sample id)

            unsigned int firstPeakPairIndex = dissimilarity.second.first;
            unsigned int secondPeakPairIndex = dissimilarity.second.second;

            //refers to index in peakGroups
            int firstContainingClusterIndex = -1;
            int secondContainingClusterIndex = -1;

            //check existing clusters
            for (unsigned int i = 0; i < peakGroups.size(); i++) {

                vector<unsigned int> cluster = peakGroups.at(i);

                if (cluster.empty()) continue;

                for (auto peakPair : cluster){
                    if (peakPair == firstPeakPairIndex){
                        firstContainingClusterIndex = static_cast<int>(i);
                    }
                    if (peakPair == secondPeakPairIndex){
                        secondContainingClusterIndex = static_cast<int>(i);
                    }
                }

                if (firstContainingClusterIndex != -1 && secondContainingClusterIndex != -1) {
                    break;
                }
            }

            if (debug) {
                cout << "FOUND firstContainingClusterIndex=" << firstContainingClusterIndex << ", " << "secondContainingClusterIndex=" << secondContainingClusterIndex << endl;
            }

            /*
             * Based on clusters retrieved, and samples already present, either
             *
             * 1. accept (i,j) pair by merging a peak to an existing cluster.
             *
             * 2. accept (i,j) pair by merging two existing clusters together.
             *
             * 3. accept (i,j) pair by creating a new cluster
             * --> only if i and j are not involved in any other clusters
             *
             * 3. reject (i,j) pair
             * --> if the merge would involve multiple peaks from the same sample joining
             * a cluster, reject the merge
             *
             * In the case of merging, need to update the peakGroups vector appropriately
             *
             * Will probably require an extensive amount of testing
             *
             * meanwhile, is this even the issue? what exactly is happening here?
            */

            //both the first and second peaks are already involved in clusters.
            //If the cluster are different, merge the clusters together.
            //If the clusters are the same, they are already merged together.
            if (firstContainingClusterIndex != -1 && secondContainingClusterIndex != -1 && firstContainingClusterIndex != secondContainingClusterIndex) {

                //retrieve clusters
                vector<unsigned int> firstContainingCluster = peakGroups.at(firstContainingClusterIndex);
                vector<unsigned int> secondContainingCluster = peakGroups.at(secondContainingClusterIndex);

                if (debug) {
                    cout << "MERGE STEP" << endl;

                    cout << "firstContainingCluster: ";
                    for (auto ind : firstContainingCluster){
                        cout << ind << " ";
                    }
                    cout << endl;

                    cout << "secondContainingCluster: ";
                    for (auto ind : secondContainingCluster){
                        cout << ind << " ";
                    }
                    cout << endl;
                }

                //check to see that merging the two clusters would not lead to a cluster with duplicate sample ids.

                vector<unsigned int> intersection;

                vector<unsigned int> firstContainingClusterSamples = vector<unsigned int>(firstContainingCluster.size());
                vector<unsigned int> secondContainingClusterSamples = vector<unsigned int>(secondContainingCluster.size());

                for (unsigned int i = 0; i < firstContainingCluster.size(); i++){
                    firstContainingClusterSamples.at(i) = peakSamplePairs.at(firstContainingCluster.at(i)).first;
                }

                for (unsigned int i = 0; i < secondContainingCluster.size(); i++){
                    secondContainingClusterSamples.at(i) = peakSamplePairs.at(secondContainingCluster.at(i)).first;
                }

                sort(firstContainingClusterSamples.begin(), firstContainingClusterSamples.end());
                sort(secondContainingClusterSamples.begin(), secondContainingClusterSamples.end());

                set_intersection(firstContainingClusterSamples.begin(),
                                 firstContainingClusterSamples.end(),
                                 secondContainingClusterSamples.begin(),
                                 secondContainingClusterSamples.end(),
                                 back_inserter(intersection));

                if (intersection.empty()) {
                    //merge will not lead to multiple peaks from the same sample

                    firstContainingCluster.insert(firstContainingCluster.end(), secondContainingCluster.begin(), secondContainingCluster.end());

                    secondContainingCluster.clear();
                    secondContainingCluster.shrink_to_fit();

                    peakGroups.at(firstContainingClusterIndex) = firstContainingCluster;
                    peakGroups.at(secondContainingClusterIndex) = secondContainingCluster;
                }

                if (debug) {
                    cout << "Updated firstContainingCluster: ";
                    for (auto ind : peakGroups.at(firstContainingClusterIndex)){
                        cout << ind << " ";
                    }
                    cout << endl;

                    cout << "Updated secondContainingCluster: ";
                    for (auto ind: peakGroups.at(secondContainingClusterIndex)){
                        cout << ind << " ";
                    }
                    cout << endl;
                }

            //only the first peak is involved in a cluster already.
            //secondPeakPair joins firstContainingCluster
            } else if (firstContainingClusterIndex != -1 && secondContainingClusterIndex == -1) {

                vector<unsigned int> firstContainingCluster  = peakGroups.at(firstContainingClusterIndex);

                if (debug) {
                    cout << "JOIN FIRST CONTAINING CLUSTER STEP" << endl;

                    cout << "firstContainingCluster: ";
                    for (auto ind : firstContainingCluster){
                        cout << ind << " ";
                    }
                    cout << endl;
                }

                firstContainingCluster.push_back(secondPeakPairIndex);

                peakGroups.at(firstContainingClusterIndex) = firstContainingCluster;

                if (debug) {
                    cout << "Updated firstContainingCluster: ";
                    for (auto ind : peakGroups.at(firstContainingClusterIndex)){
                        cout << ind << " ";
                    }
                    cout << endl;
                }

            //only the second peak is involved in a cluster already.
            //firstPeakPair joins secondContainingCluster
            } else if (firstContainingClusterIndex == -1 && secondContainingClusterIndex != -1) {

                vector<unsigned int> secondContainingCluster = peakGroups.at(secondContainingClusterIndex);

                if (debug) {
                    cout << "JOIN SECOND CONTAINING CLUSTER STEP" << endl;

                    cout << "secondContainingCluster: ";
                    for (auto ind : secondContainingCluster){
                        cout << ind << " ";
                    }
                    cout << endl;
                }

                secondContainingCluster.push_back(firstPeakPairIndex);

                peakGroups.at(secondContainingClusterIndex) = secondContainingCluster;

                if (debug) {
                    cout << "Updated secondContainingCluster: ";
                    for (auto ind: peakGroups.at(secondContainingClusterIndex)){
                        cout << ind << " ";
                    }
                    cout << endl;
                }

            //both the first and second peak are not part of any extant cluster, they merge together to create a new cluster.
            } else if (firstContainingClusterIndex == -1 && secondContainingClusterIndex == -1){

                if (debug) {
                    cout << "NEW CLUSTER STEP" << endl;
                }

                //no existing clusters involving i or j - add a new cluster
                vector<unsigned int> newCluster = {firstPeakPairIndex, secondPeakPairIndex};
                peakGroups.push_back(newCluster);
            }

            if (debug) {
                cout << endl;
                cout << "peakGroups status:" << endl;
                for (auto peakGroup : peakGroups) {
                    if (peakGroup.empty()) continue;
                    cout << "peakGroup Indexes: ";
                    for (auto index : peakGroup){
                        cout << index << " ";
                    }
                    cout << endl;
                }

                cout << "***********************" << endl;
            }
        }

        //Translate results and return
        for (unsigned int i = 0; i < peakGroups.size(); i++){

            if (peakGroups.at(i).empty()) continue;

            PeakGroup grp;
            grp.groupId = static_cast<int>(i);

            for (auto peakPairIndex : peakGroups.at(i)) {
                grp.addPeak(peakSamplePairs.at(peakPairIndex).second);
            }

            grp.groupStatistics();
            pgroups.push_back(grp);
        }

        if (debug) {
            cout << "Returning " << pgroups.size() << " peak groups." << endl;
        }

        return(pgroups);
}

/**
 * @brief EIC::groupPeaks
 * @param eics
 * @param smoothingWindow
 * @param maxRtDiff
 * @return
 *
 *
 * @deprecated in favor of EIC::groupPeaksC()
 */
vector<PeakGroup> EIC::groupPeaks(vector<EIC*>& eics, int smoothingWindow, float maxRtDiff) { 

	//list filled and return by this function
	vector<PeakGroup> pgroups;

	//case there is only a single EIC, there is nothing to group
	if ( eics.size() == 1 && eics[0] != NULL) {
		EIC* m=eics[0];
		for(unsigned int i=0; i< m->peaks.size(); i++ ) {
			PeakGroup grp;
			grp.groupId = i;
			grp.addPeak(m->peaks[i]);
			grp.groupStatistics();
			pgroups.push_back(grp);
		}
		return pgroups;
	}

    //create EIC composed from all sample eics
    EIC* m = EIC::eicMerge(eics);
    if (!m) return pgroups;

    //find peaks in merged eic
    m->getPeakPositions(smoothingWindow);

    //return a.rt < b.rt;
    sort(m->peaks.begin(), m->peaks.end(), Peak::compRt);

	for(unsigned int i=0; i< m->peaks.size(); i++ ) {
		PeakGroup grp;
		grp.groupId = i;
		pgroups.push_back(grp);
  	}

//    cout << "EIC::groupPeaks() eics=" << eics.size() << endl;
//    cout << "EIC::groupPeaks() peakgroups pre-processing=" << pgroups.size() << endl;

	for(unsigned int i=0; i < eics.size(); i++ ) {	//for every sample
		for(unsigned int j=0; j < eics[i]->peaks.size(); j++ ) { //for every peak in the sample
			Peak& b = eics[i]->peaks[j]; 
            b.groupNum=-1;  
            b.groupOverlap=FLT_MIN; 

            vector<Peak>::iterator itr = lower_bound(m->peaks.begin(), m->peaks.end(), b, Peak::compRtMin);
            int lb = (itr-(m->peaks.begin()))-1; if (lb < 0) lb=0;
            //cout << "\tb=" << b.rtmin << "<=>" << b.rtmax << " lb=" << lb << endl;

            //Find best matching group
            for(unsigned int k=lb; k< m->peaks.size(); k++ ) {
                Peak& a = m->peaks[k];
                float overlap = checkOverlap(a.rtmin,a.rtmax,b.rtmin,b.rtmax); //check for overlap
                //cout << "\t\ta=" << a.rtmin << "<=>" << a.rtmax  << " overlap=" << overlap << endl;

                if(overlap == 0 and a.rtmax < b.rtmin) continue;
                if(overlap == 0 and a.rtmin > b.rtmax) break;

                float distx = abs(b.rt-a.rt);
			    if ( distx > maxRtDiff && overlap < 0.2 ) continue;

				float disty = abs(b.peakIntensity-a.peakIntensity);
				//float score= overlap+1/(distx+0.01)+1/(disty+0.01);
				float score = 1.0/(distx+0.01)/(disty+0.01)*overlap;
				//Feng note: the new score function above makes sure that the three terms are weighted equally.
			    if ( score > b.groupOverlap) { b.groupNum=k; b.groupOverlap=score; }
                //cout << "x" << b.rt << " " << b.peakIntensity;
            }

            /*
            cout << b->peakMz <<  " " << b->rtmin << " " << b->rtmax << "->"  << b->groupNum <<
                    " " << b->groupOverlap << endl;
            */


            if (b.groupNum != -1 ) {
                PeakGroup& bestPeakGroup = pgroups[ b.groupNum ];
                bestPeakGroup.addPeak(b);
            } else {
				PeakGroup grp;
				pgroups.push_back(grp);
				grp.groupId = pgroups.size()+1;
				grp.addPeak(b);
				b.groupOverlap=0;
			}
        }
    }

//    cout << "EIC::groupPeaks() peakgroups pre sample cleaning=" << pgroups.size() << endl;

	//clean up peakgroup such that there is only one peak for each sample
    for(unsigned int i=0; i< pgroups.size(); i++) {
            PeakGroup& grp = pgroups[i];
            if (grp.peaks.size() > 1) {
                 grp.reduce();

                 //grp.fillInPeaks(eics);
                 //Feng note: fillInPeaks is unnecessary
                 //Phil note: fillInPeaks should probably be a configurable option

                 //arbitrary order ensures consistency between runs
                 sort(grp.peaks.begin(), grp.peaks.end(), Peak::compSampleName);

                 grp.groupStatistics();
            } else {	//empty group..
                pgroups.erase(pgroups.begin()+i);
                i--;
            }
    }

 //   cout << "EIC::groupPeaks() peakgroups post-processing=" << pgroups.size() << endl;

	//now merge overlapping groups
	//EIC::mergeOverlapingGroups(pgroups);
    //cout << "Found " << pgroups.size() << "groups" << endl;

	if(m) delete(m);
	return(pgroups);
}

/**
 * Group peaks using merged EIC, intensity, RT proximity
 *
 * From MAVEN Issue 381
 *
 * Try to avoid stochasticity/OS-specific issues
 *
 * @brief groupPeaksC
 * @param eics
 * @param smoothingWindow
 * @param maxRtDiff
 * @param baselineSmoothingWindow
 * @param baselineDropTopX
 * @return
 */
vector<PeakGroup> EIC::groupPeaksC(vector<EIC*>& eics, int smoothingWindow, float maxRtDiff, int baselineSmoothingWindow, int baselineDropTopX) {

    //list filled and return by this function
    vector<PeakGroup> pgroups{};

    //case there is only a single EIC, there is nothing to group
    if ( eics.size() == 1 && eics[0]) {
        EIC* m=eics[0];
        for(unsigned int i=0; i< m->peaks.size(); i++ ) {
            PeakGroup grp;
            grp.groupId = static_cast<int>(i);
            grp.addPeak(m->peaks[i]);
            grp.groupStatistics();
            pgroups.push_back(grp);
        }
        return pgroups;
    }

    //create EIC composed from all sample eics
    EIC* m = EIC::eicMerge(eics);
    if (!m) return pgroups;

    m->setBaselineSmoothingWindow(baselineSmoothingWindow);
    m->setBaselineDropTopX(baselineDropTopX);

    //find peaks in merged eic
    m->getPeakPositionsC(smoothingWindow, false, false);

    sort(m->peaks.begin(), m->peaks.end(), Peak::compRt);

//    //debugging
//    cout << "group peak RTs:" << endl;
//    for (auto peak : m->peaks) {
//        cout << to_string(peak.rt) << endl;
//    }
//    cout << endl;

    //m->peaks.pos, sample peaks
    map<int, vector<Peak>> peakGroupData{};

    map<int, set<mzSample*>> peakGroupSamples{};
    for (unsigned int i = 0; i < m->peaks.size(); i++) {
        peakGroupData.insert(make_pair(i, vector<Peak>{}));
        peakGroupSamples.insert(make_pair(i, set<mzSample*>{}));
    }

    vector<Peak> allPeaks{};
    for (auto& eic : eics) {
        for (auto peak : eic->peaks) {
            allPeaks.push_back(peak);
        }
    }

    //try to annotate most intense peaks first
    sort(allPeaks.begin(), allPeaks.end(), [](Peak& lhs, Peak& rhs){
        if (lhs.peakIntensity == rhs.peakIntensity) { //corner case
            return lhs.sample->sampleName < rhs.sample->sampleName; //sample names must differ
        } else {
            return lhs.peakIntensity > rhs.peakIntensity;
        }
    });

    for (auto peak : allPeaks) {

        float peakRt = peak.rt;
        float minPeakRt = max(0.0f, peakRt-maxRtDiff);

        auto lb = lower_bound(m->peaks.begin(), m->peaks.end(), minPeakRt, [](const Peak& p, float rt){
            return p.rt < rt;
        });

        float deltaRt = 999999.0f;
        float bestDeltaRt = deltaRt;
        int bestGroupIndex = -1;

        for (long k = lb - m->peaks.begin(); k < static_cast<long>(m->peaks.size()); k++){

            //corresponds to a merged EIC peak group.
            int groupIndex = static_cast<int>(k);

            //Check if RT for sample query peak is within range of merged eic peak.
            deltaRt = abs(m->peaks[static_cast<unsigned long>(groupIndex)].rt-peakRt);

            //once out of range, always out of range for rest of the loop.
            if (deltaRt > maxRtDiff) break;

            //if a merged eic peak group already contains a peak from the sample query peak's sample,
            //unable to add sample query peak to this merged eic peak group.
            //use mzSample* name to determine uniqueness.
            set<mzSample*> samples = peakGroupSamples[groupIndex];

            if (std::find_if(samples.begin(), samples.end(),[peak](mzSample* sample) {
                         return sample->sampleName == peak.sample->sampleName;}) != samples.end()) continue;

            //only update groupIndex when a closer RT comes along.
            if (deltaRt < bestDeltaRt || bestGroupIndex == -1) {
                bestDeltaRt = deltaRt;
                bestGroupIndex = static_cast<int>(groupIndex);
            }

        }

        if (bestGroupIndex != -1) {
            peakGroupSamples[bestGroupIndex].insert(peak.sample);
            peakGroupData[bestGroupIndex].push_back(peak);
        }

    }

    for (auto it = peakGroupData.begin(); it != peakGroupData.end(); ++it) {
        int groupIndex = it->first;
        vector<Peak> peaks = it->second;

        if (peaks.empty()) continue;

        PeakGroup grp;
        grp.groupId = groupIndex;
        for (auto peak : peaks) {
            grp.addPeak(peak);
        }
        sort(grp.peaks.begin(), grp.peaks.end(), Peak::compSampleName);

        grp.groupStatistics();

        pgroups.push_back(grp);
    }

    if (m) delete(m);
    return pgroups;
}

/**
 * @brief groupPeaksD
 * Group peaks using merged EIC, intensity, RT proximity
 *
 * From mass_spec Issue #692 / MAVEN Issue #482
 *
 * Combine similar peak groups together, when appropriate.
 * This replaces peakdetector method reduceGroups()
 *
 * Comment: This approach will select the most intense peak within RT tolerance for a given sample,
 * and associate this peak with the peak group.
 * This means that any less intense peaks from the same sample that have closer RT will be excluded.
 * This tends to happen when the maxRtDiff is large enough to permit multiple peaks within tolerance
 * to match to a given peakgroup RT, which, in practice, are usually noise peaks.
 * A series of quant approaches are compared. If none of these are different, the sample name is considered,
 * and finally, if multiple peaks from the same sample have identical quant, the RT value is used as an absolute tiebreaker.
 *
 * In some cases, it may be preferable to prefer RT over every measure of intensity.  In that case,
 * a different algorithm should be used.
 *
 * @param eics
 * @param smoothingWindow
 * @param maxRtDiff
 * @param baselineSmoothingWindow
 * @param baselineDropTopX
 * @param mergeOverlap
 * @return
 */
vector<PeakGroup> EIC::groupPeaksD(vector<EIC*>& eics, int smoothingWindow, float maxRtDiff, int baselineSmoothingWindow, int baselineDropTopX, float mergeOverlap, bool debug) {

    shared_ptr<PeakPickingAndGroupingParameters> params = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());

    params->mergedBaselineSmoothingWindow = baselineSmoothingWindow;
    params->mergedBaselineDropTopX = baselineDropTopX;
    params->mergedSmoothingWindow = smoothingWindow;
    params->mergedIsComputeBounds = false;
    params->mergedPeakRtBoundsMaxIntensityFraction = -1.0f;
    params->mergedPeakRtBoundsSlopeThreshold = -1.0f;
    params->groupMaxRtDiff = maxRtDiff;
    params->groupMergeOverlap = mergeOverlap;

    return groupPeaksE(eics, params, debug);
}

vector<PeakGroup> EIC::groupPeaksE(vector<EIC*>& eics, shared_ptr<PeakPickingAndGroupingParameters> params, bool debug){

    //list filled and return by this function
    vector<PeakGroup> pgroups{};

    //case there is only a single EIC, there is nothing to group
    if ( eics.size() == 1 && eics[0]) {
        EIC* m=eics[0];
        for(unsigned int i=0; i< m->peaks.size(); i++ ) {
            PeakGroup grp;
            grp.groupId = static_cast<int>(i);
            grp.addPeak(m->peaks[i]);
            grp.groupStatistics();
            pgroups.push_back(grp);
        }
        return pgroups;
    }

    //create EIC composed from all sample eics
    EIC* m = EIC::eicMerge(eics);
    if (!m) return pgroups;

    shared_ptr<PeakPickingAndGroupingParameters> mergedEICParams =
        PeakPickingAndGroupingParameters::getMergedAsPeakParams(params);

    m->setBaselineSmoothingWindow(params->mergedBaselineSmoothingWindow);
    m->setBaselineDropTopX(params->mergedBaselineDropTopX);

    //find peaks in merged eic
    m->getPeakPositionsD(mergedEICParams, debug);

    //Issue 597: Remove peaks with insufficient ratios
    if (params->mergedSmoothedMaxToBoundsMinRatio > 0) {

        if (debug) cout << "mergedSmoothedMaxToBoundsMinRatio filtering applied." << endl;

        vector<Peak> passingPeaks{};
        for (auto p : m->peaks) {

            float leftSmoothedIntensity = m->spline[p.minpos];
            float rightSmoothedIntensity = m->spline[p.maxpos];
            float maxSmoothedIntensity = m->spline[p.pos];

            if (debug) cout << "max: " << maxSmoothedIntensity << ", left=" << leftSmoothedIntensity << ", right=" << rightSmoothedIntensity;

            float boundsIntensity;
            if (params->mergedSmoothedMaxToBoundsIntensityPolicy == SmoothedMaxToBoundsIntensityPolicy::MINIMUM) {
                boundsIntensity = min(leftSmoothedIntensity, rightSmoothedIntensity);
            } else if (params->mergedSmoothedMaxToBoundsIntensityPolicy == SmoothedMaxToBoundsIntensityPolicy::MAXIMUM) {
                boundsIntensity = max(leftSmoothedIntensity, rightSmoothedIntensity);
            } else if (params->mergedSmoothedMaxToBoundsIntensityPolicy == SmoothedMaxToBoundsIntensityPolicy::MEDIAN) {
                boundsIntensity = 0.5f * (leftSmoothedIntensity + rightSmoothedIntensity);
            } else {
                cerr << "Unhandled SmoothedMaxToBoundsIntensityPolicy policy in EIC::groupPeaksE()! exiting." << endl;
                abort();
            }

            if (debug) cout << ", boundsIntensity=" << boundsIntensity;

            if (boundsIntensity <= 0.0f) {
                passingPeaks.push_back(p);

                if (debug) cout << ", peak retained.";

            } else {

                float smoothedBoundsRatio = maxSmoothedIntensity / boundsIntensity;

                if (debug) cout << ", ratio=" << smoothedBoundsRatio;

                if (smoothedBoundsRatio >= params->mergedSmoothedMaxToBoundsMinRatio) {
                    passingPeaks.push_back(p);

                    if (debug) cout << ", peak retained.";
                }
            }

            if (debug) cout << endl;
        }
        m->peaks = passingPeaks;
    }

    //RT values should never be identical for merged EIC.
    sort(m->peaks.begin(), m->peaks.end(), [](Peak& lhs, Peak& rhs){
        return lhs.rt < rhs.rt;
    });

    if (debug) {
        for (unsigned int i = 1; i < m->peaks.size(); i++) {
            if (m->peaks.at(i).rt == m->peaks.at(i-1).rt) {
                cerr << "Merged EIC Peaks found with identical RT!";
                cerr << "RT="<< m->peaks.at(i).rt << endl;
                abort();
            }
        }
    }

    //calls PeakGroups::groupStatistics()
    pgroups = mergedEICToGroups(eics, m, params->groupMaxRtDiff, params->groupMergeOverlap, debug);

    if (m) delete(m);
    return pgroups;
}

vector<PeakGroup> EIC::mergedEICToGroups(vector<EIC*>& eics, EIC* m, float groupMaxRtDiff, float groupMergeOverlap, bool debug){

    vector<PeakGroup> pgroups{};

    if (!m) return pgroups;

    //m->peaks.pos, sample peaks
    map<int, PeakContainer> peakGroupData{};

    map<int, set<mzSample*>> peakGroupSamples{};
    for (unsigned int i = 0; i < m->peaks.size(); i++) {
        peakGroupData.insert(make_pair(i, PeakContainer()));
        peakGroupSamples.insert(make_pair(i, set<mzSample*>{}));
    }

    vector<Peak> allPeaks{};
    for (auto& eic : eics) {
        for (auto peak : eic->peaks) {
            allPeaks.push_back(peak);
        }
    }

    if (debug) {
       cout << "EIC::mergedEICToGroups(): START allPeaks" << endl;
       for (auto peak : allPeaks) {
           cout << fixed << setprecision(5)
                << "(" << peak.peakMz << ", " << peak.rt << ", " << peak.peakIntensity << ") " << peak.sample->sampleName
                << endl;
       }
       cout << "EIC::mergedEICToGroups(): END allPeaks" << endl;
    }

    //try to annotate most intense peaks first
    // fall back to other metrics of intensity
    // use, samples, rt to break ties
    sort(allPeaks.begin(), allPeaks.end(), [](Peak& lhs, Peak& rhs){

        if (lhs.peakIntensity != rhs.peakIntensity) {
            return lhs.peakIntensity > rhs.peakIntensity;
        }

        if (lhs.peakAreaTop != rhs.peakAreaTop) {
            return lhs.peakAreaTop > rhs.peakAreaTop;
        }

        if (lhs.peakAreaCorrected != rhs.peakAreaCorrected) {
            return lhs.peakAreaCorrected > rhs.peakAreaCorrected;
        }

        if (lhs.peakArea != rhs.peakArea) {
            return lhs.peakArea > rhs.peakArea;
        }

        //sample name is arbitrary tiebreaker if every measure of quant is equal
        if (lhs.sample != rhs.sample) {
            return lhs.sample < rhs.sample;
        }

        //rt is last resort tiebreaker to produce deterministic ordering
        return lhs.rt < rhs.rt;
    });

    for (auto peak : allPeaks) {

        float peakRt = peak.rt;
        float minPeakRt = max(0.0f, peakRt-groupMaxRtDiff);

        auto lb = lower_bound(m->peaks.begin(), m->peaks.end(), minPeakRt, [](const Peak& p, float rt){
            return p.rt < rt;
        });

        float deltaRt = 999999.0f;
        float bestDeltaRt = deltaRt;
        int bestGroupIndex = -1;

        //TODO: How to handle cases where there are multiple peaks from the same sample
        //with identical intensities, one should be associated with a peak group, while another gets associated with no peak?
        //Currently, this is peak-centric assignment.
        //Probably, instead of skipping over cases where the peak is already associated,
        //need to record multiple associations, and take only the one that is closest to the group Rt.

        for (long k = lb - m->peaks.begin(); k < static_cast<long>(m->peaks.size()); k++){

            //corresponds to a merged EIC peak group.
            int groupIndex = static_cast<int>(k);

            //Check if RT for sample query peak is within range of merged eic peak.
            deltaRt = abs(m->peaks[static_cast<unsigned long>(groupIndex)].rt-peakRt);

            //once out of range, always out of range for rest of the loop.
            if (deltaRt > groupMaxRtDiff) break;

            //if a merged eic peak group already contains a peak from the sample query peak's sample,
            //unable to add sample query peak to this merged eic peak group.
            //use mzSample* name to determine uniqueness.
            set<mzSample*> samples = peakGroupSamples[groupIndex];

            if (std::find_if(samples.begin(), samples.end(),[peak](mzSample* sample) {
                         return sample->sampleName == peak.sample->sampleName;}) != samples.end()) continue;

            //only update groupIndex when a closer RT comes along.
            if (deltaRt < bestDeltaRt || bestGroupIndex == -1) {
                bestDeltaRt = deltaRt;
                bestGroupIndex = static_cast<int>(groupIndex);
            }

        }

        if (bestGroupIndex != -1) {
            peakGroupSamples[bestGroupIndex].insert(peak.sample);
            peakGroupData[bestGroupIndex].peaks.insert(make_pair(peak.sample, peak));
        }

    }

    for (auto& it : peakGroupData) {
        it.second.mergedEICPeakIndexes.insert(it.first);
        it.second.recomputeProperties();
    }

    if (debug) {
        cout << "# allPeaks: " << allPeaks.size()
             << ", # peakGroupData: " << peakGroupData.size()
             << endl;
    }

    unsigned long iterationCounter = 0;

    //Progressively merge peaks until group overlap issues are resolved
    while (true) {

        if (debug) {
            cout << "ITERATION #" << iterationCounter << ":" << endl;
        }

        IntegerSetContainer merges;

        //initialize set container
        for (unsigned int i = 0; i < m->peaks.size(); i++) {
            merges.containerBySet.insert(make_pair(i, set<int>{static_cast<int>(i)}));
        }

        if (debug) {
            cout << "# merges.containerBySet(): " << merges.containerBySet.size() << endl;
        }

        //merge data based on existing containers
        for (unsigned int i = 0; i < m->peaks.size(); i++) {

            PeakContainer peaksI = peakGroupData[static_cast<int>(i)];

            if (peaksI.peaks.empty()) continue;

            float minPeakRtI = peaksI.minPeakRt;
            float maxPeakRtI = peaksI.maxPeakRt;

            for (unsigned int j = i+1; j < m->peaks.size(); j++) {

                PeakContainer peaksJ = peakGroupData[static_cast<int>(j)];

                if (peaksJ.peaks.empty()) continue;

                float minPeakRtJ = peaksJ.minPeakRt;
                float maxPeakRtJ = peaksJ.maxPeakRt;

                //Implies an overlap of 0, and all subsequent peaks should have an overlap of 0 compared to peaksI.
                if (minPeakRtJ > maxPeakRtI) break;

                float rtOverlap = mzUtils::checkOverlap(minPeakRtI, maxPeakRtI, minPeakRtJ, maxPeakRtJ);

                if (rtOverlap >= groupMergeOverlap) {
                    merges.addMerge((make_pair(i, j)));
                }
            }
        }

        merges.combineContainers(debug);

        //stay in the while loop until no more merges need to be made.
        if (merges.isAllContainersSize(1)) {
            if (debug) cout << "Peak group convergence reached." << endl;
            break;

        //recreate peakGroupData using the containers indicated by merges.
        } else {
            map<int, PeakContainer> updatedPeakGroupData{};

            for (unsigned int i = 0; i < m->peaks.size(); i++) {

                int intKey = static_cast<int>(i);

                //ensure that every peak index is reassociated in the map
                if (updatedPeakGroupData.find(intKey) == updatedPeakGroupData.end()) {

                    //Case 1: create a new peak container from the merged data
                    if (merges.containerBySet.find(intKey) != merges.containerBySet.end() && merges.containerBySet[intKey].size() > 1) {

                        set<int> container = merges.containerBySet[intKey];

                        PeakContainer mergedContainer;
                        for (auto index : container) {
                            mergedContainer.mergePeakContainer(peakGroupData[index]);
                        }
                        mergedContainer.recomputeProperties();

                        bool isAddedMerged = false;
                        for (auto index : container) {
                            if (isAddedMerged) {
                                updatedPeakGroupData.insert(make_pair(index, PeakContainer()));
                            } else {
                                updatedPeakGroupData.insert(make_pair(index, mergedContainer));
                                isAddedMerged = true;
                            }
                        }

                    //Case 2: peak container is unchanged from previous set
                    } else {
                        updatedPeakGroupData.insert(make_pair(i, peakGroupData[intKey]));
                    }

                }
            }

            iterationCounter++;
            peakGroupData = updatedPeakGroupData;
        }
    }

    for (auto it = peakGroupData.begin(); it != peakGroupData.end(); ++it) {

        int groupIndex = it->first;
        map<mzSample*, Peak> peaks = it->second.peaks;

        if (debug) {
            cout << "groupIndex=" << groupIndex << ", # Peaks=" << peaks.size() << endl;
        }

        if (peaks.empty()) continue;

        PeakGroup grp;
        grp.groupId = groupIndex;
        for (auto& peak : peaks) {
            grp.addPeak(peak.second);
        }
        sort(grp.peaks.begin(), grp.peaks.end(), Peak::compSampleName);

        grp.groupStatistics();

        if (debug) {
            grp.summary();
        }

        BlankSingleIntensities intensities = EIC::calculateBlankBackground(eics, grp.minRt, grp.maxRt, debug);

        grp.blankMaxHeight = intensities.maxSingleIntensity;
        grp.blankMedianHeight = intensities.medianSingleIntensity;
        grp.mergedEICSummaryData = EIC::calculateMergedEICSummaryData(m, it->second.mergedEICPeakIndexes, debug);
        grp.maxBlankRawSignal = EIC::calculateMaxBlankSignalBackground(m, eics, it->second.mergedEICPeakIndexes, false, debug);
        grp.maxBlankSmoothedSignal = EIC::calculateMaxBlankSignalBackground(m, eics, it->second.mergedEICPeakIndexes, true, debug);

        pgroups.push_back(grp);
    }

    return pgroups;
}


map<int, PeakContainer> EIC::mergePeakContainers(map<int, PeakContainer> peakGroupData, float groupMergeOverlap, bool debug) {

    unsigned long iterationCounter = 0;
    unsigned int N = peakGroupData.size();

    //Progressively merge peaks until group overlap issues are resolved
    while (true) {

        if (debug) {
            cout << "ITERATION #" << iterationCounter << ":" << endl;
        }

        IntegerSetContainer merges;

        //initialize set container
        for (unsigned int i = 0; i < N; i++) {
            merges.containerBySet.insert(make_pair(i, set<int>{static_cast<int>(i)}));
        }

        if (debug) {
            cout << "# merges.containerBySet(): " << merges.containerBySet.size() << endl;
        }

        //merge data based on existing containers
        for (unsigned int i = 0; i < N; i++) {

            PeakContainer peaksI = peakGroupData[static_cast<int>(i)];

            if (peaksI.peaks.empty()) continue;

            float minPeakRtI = peaksI.minPeakRt;
            float maxPeakRtI = peaksI.maxPeakRt;

            for (unsigned int j = i+1; j < N; j++) {

                PeakContainer peaksJ = peakGroupData[static_cast<int>(j)];

                if (peaksJ.peaks.empty()) continue;

                float minPeakRtJ = peaksJ.minPeakRt;
                float maxPeakRtJ = peaksJ.maxPeakRt;

                //Implies an overlap of 0, and all subsequent peaks should have an overlap of 0 compared to peaksI.
                if (minPeakRtJ > maxPeakRtI) break;

                float rtOverlap = mzUtils::checkOverlap(minPeakRtI, maxPeakRtI, minPeakRtJ, maxPeakRtJ);

                if (rtOverlap >= groupMergeOverlap) {
                    merges.addMerge((make_pair(i, j)));
                }
            }
        }

        merges.combineContainers(debug);

        //stay in the while loop until no more merges need to be made.
        if (merges.isAllContainersSize(1)) {
            if (debug) cout << "Peak group convergence reached." << endl;
            break;

            //recreate peakGroupData using the containers indicated by merges.
        } else {
            map<int, PeakContainer> updatedPeakGroupData{};

            for (unsigned int i = 0; i < N; i++) {

                int intKey = static_cast<int>(i);

                //ensure that every peak index is reassociated in the map
                if (updatedPeakGroupData.find(intKey) == updatedPeakGroupData.end()) {

                    //Case 1: create a new peak container from the merged data
                    if (merges.containerBySet.find(intKey) != merges.containerBySet.end() && merges.containerBySet[intKey].size() > 1) {

                        set<int> container = merges.containerBySet[intKey];

                        PeakContainer mergedContainer;
                        for (auto index : container) {
                            mergedContainer.mergePeakContainer(peakGroupData[index]);
                        }
                        mergedContainer.recomputeProperties();

                        bool isAddedMerged = false;
                        for (auto index : container) {
                            if (isAddedMerged) {
                                updatedPeakGroupData.insert(make_pair(index, PeakContainer()));
                            } else {
                                updatedPeakGroupData.insert(make_pair(index, mergedContainer));
                                isAddedMerged = true;
                            }
                        }

                        //Case 2: peak container is unchanged from previous set
                    } else {
                        updatedPeakGroupData.insert(make_pair(i, peakGroupData[intKey]));
                    }

                }
            }

            iterationCounter++;
            peakGroupData = updatedPeakGroupData;
        }
    }

    return peakGroupData;
}

void EIC::interpolate() {

    unsigned int lastNonZero=0;
    for(unsigned int posi=0; posi < intensity.size(); posi++ ) {

        if (intensity[posi] != 0 ) {
             lastNonZero=posi;  //if this position has nonzero intensity, mark it as lastNonZero position
        }

        if (intensity[posi] == 0 and lastNonZero > 0) { //interpolate
            unsigned int nextNonZero=0;
            for(unsigned int j=posi; j<intensity.size();j++) { if (intensity[j] != 0 ) { nextNonZero=j; }}
            if( nextNonZero == 0) continue;

            //start at first empty position and until next non empty
            for(unsigned int j=posi; j<nextNonZero; j++)  {
                float fracDist = (j-lastNonZero)/ (float) (nextNonZero-lastNonZero);
                float newIntensity = intensity[lastNonZero] + fracDist*intensity[nextNonZero];
                intensity[j] = newIntensity;
                lastNonZero=j; posi++;
            }
        }
    }
}

/*
void EIC::cubicSplineFit()  {
	unsigned int n = size();
	float* x = new float[n];
	float* f = new float[n];
	float* b = new float[n];
	float* c = new float[n];
	float* d = new float[n];

	int N=0;
	for(int j=0; j<n; j++) {
		x[j]=rt[j];
		f[j]=intensity[j];
		b[j]=c[j]=d[j]=0; //init all elements to 0

		if(spline[j]>baseline[j] and intensity[j]>0) {
			x[N]=rt[j]; f[N]=intensity[j];
			N++;
		} else if (spline[j] <= baseline[j]*1.1) {
			x[N]=rt[j]; f[N]=baseline[j];
			N++;
		}
	}

	if(N <= 2) continue;
	mzUtils::cubic_nak(N,x,f,b,c,d);

	for(int j=1; j<N; j++) {
		float rtstep = (x[j]-x[j-1])/10;
		for(int k=0; k<10; k++) {
			float dt = rtstep*k;
			float y = f[j-1] + ( dt ) * ( b[j-1] + ( dt ) * ( c[j-1] + (dt) * d[j-1] ) );
			//float y = mzUtils::spline_eval(n,x,f,b,c,d,x[j]+dt);
			if(y < 0) y= 0;
		}
	}
	delete[] x;
	delete[] f;
	delete[] b;
	delete[] c;
	delete[] d;
}
*/



vector<Scan*> EIC::getFragmentationEvents() {
    vector<Scan*>matchedscans;
	if(!sample) return matchedscans;

    for( unsigned int j=0; j < sample->scans.size(); j++ ) {
            Scan* scan = sample->scans[j];
            if (!scan or scan->mslevel <= 1 or scan->rt < rtmin) continue; //skip ms1 events
            if (scan->rt > rtmax) break;
            if (scan->precursorMz >= mzmin and scan->precursorMz <= mzmax) {
                matchedscans.push_back(scan);
            }
     }
    return matchedscans;
}   

void EIC::removeOverlapingPeaks() { 
    for(unsigned int i=0; i <peaks.size(); i++ ) peaks[i].localMaxFlag=true;

    for(unsigned int i=0; i < peaks.size(); i++ ) {
		Peak& a = peaks[i];
		if (a.peakIntensity == 0) { a.localMaxFlag=false; continue; }

        for(unsigned int j=i+1; j < peaks.size(); j++ ) {
			Peak& b = peaks[j];
			if ( mzUtils::checkOverlap(a.rtmin,a.rtmax,b.rtmin,b.rtmax)> 0.95) { //overlap
				if (a.peakArea >= b.peakArea) { b.localMaxFlag=false;}
			}
		}
	}
	vector<Peak>reduced;
    for(unsigned int i=0; i < peaks.size(); i++ ) {
		if (peaks[i].localMaxFlag == true) reduced.push_back(peaks[i]);
	}
	peaks = reduced;
}

PeakGroupBaseline EIC::calculateMergedEICSummaryData(EIC* mergedEIC, set<int> mergedEICPeakIndexes, bool debug) {
    PeakGroupBaseline mergedEICSummaryData;

    if (mergedEIC) {

        //Find representative index
        int representativeIndex = -1;
        float representativeIntensity = -1.0f;

        for (auto peakIndex : mergedEICPeakIndexes) {
            Peak p = mergedEIC->peaks.at(peakIndex);
            if (mergedEIC->spline[p.pos] > representativeIntensity) {
                representativeIntensity = mergedEIC->spline[p.pos];
                representativeIndex = peakIndex;
            }
        }

        if (debug) cout << "EIC::calculateMergedEICSummaryData(): peakIndex=" << representativeIndex << endl;

        if (representativeIndex != -1) {
            // compute values based on representative index
            Peak p = mergedEIC->peaks.at(representativeIndex);

            if (debug) {
                cout << "EIC::calculateMergedEICSummaryData(): "
                     << "( " << p.minpos << " - ["
                     << p.minPosFWHM << " - {" << (p.pos-1) << " " << p.pos << " " << (p.pos+1) << "} - "
                     << p.maxPosFWHM << " ] - " << p.maxpos << " )"
                     << endl;
            }

            for (auto i = p.minpos; i <= p.maxpos; i++) {

                mergedEICSummaryData.fullRangeBaseline += mergedEIC->baseline[i];

                //Only compute FWHM baseline if the FWHM could be computed.
                if (p.maxPosFWHM > p.minPosFWHM && i >= p.minPosFWHM && i <= p.maxPosFWHM) {
                    mergedEICSummaryData.FWHMBaseline += mergedEIC->baseline[i];
                }

                if (i >= p.pos-1 && i <= p.pos+1) {
                    mergedEICSummaryData.threePointBaseline += mergedEIC->baseline[i];
                }

                if (i == p.pos) {
                    mergedEICSummaryData.pickedPeakBaseline = mergedEIC->baseline[i];
                }
            }

            if (debug) {
                cout << "EIC::calculateMergedEICSummaryData(): "
                     << " FULL=" << mergedEICSummaryData.fullRangeBaseline
                     << " FWHM=" << mergedEICSummaryData.FWHMBaseline
                     << " THREE=" << mergedEICSummaryData.threePointBaseline
                     << " ONE=" << mergedEICSummaryData.pickedPeakBaseline
                     << endl;
            }
        }
    }

    return mergedEICSummaryData;
}

BlankSingleIntensities EIC::calculateBlankBackground(vector<EIC *>& eics, float rtMin, float rtMax, bool debug){
    if (debug) cout << "EIC::calculateBlankBackground() started." << endl;

    vector<float> sampleMaxIntensities{};

    for (auto eic : eics) {
        if (eic->sample->isBlank) {

            float sampleMaxIntensity = 0;

            for (unsigned int i = 0; i < eic->size(); i++) {
                if (eic->rt[i] >= rtMin && eic->rt[i] <= rtMax) {
                    if (eic->intensity[i] > sampleMaxIntensity) {
                        sampleMaxIntensity = eic->intensity[i];
                    }
                }
            }

            sampleMaxIntensities.push_back(sampleMaxIntensity);
        }
    }

    //Issue 678: Handle case where there are no blank samples in the sample set.
    float blankMaxIntensity = 0.0f;
    if (!sampleMaxIntensities.empty()){
        blankMaxIntensity = *max_element(sampleMaxIntensities.begin(), sampleMaxIntensities.end());
    }

    BlankSingleIntensities blankSingleIntensities;
    blankSingleIntensities.maxSingleIntensity = blankMaxIntensity;
    blankSingleIntensities.medianSingleIntensity = mzUtils::median(sampleMaxIntensities);

    if (debug) cout << "EIC::calculateBlankBackground() finished." << endl;
    return blankSingleIntensities;
}

float EIC::getAnalogousIntensitySum(EIC* eic, float rtAnchor, unsigned int numPoints, bool isUseSmoothedIntensity, bool debug) {
    if (debug) cout << "EIC::getAnalogousIntensitySum() started." << endl;

    float intensitySum = 0.0f;

    if (!eic || eic->size() == 0) return intensitySum;

    if (debug){
        cout << "EIC::getAnalogousIntensitySum(): "
             << " eic->size()=" << eic->size()
             << " spline.size()=" << eic->spline.size()
             << " rt.size()=" << eic->rt.size()
             << endl;
    }
    //find closest point in rt to rtAnchor
    auto lb = lower_bound(eic->rt.begin(), eic->rt.end(), rtAnchor);
    auto pos = lb - eic->rt.begin();

    //lower_bound returns first value not less than the rtAnchor,
    //check the position before the lower_bound is closer in RT to rtAnchor than the lower_bound.
    auto rtPos = (pos > 0 && abs(rtAnchor - eic->rt[pos-1]) < abs(rtAnchor - eic->rt[pos])) ? (pos-1) : pos;

    if (debug) {
        cout << "EIC::getAnalogousIntensitySum(): "
             << "rtPos=" << rtPos
             << ", eic->rt[rtPos]=" << eic->rt[rtPos]
             << endl;
    }

    intensitySum += isUseSmoothedIntensity ? eic->spline[rtPos] : eic->intensity[rtPos];

    if (debug) {
        cout << "EIC::getAnalogousIntensitySum(): "
             << "Initial intensitySum=" << intensitySum
             << endl;
    }

    unsigned int pointsAroundMax = numPoints-1;

    if (pointsAroundMax % 2 == 1){
        pointsAroundMax = pointsAroundMax + 1;
    }

    pointsAroundMax /= 2;

    if (debug) {
        cout << "EIC::getAnalogousIntensitySum(RT = " << rtAnchor << ", " << numPoints << " points):"
             << " pointsAroundMax= " << pointsAroundMax << ", eic->size()=" << eic->size() << "\n"
             << "rt[i-2]: " << "rt[" << (rtPos-2) << "]: " << eic->rt[rtPos-2] << ", diff = " << abs(eic->rt[rtPos-2]-rtAnchor) << "\n"
             << "rt[i-1]: " << "rt[" << (rtPos-1) << "]: " << eic->rt[rtPos-1] << ", diff = " << abs(eic->rt[rtPos-1]-rtAnchor) << "\n"
             << "** rt[i]: " << "rt[" << (rtPos) << "]: " <<  eic->rt[rtPos] << ", diff = " << abs(eic->rt[rtPos]-rtAnchor) << " **\n"
             << "rt[i+1]: " << "rt[" << (rtPos+1) << "]: " <<  eic->rt[rtPos+1] << ", diff = " << abs(eic->rt[rtPos+1]-rtAnchor) << "\n"
             << "rt[i+2]: " << "rt[" << (rtPos+2) << "]: " <<  eic->rt[rtPos+2] << ", diff = " << abs(eic->rt[rtPos+2]-rtAnchor) << "\n"
             << endl;
    }

    unsigned int pointsLeft = 0;
    while (pointsLeft < pointsAroundMax) {
        pointsLeft++;
        if (rtPos-pointsLeft < 0) break;
        intensitySum += isUseSmoothedIntensity ? eic->spline[rtPos-pointsLeft] : eic->intensity[rtPos-pointsLeft];

        if (debug) {
            cout << "EIC::getAnalogousIntensitySum() left:"
                 <<" pointsLeft=" << pointsLeft
                 << ", pointsAroundMax=" << pointsAroundMax
                 << ", intensitySum=" << intensitySum
                 << endl;
        }
    }

    unsigned int pointsRight = 0;
    while (pointsRight < pointsAroundMax) {
        pointsRight++;
        if (rtPos+pointsRight > eic->size()-1) break;
        intensitySum += isUseSmoothedIntensity ? eic->spline[rtPos+pointsRight] : eic->intensity[rtPos+pointsRight];

        if (debug) {
            cout << "EIC::getAnalogousIntensitySum() left:"
                 <<" pointsLeft=" << pointsRight
                 << ", pointsAroundMax=" << pointsAroundMax
                 << ", intensitySum=" << intensitySum
                 << endl;
        }
    }

    if (debug) cout << "EIC::getAnalogousIntensitySum(): final intensitySum=" << intensitySum << endl;

    return intensitySum;
}

/**
 * Issue 679: deprecated in favor of Peakgroup::getMaxBlankCorrespondingQuant
 * @deprecated
 */
PeakGroupBaseline EIC::calculateMaxBlankSignalBackground(
    EIC* mergedEIC,
    vector<EIC *>& eics,
    set<int> mergedEICPeakIndexes,
    bool isUseSmoothedIntensity,
    bool debug) {
    if (debug) cout << "EIC::calculateMaxBlankSignalBackground() started." << endl;

    PeakGroupBaseline maxBlankSignalBackground;

    if (mergedEIC) {

        //Find representative index
        int representativeIndex = -1;
        float representativeIntensity = -1.0f;

        if (debug) {
            cout << "EIC::calculateMaxBlankSignalBackground():"
                 << " Peaks: {";
            for (unsigned int i = 0; i < mergedEIC->peaks.size(); i++) {
                if (i > 0) cout << ", ";
                cout << mergedEIC->peaks[i].pos;
            }
            cout << "}" << endl;
        }

        for (auto peakIndex : mergedEICPeakIndexes) {
            Peak p = mergedEIC->peaks.at(peakIndex);
            if (debug) {
                cout << "EIC::calculateMaxBlankSignalBackground():"
                     << " mergedEICPeakIndex=" << peakIndex << ", pos=" << p.pos
                     << " (max pos = " << (mergedEIC->size()-1) << ")"
                     << endl;
            }
            if (p.pos < 0 || p.pos > (mergedEIC->size()-1)) continue;
            if (mergedEIC->spline[p.pos] > representativeIntensity) {
                representativeIntensity = mergedEIC->spline[p.pos];
                representativeIndex = peakIndex;
            }
        }

        if (debug) cout << "EIC::calculateMaxBlankSignalBackground(): peakIndex=" << representativeIndex << endl;

        if (representativeIndex != -1) {
            // compute values based on representative index
            Peak p = mergedEIC->peaks.at(representativeIndex);

            if (debug) {
                cout << "EIC::calculateMaxBlankSignalBackground(): "
                     << "( " << p.minpos << " - ["
                     << p.minPosFWHM << " - {" << (p.pos-1) << " " << p.pos << " " << (p.pos+1) << "} - "
                     << p.maxPosFWHM << "] - " << p.maxpos << " )"
                     << endl;
            }

            float maxFullRangeBaseline = 0.0f;
            float maxFWHMBaseline = 0.0f;
            float maxThreePointBaseline = 0.0f;
            float maxOnePointBaseline = 0.0f;

            for (auto eic : eics) {
                if (eic->sample && eic->sample->isBlank) {
                    float fullRangeBaseline = EIC::getAnalogousIntensitySum(eic, p.rt, (p.maxpos-p.minpos+1), isUseSmoothedIntensity, debug);

                    //Issue 668: Only try to compute a FWHM baseline if the FWHM could be computed for the merged EIC
                    float fwhmBaseline = 0.0f;
                    if (p.maxPosFWHM > p.minPosFWHM) {
                        fwhmBaseline = EIC::getAnalogousIntensitySum(eic, p.rt, (p.maxPosFWHM-p.minPosFWHM+1), isUseSmoothedIntensity, debug);
                    }

                    float threePointBaseline = EIC::getAnalogousIntensitySum(eic, p.rt, 3, isUseSmoothedIntensity, debug);
                    float onePointBaseline = EIC::getAnalogousIntensitySum(eic, p.rt, 1, isUseSmoothedIntensity, debug);

                    if (fullRangeBaseline > maxFullRangeBaseline) maxFullRangeBaseline = fullRangeBaseline;
                    if (fwhmBaseline > maxFWHMBaseline) maxFWHMBaseline = fwhmBaseline;
                    if (threePointBaseline > maxThreePointBaseline) maxThreePointBaseline = threePointBaseline;
                    if (onePointBaseline > maxOnePointBaseline) maxOnePointBaseline = onePointBaseline;
                }
            }

            maxBlankSignalBackground.fullRangeBaseline = maxFullRangeBaseline;
            maxBlankSignalBackground.FWHMBaseline = maxFWHMBaseline;
            maxBlankSignalBackground.threePointBaseline = maxThreePointBaseline;
            maxBlankSignalBackground.pickedPeakBaseline = maxOnePointBaseline;

            if (debug) {
                cout << "EIC::calculateMaxBlankSignalBackground(): "
                     << " FULL=" << maxBlankSignalBackground.fullRangeBaseline
                     << " FWHM=" << maxBlankSignalBackground.FWHMBaseline
                     << " THREE=" << maxBlankSignalBackground.threePointBaseline
                     << " ONE=" << maxBlankSignalBackground.pickedPeakBaseline
                     << endl;
            }
        }
    }

    return maxBlankSignalBackground;
}

// Useful for debugging
void EIC::encodeToFile(vector<EIC*> eics, string filePath){

    ofstream stream(filePath, ios::out);

    stream << "START_SET" << endl;
    stream << "NUM_EICS " << eics.size() << endl;

    for (auto eic : eics) {
        stream << "START_EIC" << endl;
        stream << eic->getSample()->sampleName << endl;
        stream << "START_COORDS" << endl;
        stream << "NUM_COORDS " << eic->size() << endl;
        for (unsigned int i = 0; i < eic->size(); i++) {
            stream << eic->mz[i] << "\t"
                   << eic->rt[i] << "\t"
                   << eic->intensity[i] << "\t"
                   << eic->spline[i] << "\t"
                   << eic->baseline[i]
                   << endl;
        }
        stream << "END_COORDS" << endl;
        stream << "START_PEAKS" << endl;
        stream << "NUM_PEAKS " << eic->peaks.size() << endl;
        for (unsigned int i = 0; i < eic->peaks.size(); i++) {
            Peak p = eic->peaks[i];
            stream << p.minpos << "\t" << p.pos << "\t" << p.maxpos << endl;
        }
        stream << "END_PEAKS" << endl;
        stream << "END_EIC" << endl;
    }

    stream << "END_SET" << endl;

    stream.flush();
    stream.close();
}

vector<EIC*> EIC::decode(string filePath, vector<mzSample*> samples) {
    vector<EIC*> eics;

    ifstream stream(filePath, ios::in);
    string line;

    bool isInCoords = false;
    bool isInPeaks = false;
    bool isInSampleName = false;
    int numEICs = 0;
    int numCoords = 0;
    int numPeaks = 0;
    int eicNum = 0;
    int coordNum = 0;
    int peakNum = 0;

    EIC *eic  = nullptr;

    while (getline(stream, line)) {
        // cout << line << endl;

        if (line.find("END_COORDS") == 0) {
            isInCoords = false;
            coordNum = 0;
        } else if (line.find("END_PEAKS") == 0) {
            isInPeaks = false;
            peakNum = 0;
        } else if (isInCoords) {
            vector<string> bits{};
            mzUtils::split(line, "\t", bits);

            eic->mz[coordNum] = stof(bits[0]);
            eic->rt[coordNum] = stof(bits[1]);
            eic->intensity[coordNum] = stof(bits[2]);
            eic->spline[coordNum] = stof(bits[3]);
            eic->baseline[coordNum] = stof(bits[4]);

            coordNum++;

        } else if (isInPeaks) {
            vector<string> bits{};
            mzUtils::split(line, "\t", bits);

            Peak p;
            p.minpos = stof(bits[0]);
            p.pos = stof(bits[1]);
            p.maxpos = stof(bits[2]);

            eic->peaks[peakNum] = p;

            peakNum++;

        } else if (line.find("NUM_EICS") == 0) {
            vector<string> bits{};
            mzUtils::split(line, " ", bits);
            numEICs = stoi(bits[1]);
            eics = vector<EIC*>(numEICs);
        } else if (line.find("NUM_COORDS") == 0) {
            vector<string> bits{};
            mzUtils::split(line, " ", bits);
            numCoords = stoi(bits[1]);

            eic->mz = vector<float>(numCoords);
            eic->rt = vector<float>(numCoords);
            eic->intensity = vector<float>(numCoords);
            eic->spline = vector<float>(numCoords);
            eic->baseline = vector<float>(numCoords);

            isInCoords = true;
        } else if (line.find("NUM_PEAKS") == 0) {
            vector<string> bits{};
            mzUtils::split(line, " ", bits);
            numPeaks = stoi(bits[1]);

            eic->peaks = vector<Peak>(numPeaks);

            isInPeaks = true;
        } else if (line.find("END_EIC") == 0) {
            eics[eicNum] = eic->clone();
            eicNum++;
        } else if (isInSampleName && eic) {
            isInSampleName = false;
            for (auto sample : samples) {
                if (sample->sampleName.find(line) == 0) {
                    eic->sample = sample;
                    break;
                }
            }
        } else if (line.find("START_EIC") == 0) {
            eic = new EIC();
            isInSampleName = true;
        }
    }

    return eics;
}
