#include "mzSample.h"
#include "parallelMassSlicer.h"
#include <cmath>
#include "PolyAligner.h"


Aligner::Aligner() {
	maxItterations=10;
	polynomialDegree=3;
}

void Aligner::doAlignment(vector<PeakGroup*>& peakgroups) {
	if (peakgroups.size() == 0) return;

	//store groups into private variable
	allgroups = peakgroups;

    samples.clear();
	set<mzSample*> samplesSet;
	for (unsigned int i=0; i < peakgroups.size();  i++ ) {
			for ( unsigned int j=0; j < peakgroups[i]->peakCount(); j++ ) {
					Peak& p = peakgroups[i]->peaks[j];
					mzSample* sample = p.getSample();
					if (sample) samplesSet.insert(sample);
			}
	}

	//unique list of samples
	samples.resize(samplesSet.size());
	copy(samplesSet.begin(), samplesSet.end(),samples.begin());

    for(unsigned int i=0; i < samples.size(); i++ ) {
        samples[i]->saveOriginalRetentionTimes();
    }

	 saveFit();
	 double R2_before = checkFit();


     cerr << "Max Iterations: " << maxItterations << endl;
     for(int iter=0; iter < maxItterations; iter++) {

       PolyFit(polynomialDegree);
        double R2_after = checkFit();
        cerr << "Iteration:" << iter << " R2_before" << R2_before << " R2_after=" << R2_after << endl;

		if (R2_after > R2_before) {
            cerr << "done...restoring previous fit.." << endl;
			restoreFit();
			break;
		} else {
			saveFit();
		}
		R2_before = R2_after;
	 }
}


void Aligner::saveFit() {
	cerr << "saveFit()" << endl;
	fit.clear();
	fit.resize(samples.size());
	for(unsigned int i=0; i < samples.size(); i++ ) {
		fit[i].resize(samples[i]->scans.size());
		for(unsigned int ii=0; ii < samples[i]->scans.size(); ii++ ) {
			fit[i][ii]=samples[i]->scans[ii]->rt;
		}
	}
}

void Aligner::restoreFit() {
	cerr << "restoreFit() " << endl;
	for(unsigned int i=0; i < samples.size(); i++ ) {
		for(unsigned int ii=0; ii < samples[i]->scans.size(); ii++ ) {
			samples[i]->scans[ii]->rt = fit[i][ii];
		}
	}
}
vector<double> Aligner::groupMeanRt() {
		//find retention time deviation
		vector<double> groupRt(allgroups.size());
		for (unsigned int i=0; i < allgroups.size(); i++ ) groupRt[i]=allgroups[i]->medianRt();
		return(groupRt);
}

double Aligner::checkFit() { 
	vector<double> groupRt  = groupMeanRt();

	double sumR2=0;
	for(unsigned int i=0; i < allgroups.size(); i++ ) {
		for(unsigned int j=0; j < allgroups[i]->peakCount(); j++ ) {
			sumR2 += POW2(groupRt[i]-allgroups[i]->peaks[j].rt);
		}
	}
	cerr << "groups=" << allgroups.size() << " checkFit() " << sumR2 << endl;
	return sumR2;
}

void Aligner::PolyFit(int poly_align_degree) {

	if (allgroups.size() < 2 ) return;
	cerr << "Align: " << allgroups.size() << endl;

	vector<double> allGroupsMeansRt  = groupMeanRt();

	for (unsigned int s=0; s < samples.size(); s++ ) {
			mzSample* sample = samples[s];
			if (sample == NULL) continue;

			StatisticsVector<float>subj;
			StatisticsVector<float>ref;
			int n=0;

            map<int,int>duplicates;
			for(unsigned int j=0; j < allgroups.size(); j++ ) {
				Peak* p = allgroups[j]->getPeak(sample);
				if (!p) continue;
                if (!p || p->rt <= 0 || allGroupsMeansRt[j] <=0 ) continue;

                int intTime = (int) p->rt*100;
                duplicates[intTime]++;
                if ( duplicates[intTime] > 5 ) continue;

                ref.push_back(allGroupsMeansRt[j]);
				subj.push_back(p->rt);
				n++; 
			}
			if ( n < 10 ) continue;

			PolyAligner polyAligner(subj,ref);
            AlignmentStats* stats = polyAligner.optimalPolynomial(1,poly_align_degree,10);

           if (stats->transformImproved()) {

                bool failedTransformation=false;
                for(unsigned int ii=0; ii < sample->scans.size(); ii++ ) {
                    double newrt =  stats->predict(sample->scans[ii]->rt);
                    if (std::isnan(newrt) || std::isinf(newrt))  failedTransformation = true;
                    break;
                }

                if (!failedTransformation) {
                    for(unsigned int ii=0; ii < sample->scans.size(); ii++ ) {
                        sample->scans[ii]->rt = stats->predict(sample->scans[ii]->rt);
                    }

                    for(unsigned int ii=0; ii < allgroups.size(); ii++ ) {
                        Peak* p = allgroups[ii]->getPeak(sample);
                        if (p)  p->rt = stats->predict(p->rt);
                    }
                }
            } else 	{
                cerr << "APPLYTING TRANSFORM FAILED! " << endl;
            }
    }
}

void Aligner::Fit(int ideg) {

	if (allgroups.size() < 2 ) return;
	cerr << "Align: " << allgroups.size() << endl;

	double* x  = new double[allgroups.size()];
	double* ref =  new double[allgroups.size()];
    double* b =  new double[allgroups.size()];
	double* c =  new double[allgroups.size()];
	double* d =  new double[allgroups.size()];

//	ofstream temp;
//	temp.open("/tmp/fit.csv");

    //polynomial fit with maximum possible degree 5
    int maxdeg=5;
    if(ideg > maxdeg) ideg=maxdeg;
	double* result = new double[maxdeg];
	double* w = new double[maxdeg*maxdeg];

	//vector<double> groupRt  = groupMeanRt();

	for (unsigned int s=0; s < samples.size(); s++ ) {
			mzSample* sample = samples[s];
			if (sample == NULL) continue;
			map<int,int>duplicates;

			int n=0;
			StatisticsVector<float>diff;
			for(unsigned int j=0; j < allgroups.size(); j++ ) {
				Peak* p = allgroups[j]->getPeak(sample);
				if (!p) continue;
				if (p->rt <= 0) continue;
                //if (p->quality < 0.5 ) continue;
                int intTime = (int) p->rt*100;
                duplicates[intTime]++;
                if ( duplicates[intTime] > 5 ) continue;

				ref[n]=allgroups[j]->medianRt();
				x[n]=p->rt; 

                diff.push_back(POW2(x[n]-ref[n]));
				n++; 
			}
			if ( n == 0 ) continue;

            //double meanDiv = diff.mean();
            double stdDiv  = diff.stddev();
            if ( stdDiv == 0) continue;

            //REMOVE OUTLIERS
            int cutpos = n*0.95;
            double cut=diff[n-1]; if(cutpos>=0) cut=diff[cutpos];

			int removedCount=0;
			for(int ii=0; ii < n; ii++ ) {
                double deltaX = POW2(x[ii]-ref[ii]);
                if(deltaX > cut) {
                    //cerr << deltaX << " " << x[ii] << " " << ref[ii] << " " << meanDiv << " " << stdDiv << endl;
					x[ii]=0; 
					ref[ii]=0; 
					removedCount++;
				}
			}
			if (n - removedCount < 10) {
				cerr << "\t Can't align.. too few peaks n=" << n << " removed=" << removedCount << endl;
				continue;
			}

            //SORT AND ALIGN
			double R_before=0;
			for(int ii=0; ii < n; ii++)  R_before += POW2(ref[ii] - x[ii]);

			double R_after=0;   
            int transformedFailed=0;
            sort_xy(x, ref, n, 1, 0);
            leasqu(n, x, ref, ideg, w, maxdeg, result);	//polynomial fit
            //for(int ii=0; ii < ideg; ii++ ) { cerr << "Fit: " << ii << " " << result[ii] << endl; }


            //for(int ii=0; ii < n; ii++) d[ii]=c[ii]=b[ii]=0;
            //spline(n, x, ref, b, c, d);

			for(int ii=0; ii < n; ii++)  { 
                double newrt = leasev(result, ideg, x[ii]);
                //float newrt = seval(n, x[ii], x, ref, b, c, d);
				//temp << s << ", " << x[ii] << "," << ref[ii] << "," <<  newrt << endl;
                if (newrt != newrt || !std::isinf(newrt)) {
                    //cerr << "Polynomical transform failed! (A) " << x[ii] << "-->" << newrt <<  endl;
                    transformedFailed++;
				}  else {
					R_after  += POW2(ref[ii] - newrt);
				}
			}

            if(R_after > R_before ) {
                cerr << "Skipping alignment of " << sample->sampleName << " failed=" << transformedFailed << endl;
                 continue;
            }

            int failedTransformation=0;
            double zeroOffset =  leasev(result, ideg, 0);

            //cerr << "Alignment Sample: n=" << n << " mean=" << meanDiv << " std=" << stdDiv;
            //cerr << "\t improved R2" << R_before-R_after << " BAD=" << transformedFailed << endl;
            // cerr << "zeroOffset=" << zeroOffset << endl;
            for(unsigned int ii=0; ii < sample->scans.size(); ii++ ) {
                //float newrt = seval(n, sample->scans[ii]->rt, x, ref, b, c, d);
                double newrt =  leasev(result, ideg, sample->scans[ii]->rt)-zeroOffset;
                if (!std::isnan(newrt) && !std::isinf(newrt)) { //nan check
                    sample->scans[ii]->rt = newrt;
                } else {
                    cerr << "error: " << sample->scans[ii]->rt << " " << newrt << endl;
                    failedTransformation++;
                }
            }

            for(unsigned int ii=0; ii < allgroups.size(); ii++ ) {
                Peak* p = allgroups[ii]->getPeak(sample);
                if (p) {
                    //float newrt = seval(n, p->rt, x, ref, b, c, d);
                    double newrt = leasev(result, ideg, p->rt)-zeroOffset;
                    if (!std::isnan(newrt) && !std::isinf(newrt)) { //nan check
                        p->rt = newrt;
                    } else {
                        //cerr << "Polynomial transformed failed! (peak)" << p->rt << endl;
                    }
                }
            }

            if (failedTransformation) {
                cerr << "APPLYTING TRANSFORM FAILED: " << failedTransformation << endl;
            }

	}

	delete[] result;
	delete[] w;
	delete[] x;
	delete[] ref;
	//temp.close();
	delete[] b;
	delete[] c;
	delete[] d;
}

void Aligner::loadAlignmentFile(string alignmentFile) { 

	//aligment file format
	//sample  rt      rt_update
	//DOplasma-set1-b1-blank-inj1-C18TBA-neg  0.0     1.4722584874037974
	//DOplasma-set1-b1-blank-inj1-C18TBA-neg  0.024231014164164164    1.4910885552525142
    //
    //Note that this expects a block of all the same samples, then the next set of samples, etc

	ifstream myfile(alignmentFile);
	if (!myfile.is_open()) { cerr << "Can't open file " << alignmentFile; return; }

	std::string line;
	int lineNum=0;

    AlignmentSegment* lastSegment=0;

	while (getline(myfile,line) ) {
		lineNum++;
		vector<string>fields;
		mzUtils::split(line,'\t', fields);

		if (fields.size() >= 3 && lineNum > 1) {
            AlignmentSegment* seg = new AlignmentSegment();

			seg->sampleName   = fields[0];
			seg->seg_start = 0;
			seg->seg_end =  string2float(fields[1]);

			seg->new_start = 0;
			seg->new_end =  string2float(fields[2]);

			if (lastSegment and lastSegment->sampleName == seg->sampleName) { 
				seg->seg_start = lastSegment->seg_end;
				seg->new_start = lastSegment->new_end;
			}

			alignmentSegments[seg->sampleName].push_back(seg);
			lastSegment = seg;
		}
	}

	cerr << "Aligner::loadAlignmentFile() " << alignmentSegments.size() << "\t" << lineNum << endl;
}

float AlignmentSegment::updateRt(float oldRt) {
		//fractional distance from start of a segement
		if (oldRt >= seg_start and oldRt < seg_end) {
			float   frac = (oldRt-seg_start)/(seg_end - seg_start);
			return  new_start + frac*(new_end-new_start);
		} else {
			cerr << "Bad Map: " << oldRt << "\t" << seg_start << "\t" << seg_end << endl;
			return oldRt; // could not correct return old rt
		}
}

void Aligner::doSegmentedAligment() {

	cerr << "Aligner::doSegmentedAligment()" << "samples=" << samples.size() << endl;

	for (mzSample* sample: samples ) {
        if (sample == NULL) continue;
        sample->saveOriginalRetentionTimes();
	
		string sampleName = sample->sampleName;
        //mzUtils::replace(sampleName,".mzXML","");

		if ( alignmentSegments.count(sampleName)  == 0) { 
				cerr << "Can't find alignment information for sample " << sampleName << endl;
				continue;
		}

		int corcount=0;
		for(int ii=0; ii < sample->scans.size(); ii++ ) {
			Scan* scan = sample->scans[ii];

            AlignmentSegment* seg = nullptr;
            for( AlignmentSegment* x: alignmentSegments[sampleName] ) {
				if(scan->rt >= x->seg_start and scan->rt < x->seg_end) {
						seg=x; break;
				}
			}

			if(seg) {
				double newRt = seg->updateRt(scan->rt);
				scan->rt = newRt;
				corcount++;
			} else {
				cerr << "Can't find segment for: " << sampleName << "\t" << scan->rt << endl;
			}
		}

        cerr << "doSegmentedAligment: " << sampleName << "\tcorrected=" << corcount << endl;
	}
}

/**
 * MZ ALIGNER 2
 * Following this point are new developments, designed for use with new manual curation max EIC approach.
 */

bool AnchorPoint::setEICRtValue(mzSlice *slice, int eic_smoothingWindow, float minPeakIntensity){

    EIC *eic = sample->getEIC(slice->mzmin, slice->mzmax, slice->rtmin, slice->rtmax, 1);

    eic->getSingleGlobalMaxPeak(eic_smoothingWindow);

    if (!eic->peaks.empty() && eic->peaks[0].peakIntensity >= minPeakIntensity){
        this->rt = eic->peaks[0].rt;
        this->isRtFromEIC = true;
    } else {
        this->isRtFromEIC = false;
    }

    if (eic) delete(eic);

    return isRtFromEIC;
}

string AnchorPointSet::toString() {
    return "mz=["\
            + to_string(slice->mzmin)\
            + " - " + to_string(slice->mzmax)\
            + "], rt=["\
            + to_string(slice->rtmin)\
            + " - "\
            + to_string(slice->rtmax)\
            + "] VALID? " + (isValid ? "true" : "false");
}

/**
 * @brief AnchorPointSet::compute
 *
 * Samples must have sample IDs in injection order of machine, otherwise interpolation will
 * not make much sense
 *
 * @param eicSamples
 * @param allSamples
 * @param eic_smoothingWindow
 */
void AnchorPointSet::compute(const vector<mzSample*>& allSamples){

    //This flag is set in the constructor
    if (!isValid) return;

    vector<mzSample*> foundEICSamples;

    //Retrieve RTs for anchor points from samples
    for (auto &x : allSamples) {
        AnchorPoint *anchorPoint = new AnchorPoint(x);

        bool isComputeEIC = false;
        if (eicSamples.empty()) {

            //if no EIC samples are specified, try to extract an EIC from all samples.
            isComputeEIC = true;

        } else {

            // if some samples have been designated as EIC-containing samples,
            // only try to extract an EIC for these samples.
            auto it = find(eicSamples.begin(), eicSamples.end(), x);
            isComputeEIC = it != eicSamples.end();

        }

        if (isComputeEIC) {
            bool isFoundEIC = anchorPoint->setEICRtValue(slice, eic_smoothingWindow, minPeakIntensity);
            if (isFoundEIC) {
                foundEICSamples.push_back(x);
                sampleToPoints.insert(make_pair(x, anchorPoint));
            }
        }
    }

    if (foundEICSamples.size() < minNumObservedSamples) {
        isValid = false; //will not use if no signal could be extracted for any samples.
        if (slice) {
            cerr << "Did not find enough samples for anchor point: " << toString() << endl;
        }
        return;
    }

    sort(foundEICSamples.begin(), foundEICSamples.end(), [](const mzSample* lhs, const mzSample* rhs){
        return lhs->sampleId < rhs->sampleId;
    });

    //interpolate for all samples that do not have RT values from the EIC.
    for (auto &x : allSamples) {

        int sampleId = x->sampleId;
        float rt;

        bool addInterpolatedSample = false;

        if (sampleId <= foundEICSamples[0]->sampleId) {
            rt = sampleToPoints[foundEICSamples[0]]->rt;
            addInterpolatedSample = (sampleId != foundEICSamples[0]->sampleId);
        } else if (sampleId > foundEICSamples[foundEICSamples.size()-1]->sampleId){
            rt = sampleToPoints[foundEICSamples[foundEICSamples.size()-1]]->rt;
            addInterpolatedSample = true;
        } else {
            for (unsigned int i = 1; i < foundEICSamples.size(); i++){
                if (sampleId > foundEICSamples[i-1]->sampleId && sampleId < foundEICSamples[i]->sampleId) {
                    int firstId = foundEICSamples[i-1]->sampleId;
                    int secondId = foundEICSamples[i]->sampleId;

                    float firstRt = sampleToPoints[foundEICSamples[i-1]]->rt;
                    float secondRt = sampleToPoints[foundEICSamples[i]]->rt;

                    float frac = static_cast<float>(sampleId-firstId) / static_cast<float>(secondId-firstId);

                    float interpolatedRt = firstRt + (secondRt-firstRt) *frac;

                    rt = interpolatedRt;
                    addInterpolatedSample = true;
                    break;

                } else {
                    rt = sampleToPoints[foundEICSamples[i]]->rt;
                }
            }
        }

        if (addInterpolatedSample) {

            AnchorPoint *anchorPoint = new AnchorPoint(x);
            anchorPoint->setInterpolatedRtValue(rt);

            sampleToPoints.insert(make_pair(x, anchorPoint));
        }

    }
}

/**
 * @brief AnchorPointSet::setEICSamplesByFilter
 * @param stringFilter
 *
 * fill out eicSamples vector using samples that match the string filter criteria.
 */
void AnchorPointSet::setEICSamplesByFilter(const vector<mzSample*>& allSamples, string stringFilter){

    transform(stringFilter.begin(), stringFilter.end(), stringFilter.begin(), ::toupper);

    eicSamples.clear();

    for (auto &x : allSamples) {

        string sampleName = x->sampleName;
        transform(sampleName.begin(), sampleName.end(), sampleName.begin(), ::toupper);

        if (sampleName.find(stringFilter) != string::npos) {
            eicSamples.push_back(x);
        }
    }
}

/**
 * @brief groupsToAnchorPoints
 * @param samples
 * @param peakGroups
 * @param eic_smoothingWindow
 * @return a vector or AnchorPoints from a set of peak groups, which can be exported using Aligner::exportPeakGroups().
 */
vector<AnchorPointSet> Aligner::groupsToAnchorPoints(vector<mzSample*>& samples, vector<PeakGroup*>& peakGroups, int eic_smoothingWindow, float minPeakIntensity) {

    //extra position for last RT in file
    vector<AnchorPointSet> anchorPointSetVector(peakGroups.size()+1);

    for (unsigned int i = 0; i < peakGroups.size(); i++) {

        PeakGroup *group = peakGroups[i];
        AnchorPointSet anchorPointSet(*group);

        anchorPointSet.eic_smoothingWindow = eic_smoothingWindow;
        anchorPointSet.minPeakIntensity = minPeakIntensity;

        anchorPointSet.compute(samples);

        anchorPointSetVector[i] = anchorPointSet;
    }

    //last point in file
    AnchorPointSet lastAnchorPointSet = AnchorPointSet::lastRt(samples);
    anchorPointSetVector[peakGroups.size()] = lastAnchorPointSet;

    return anchorPointSetVector;
}


map<mzSample*, vector<pair<float, float>>> Aligner::anchorPointSetToUpdatedRtMap(vector<AnchorPointSet>& anchorPoints, mzSample* refSample){

    sort(anchorPoints.begin(), anchorPoints.end(), [refSample] (const AnchorPointSet& lhs, const AnchorPointSet& rhs){
        if (lhs.isValid && rhs.isValid) {
            return lhs.sampleToPoints.at(refSample)->rt < rhs.sampleToPoints.at(refSample)->rt;
        } else {
            return lhs.slice->rtmax < rhs.slice->rtmax;
        }
    });

    //                      observedRt  referenceRt
    //                          rt      rt_update
    map<mzSample*, vector<pair<float, float>>> sampleToUpdatedRts{};

    for (auto &pt : anchorPoints) {

        //invalid anchor point sets are skipped.
        if (!pt.isValid) continue;

        for (auto it = pt.sampleToPoints.begin(); it != pt.sampleToPoints.end(); ++it) {

            mzSample* sample = it->first;
            AnchorPoint* point = it->second;

            float observedRt = point->rt;
            float referenceRt = pt.sampleToPoints[refSample]->rt;

            if (sampleToUpdatedRts.find(sample) != sampleToUpdatedRts.end()) {

                //simplest check for monotonicity
                pair<float, float> lastPair = sampleToUpdatedRts[sample].at(sampleToUpdatedRts[sample].size()-1);
                float lastObserved = lastPair.first;
                float lastReference = lastPair.second;

                if (observedRt >= lastObserved && lastReference >= lastReference) {
                    sampleToUpdatedRts[sample].push_back(make_pair(observedRt, referenceRt));
                }

            } else {

                pair<float, float> rtPair = make_pair(observedRt, referenceRt);
                vector<pair<float, float>> rtInfo = vector<pair<float, float>>{};
                rtInfo.push_back(rtPair);

                sampleToUpdatedRts.insert(make_pair(sample, rtInfo));
            }

        }
    }

    return sampleToUpdatedRts;
}

/**
 * @brief exportAlignmentFile
 * @param anchorPoints
 * @return boolean flag indicating if export is successful.
 */
void Aligner::exportAlignmentFile(vector<AnchorPointSet>& anchorPoints, mzSample* refSample, string outputFile) {

    map<mzSample*, vector<pair<float, float>>> sampleToUpdatedRts = anchorPointSetToUpdatedRtMap(anchorPoints, refSample);

    ofstream outputStream;
    outputStream.open(outputFile);
    outputStream << "sample\trt\trt_update\n";

    for (auto &x : sampleToUpdatedRts) {

        string sampleName = x.first->sampleName;

        for (auto &pt : x.second) {

            outputStream
                    << sampleName << "\t"
                    << pt.first << "\t"
                    << pt.second << "\n";
        }
    }

    outputStream.close();
}

/**
 * @brief AnchorPointSet::lastRt
 * @param allSamples
 *
 * @return an AnchorPointSet that maps the last RT value collected from all samples together.
 * This highly imperfect, and so other AnchorPointSets should come as close as possible to the end of the
 * gradient.
 */
AnchorPointSet AnchorPointSet::lastRt(vector<mzSample*>& allSamples) {
    AnchorPointSet lastAnchorPoint;
    for (auto samp : allSamples) {

        AnchorPoint *anchorPoint = new AnchorPoint(samp);
        anchorPoint->setInterpolatedRtValue(samp->scans.back()->rt);

        lastAnchorPoint.sampleToPoints.insert(make_pair(samp, anchorPoint));
    }

    return lastAnchorPoint;
}
