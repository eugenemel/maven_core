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

	ifstream myfile(alignmentFile);
	if (!myfile.is_open()) { cerr << "Can't open file " << alignmentFile; return; }

	std::string line;
	int lineNum=0;

	AligmentSegment* lastSegment=0;

	while (getline(myfile,line) ) {
		lineNum++;
		vector<string>fields;
		mzUtils::split(line,'\t', fields);

		if (fields.size() >= 3 && lineNum > 1) {
			AligmentSegment* seg = new AligmentSegment(); 

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

float AligmentSegment::updateRt(float oldRt) { 
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

			AligmentSegment* seg = NULL;
			for( AligmentSegment* x: alignmentSegments[sampleName] ) { 
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

		cerr << "doSegmentedAligment: " << sampleName << "\tcorected=" << corcount << endl;
	}
}
