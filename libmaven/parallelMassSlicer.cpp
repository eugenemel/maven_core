#include "parallelMassSlicer.h"

void ParallelMassSlicer::algorithmA() {
    delete_all(slices);
    slices.clear();
    cache.clear();
    map< string, int> seen;

    for(unsigned int i=0; i < samples.size(); i++) {
        for(unsigned int j=0; j < samples[i]->scans.size(); j++ ) {
            Scan* scan = samples[i]->scans[j];
            if ( scan->filterLine.empty() ) continue;
            if ( seen.count( scan->filterLine ) ) continue;
            mzSlice* s = new mzSlice(scan->filterLine);
            slices.push_back(s);
            seen[ scan->filterLine ]=1;
        }
    }
    cout << "#algorithmA" << slices.size() << endl;
}

void ParallelMassSlicer::algorithmB(float userPPM, float minIntensity, int rtStep) { 
	delete_all(slices);
	slices.clear();
	cache.clear();

	float rtWindow=2.0;
	this->_precursorPPM=userPPM;
	this->_minIntensity=minIntensity;

	if (samples.size() > 0 and rtStep > 0 ) rtWindow = (samples[0]->getAverageFullScanTime()*rtStep);

    cout << "#Parallel algorithmB:" << endl;
    cout << " PPM=" << userPPM << endl;
    cout << " rtWindow=" << rtWindow << endl;
    cout << " rtStep=" << rtStep << endl;
    cout << " minCharge="     << _minCharge << endl;
    cout << " maxCharge="     << _maxCharge << endl;
    cout << " minIntensity="  << _minIntensity << endl;
    cout << " maxIntensity="  << _maxIntensity << endl;
    cout << " minMz="  << _minMz << endl;
    cout << " maxMz="  << _maxMz << endl;
    cout << " minRt="  << _minRt << endl;
    cout << " maxRt="  << _maxRt << endl;

#ifdef OMP_PARALLEL
     #pragma omp parallel for ordered num_threads(4) schedule(static)  
#endif
	for(unsigned int i=0; i < samples.size(); i++) {
		//if (slices.size() > _maxSlices) break;
		//
        //Scan* lastScan = NULL;
        cout << "#algorithmB:" << samples[i]->sampleName << endl;

		for(unsigned int j=0; j < samples[i]->scans.size(); j++ ) {
			Scan* scan = samples[i]->scans[j];
			if (scan->mslevel != 1 ) continue;
            if (_maxRt and !isBetweenInclusive(scan->rt,_minRt,_maxRt)) continue;
			float rt = scan->rt;

                vector<int> charges;
                if (_minCharge > 0 or _maxCharge > 0) charges = scan->assignCharges(userPPM);

                for(unsigned int k=0; k < scan->nobs(); k++ ){
                if (_maxMz and !isBetweenInclusive(scan->mz[k],_minMz,_maxMz)) continue;
                if (_maxIntensity and !isBetweenInclusive(scan->intensity[k],_minIntensity,_maxIntensity)) continue;
                if ((_minCharge or _maxCharge) and !isBetweenInclusive(charges[k],_minCharge,_maxCharge)) continue;

				float mz = scan->mz[k];
				float mzmax = mz + mz/1e6*_precursorPPM;
				float mzmin = mz - mz/1e6*_precursorPPM;

               // if(charges.size()) {
                    //cout << "Scan=" << scan->scannum << " mz=" << mz << " charge=" << charges[k] << endl;
               // }
				mzSlice* Z;
#ifdef OMP_PARALLEL
    #pragma omp critical
#endif
			    Z = sliceExists(mz,rt);
#ifdef OMP_PARALLEL
    #pragma omp end critical
#endif
		
				if (Z) {  //MERGE
                    //cout << "Merged Slice " <<  Z->mzmin << " " << Z->mzmax
					//<< " " << scan->intensity[k] << "  " << Z->ionCount << endl;

					Z->ionCount = std::max((float) Z->ionCount, (float ) scan->intensity[k]);
					Z->rtmax = std::max((float)Z->rtmax, rt+2*rtWindow);
					Z->rtmin = std::min((float)Z->rtmin, rt-2*rtWindow);
					Z->mzmax = std::max((float)Z->mzmax, mzmax);
					Z->mzmin = std::min((float)Z->mzmin, mzmin);
					//make sure that mz windown doesn't get out of control
					if (Z->mzmin < mz-(mz/1e6*userPPM)) Z->mzmin =  mz-(mz/1e6*userPPM);
					if (Z->mzmax > mz+(mz/1e6*userPPM)) Z->mzmax =  mz+(mz/1e6*userPPM);
					Z->mz =(Z->mzmin+Z->mzmax)/2; Z->rt=(Z->rtmin+Z->rtmax)/2;
                    //cout << Z->mz << " " << Z->mzmin << " " << Z->mzmax << " "
					//<< ppmDist((float)Z->mzmin,mz) << endl;
				} else { //NEW SLICE
					//				if ( lastScan->hasMz(mz, userPPM) ) {
                    //cout << "\t" << rt << "  " << mzmin << "  "  << mzmax << endl;
					mzSlice* s = new mzSlice(mzmin,mzmax, rt-2*rtWindow, rt+2*rtWindow);
					s->ionCount = scan->intensity[k];
					s->rt=scan->rt; 
					s->mz=mz;
#ifdef OMP_PARALLEL
    #pragma omp critical
#endif
					addSlice(s);
#ifdef OMP_PARALLEL
    #pragma omp end critical
#endif
				}

                //if ( slices.size() % 10000 == 0) cout << "ParallelMassSlicer count=" << slices.size() << endl;
			} //every scan m/z
			//lastScan = scan;
		} //every scan
	} //every samples

    cout << "#algorithmB:  Found=" << slices.size() << " slices" << endl;
	sort(slices.begin(),slices.end(), mzSlice::compIntensity);
}

void ParallelMassSlicer::addSlice(mzSlice* s) {
		slices.push_back(s);
		int mzRange = s->mz*10;
		cache.insert( pair<int,mzSlice*>(mzRange, s));
	}

void ParallelMassSlicer::algorithmC(float ppm, float minIntensity, float rtWindow, int topN=20, int minCharge=1) {
        delete_all(slices);
        slices.clear();
        cache.clear();

        for(unsigned int i=0; i < samples.size(); i++) {
            mzSample* s = samples[i];
            for(unsigned int j=0; j < s->scans.size(); j++) {
                Scan* scan = samples[i]->scans[j];

                if (scan->mslevel != 1 ) continue;
                vector<int>chargeState = scan->assignCharges(ppm);

                vector<int> positions = scan->intensityOrderDesc();
                for(unsigned int k=0; k< positions.size() && k<topN; k++ ) {
                    int pos = positions[k];
                    if(scan->intensity[pos] < minIntensity) continue;
                    if(scan->isMonoisotopicPrecursor(scan->mz[pos],ppm) == false) continue;
                    //if (chargeState[pos] < minCharge )  continue;

                    float rt = scan->rt;
                    float mz = scan->mz[ pos ];
                    float mzmax = mz + mz/1e6*ppm;
                    float mzmin = mz - mz/1e6*ppm;

                    mzSlice* slice = sliceExists(mz,rt);
                    if(!slice ) {
                        mzSlice* s = new mzSlice(mzmin,mzmax, rt-2*rtWindow, rt+2*rtWindow);
                        s->ionCount = scan->intensity[pos];
                        s->rt=scan->rt;
                        s->mz=mz;
                        slices.push_back(s);
                        int mzRange = mz*10;
                        cache.insert( pair<int,mzSlice*>(mzRange, s));
                    } else if ( slice->ionCount < scan->intensity[pos]) {
                            slice->ionCount = scan->intensity[pos];
                            slice->rt = scan->rt;
                            slice->mz = mz;
                    }

                }
            }
        }
        cout << "#algorithmC" << slices.size() << endl;
    }

void ParallelMassSlicer::algorithmD(float ppm, float rtWindow) {        //features that have ms2 events
        delete_all(slices);
        slices.clear();
        cache.clear();

        for(unsigned int i=0; i < samples.size(); i++) {
            mzSample* s = samples[i];
            for(unsigned int j=0; j < s->scans.size(); j++) {
                Scan* scan = samples[i]->scans[j];
                if (scan->mslevel != 2 ) continue;
                float rt = scan->rt;
                float mz = scan->precursorMz;
                float mzmax = mz + mz/1e6*ppm;
                float mzmin = mz - mz/1e6*ppm;

                if(! sliceExists(mz,rt) ) {
                    mzSlice* s = new mzSlice(mzmin,mzmax, rt-2*rtWindow, rt+2*rtWindow);
                    s->ionCount = scan->totalIntensity();
                    s->rt=scan->rt;
                    s->mz=mz;
                    slices.push_back(s);
                    int mzRange = mz*10;
                    cache.insert( pair<int,mzSlice*>(mzRange, s));
                }
            }
        }
        cout << "#algorithmD" << slices.size() << endl;
}

void ParallelMassSlicer::algorithmE(float ppm, float rtHalfWindowInMin) {        //features that have ms2 events
        delete_all(slices);
        slices.clear();
        cache.clear();

		vector<mzSlice*> sample_slices;

        for(unsigned int i=0; i < samples.size(); i++) {
            mzSample* s = samples[i];

            for(unsigned int j=0; j < s->scans.size(); j++) {
                Scan* scan = samples[i]->scans[j];
                if (scan->mslevel != 2 ) continue;
                float rt = scan->rt;
                float mz = scan->precursorMz;
                float mzmax = mz + mz/1e6f*ppm;
                float mzmin = mz - mz/1e6f*ppm;

                mzSlice* s = new mzSlice(mzmin,mzmax, rt-rtHalfWindowInMin, rt+rtHalfWindowInMin);
                s->rt=scan->rt;
                s->mz=mz;
                s->deleteFlag = false;
				sample_slices.push_back(s);
            }
		}

        unsigned long deleteCounter = mergeOverlappingSlices(sample_slices, ppm, true);

        cout << deleteCounter << " mz slices flagged for exclusion." << endl;

		for (mzSlice* x: sample_slices) { 
			if (!x->deleteFlag) slices.push_back(x); 
		}

        cout << "#algorithmE number of mz slices after merge: " << slices.size() << endl;
}

bool ParallelMassSlicer::isOverlapping(mzSlice *a, mzSlice *b){

      bool isMzOverlapping = checkOverlap(a->mzmin, a->mzmax, b->mzmin, b->mzmax) > 0.0f;
      bool isRtOverlapping = checkOverlap(a->rtmin, a->rtmax, b->rtmin, b->rtmax) > 0.0f;

    return isMzOverlapping && isRtOverlapping;
}

//Issue 584
unsigned long ParallelMassSlicer::mergeOverlappingSlices(vector<mzSlice*>& slices, float ppm, bool debug){
    unsigned long deleteCounter = 0;
    unsigned int iterationCounter = 0;

    bool isCheckMerges = true;

    while (isCheckMerges) {

        if (debug) {
            iterationCounter++;
            cout << "ParallelMassSlicer::mergeOverlappingSlices() iteration " << iterationCounter << ": " << deleteCounter << "/" << slices.size() << " slices marked for deletion." << endl;
        }

        isCheckMerges = false;

        //sort only by m/z, rely on RT tolerance downstream
        sort(slices.begin(), slices.end(), [ ](const mzSlice* lhs, const mzSlice* rhs){
            return lhs->mz < rhs->mz;
        });

        for(unsigned int i=0; i < slices.size(); i++ ) {

            mzSlice* a  = slices[i];

            if (a->deleteFlag) continue; //skip over if already marked

            for(unsigned int j=i+1; j < slices.size(); j++ ) {

                mzSlice* b  = slices[j];

                //Once the distance in m/z exceeds user-specified limit, no need to keep comparing for merges.
                //Note that b->mzmin < a->mzmax comparison is necessary
                if (b->mzmin > a->mzmax && ppmDist(a->mzmax, b->mzmin) > ppm) break;

                //skip over mz slices that have already been merged
                if (b->deleteFlag) continue;

                if (ParallelMassSlicer::isOverlapping(a, b)) {

                    //b swallows up a
                    b->rtmin = min(a->rtmin, b->rtmin);
                    b->rtmax = max(a->rtmax, b->rtmax);

                    b->mzmin = min(a->mzmin, b->mzmin);
                    b->mzmax = max(a->mzmax, b->mzmax);

                    b->mz  = b->mzmin + ((b->mzmax - b->mzmin)/2.0f);
                    b->rt  = b->rtmin + ((b->rtmax - b->rtmin)/2.0f);

                    //a is marked to be ignored in the future
                    a->deleteFlag = true;
                    deleteCounter++;

                    //Issue 584: slices need to be checked after each round of merging for additional merges to avoid transitivity issues
                    isCheckMerges = true;
                }

                //b now contains all of the information in a, so proceed to next index
                if (a->deleteFlag) break;
            }
        }

    }

    return deleteCounter;
}


mzSlice*  ParallelMassSlicer::sliceExists(float mz, float rt) {
	pair< multimap<int, mzSlice*>::iterator,  multimap<int, mzSlice*>::iterator > ppp;
	ppp = cache.equal_range( (int) (mz*10) );
	multimap<int, mzSlice*>::iterator it2 = ppp.first;

	float bestDist=FLT_MAX; 
	mzSlice* best=NULL;

	for( ;it2 != ppp.second; ++it2 ) {
		mzSlice* x = (*it2).second; 
		if (mz >= x->mzmin && mz <= x->mzmax && rt >= x->rtmin && rt <= x->rtmax) {
			float mzc = (x->mzmax - x->mzmin)/2;
			float rtc = (x->rtmax - x->rtmin)/2;
			float d = sqrt((POW2(mz-mzc) + POW2(rt-rtc)));
			if ( d < bestDist ) { best=x; bestDist=d; }
		}
	}
	return best;
}
