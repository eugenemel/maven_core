#include "mzSample.h"
#include "directinfusionprocessor.h"

PeakGroup::PeakGroup()  { 
    groupId=0;
    metaGroupId=0;
    groupRank=1000;

    maxIntensity=0;
    meanRt=0;
    meanMz=0;

    blankMax=0;
    blankSampleCount=0;
    blankMean=0;

    sampleMax=0;
    sampleCount=0;
    sampleMean=0;

    ms2EventCount=0;
    maxNoNoiseObs=0;
    maxPeakFracionalArea=0;
    maxSignalBaseRatio=0;
    maxSignalBaselineRatio=0;
    maxPeakOverlap=0;
    maxQuality=0;

    expectedRtDiff=-1;
    expectedAbundance=0;
    isotopeC13count=0;

    minRt=0;
    maxRt=0;

    minMz=0;
    maxMz=0;

    parent = nullptr;
    adduct = nullptr;
    compound = nullptr;
    deletedFlag=false;

    isFocused=false;
    displayName="";
    importedCompoundName="";
    compoundId="";
    compoundDb="";

    goodPeakCount=0;
    _type = None;

    chargeState=0;
    isotopicIndex=0;

    changePValue=0;
    changeFoldRatio=0;
    //children.reserve(0);
    peaks.reserve(0);

    srmPrecursorMz = 0.0f;
    srmProductMz = 0.0f;

}      

void PeakGroup::copyObj(const PeakGroup& o)  { 
    groupId= o.groupId;
    metaGroupId= o.metaGroupId;
    groupRank= o.groupRank;

    maxIntensity= o.maxIntensity;
    meanRt=o.meanRt;
    meanMz=o.meanMz;

    blankMax=o.blankMax;
    blankSampleCount=o.blankSampleCount;
    blankMean=o.blankMean;

    sampleMax=o.sampleMax;
    sampleCount=o.sampleCount;
    sampleMean=o.sampleMean;

    ms2EventCount=o.ms2EventCount;
    maxNoNoiseObs=o.maxNoNoiseObs;
    maxPeakFracionalArea=o.maxPeakFracionalArea;
    maxSignalBaseRatio=o.maxSignalBaseRatio;
    maxSignalBaselineRatio=o.maxSignalBaselineRatio;
    maxPeakOverlap=o.maxPeakOverlap;
    maxQuality=o.maxQuality;
    expectedRtDiff=o.expectedRtDiff;
    expectedAbundance = o.expectedAbundance;
    isotopeC13count=o.isotopeC13count;

    minRt=o.minRt;
    maxRt=o.maxRt;

    minMz=o.minMz;
    maxMz=o.maxMz;

    parent = o.parent;
    compound = o.compound;
    adduct = o.adduct;

    compoundId = o.compoundId;
    compoundDb = o.compoundDb;

    srmId=o.srmId;
    isFocused=o.isFocused;
    displayName=o.displayName;
    labels = o.labels;
    importedCompoundName=o.importedCompoundName;

    goodPeakCount=o.goodPeakCount;
    _type = o._type;
    tagString = o.tagString;
    searchTableName = o.searchTableName;
    deletedFlag = o.deletedFlag;


    fragMatchScore = o.fragMatchScore;
    fragmentationPattern = o.fragmentationPattern;

    changeFoldRatio = o.changeFoldRatio;
    changePValue    = o.changePValue;
    peaks = o.peaks;


    chargeState=o.chargeState;
    isotopicIndex=o.isotopicIndex;

    srmPrecursorMz = o.srmPrecursorMz;
    srmProductMz = o.srmProductMz;

    copyChildren(o);
}

PeakGroup::~PeakGroup() { clear(); }

void PeakGroup::copyChildren(const PeakGroup& o) {
    children = o.children;
    for(unsigned int i=0; i < children.size(); i++ ) children[i].parent = this;
}

bool PeakGroup::isPrimaryGroup() { 
    if(compound && compound->getPeakGroup() == this) return true;
    return false;
}

void PeakGroup::clear() { 
    deletePeaks();
    deleteChildren();
    meanMz  = 0;
    groupRank=INT_MAX;

    //Issue 235: try forcing these connections to nullptr
    parent = nullptr;
    compound = nullptr;
    adduct = nullptr;
}

Peak* PeakGroup::getSamplePeak(mzSample* sample) { 
    for (unsigned int i=0; i< peaks.size(); i++ ) {
        if (peaks[i].getSample() == sample )  return &peaks[i];
    }
    return nullptr;
}

void PeakGroup::deletePeaks() { 
    peaks.clear();
}

bool PeakGroup::deletePeak(unsigned int index) { 
    if ( index < children.size() ) {
        peaks.erase(peaks.begin()+index);
        return true;
    }
    return false;
}

float PeakGroup::meanRtW() { 
    if (peakCount() == 0) return 0;

    float mean=0; float Wtotal=0;
    for(unsigned int i=0; i < peakCount(); i++ ) Wtotal +=peaks[i].peakIntensity;

    if (Wtotal > 0 ) {
        for(unsigned int i=0; i < peakCount(); i++ ) mean +=  peaks[i].peakIntensity/Wtotal * peaks[i].rt;
        return mean;
    } else {
        for(unsigned int i=0; i < peakCount(); i++ ) mean += peaks[i].rt;
        return mean / peakCount();
    }
}

float PeakGroup::medianRt() { 
    vector<float> rts(peaks.size(),0);
    for(unsigned int i=0; i < peakCount(); i++ ) rts[i]=peaks[i].rt;
    return mzUtils::median(rts);
}

void PeakGroup::deleteChildren() { 
    children.clear();
}

bool PeakGroup::deleteChild(unsigned int index) { 
    if ( index < children.size() ) {
        children.erase(children.begin()+index);
        return true;
    }
    return false;
}

bool PeakGroup::deleteChild(PeakGroup* child ) {
    if (!child) return false;

    vector<PeakGroup>::iterator it;
    it = find(children.begin(),children.end(),child);
    if ( *it == child ) {
        cerr << "deleteChild: setting child to empty";
        child->clear();
        return true;
        //sort(children.begin(), children.end(),PeakGroup::compIntensity);
        //for(int i=0; i < children.size(); i++ ) { cerr << &children[i] << endl; }
    }

    return false;
}

//return intensity vectory ordered by samples 
vector<float> PeakGroup::getOrderedIntensityVector(vector<mzSample*>& samples, QType type) {

    if (samples.size() == 0) { vector<float>x; return x; } //empty vector;

    map<mzSample*,float> sampleOrder;
    vector<float>maxIntensity(samples.size(),0);

    for( unsigned int j=0; j < samples.size(); j++) {
        sampleOrder[samples[j]]=j;
        maxIntensity[j]=0;
    }

    for( unsigned int j=0; j < peaks.size(); j++) {
        Peak& peak = peaks.at(j);
        mzSample* sample = peak.getSample();

        if ( sampleOrder.count(sample) > 0 ) {
            int s  = sampleOrder[ sample ];
            float y = 0;
            switch (type)  {
            case AreaTop: y = peak.peakAreaTop; break;
            case Area: y = peak.peakAreaCorrected; break;
            case Height: y = peak.peakIntensity; break;
            case AreaNotCorrected: y = peak.peakArea; break;
            case RetentionTime: y = peak.rt; break;
            case Quality: y = peak.quality; break;
            case SNRatio: y = peak.signalBaselineRatio; break;
            case MS2Count: y = peak.ms2EventCount; break;
            default: y = peak.peakAreaTop; break;
            }

            //normalize
            if(sample) y *= sample->getNormalizationConstant();
            if(maxIntensity[s] < y) { maxIntensity[s] = y; }
        }
    }
    return maxIntensity;
}

void PeakGroup::computeAvgBlankArea(const vector<EIC*>& eics) { 

    if (peaks.size() == 0 ) return;

    //find range to fill in
    float rtmin = peaks[0].rtmin;
    float rtmax = peaks[0].rtmax;

    for (unsigned int i=1; i < peaks.size(); i++ ) {
        if (peaks[i].rtmin < rtmin) rtmin = peaks[i].rtmin;
        if (peaks[i].rtmax > rtmax) rtmax = peaks[i].rtmax;
    }
    rtmin = rtmin-0.25;
    rtmax = rtmax+0.25;

    float sum=0; int len=0;
    for(unsigned int i=0; i < eics.size(); i++ ) {
        EIC* eic = eics[i];
        if(eic->sample != nullptr && eic->sample->isBlank == false) continue;
        for(unsigned int pos=0; pos < eic->intensity.size(); pos++ ) {
            if ( eic->rt[pos] >= rtmin && eic->rt[pos] <= rtmax
                 && eic->intensity[pos] > 0) {
                sum += eic->intensity[pos];
                len++;
            }
        }
    }
    this->blankMean = 0; 	//default zero
    if ( len > 0 ) this->blankMean = (float) sum / len;
}

void PeakGroup::fillInPeaks(const vector<EIC*>& eics) {

    if (peaks.size() == eics.size()) return;
    if (peaks.size() == 0 ) return;

    //find range to fill in
    float rtmin = peaks[0].rtmin;
    float rtmax = peaks[0].rtmax;

    for (unsigned int i=1; i < peaks.size(); i++ ) {
        if (peaks[i].rtmin < rtmin) rtmin = peaks[i].rtmin;
        if (peaks[i].rtmax > rtmax) rtmax = peaks[i].rtmax;
    }

    int filledInCount=0;

    for(unsigned int i=0; i < eics.size(); i++ ) {
        EIC* eic = eics[i];
        if (eic == nullptr ) continue;
        if (eic->spline.size() ) continue;
        if (eic->intensity.size() == 0) continue;

        bool missing=true;

        for(unsigned int j=0; j < peaks.size(); j++ ) {
            if ( peaks[j].getEIC() == eic) {
                missing = false;
                break;
            }
        }

        if (missing) { //fill in peak
            int maxpos = 0;
            for(unsigned int pos=1; pos < eic->intensity.size()-1; pos++ ) {
                if ( eic != nullptr && eic->intensity[pos] != 0 && eic->mz[pos] != 0 &&
                     eic->rt[pos] >= rtmin && eic->rt[pos] <= rtmax
                     && eic->spline[pos] > eic->spline[pos-1] && eic->spline[pos] > eic->spline[pos+1]
                     ) {
                    if (maxpos != 0 && eic->intensity[pos] > eic->intensity[maxpos]) {
                        maxpos=pos;
                    } else {
                        maxpos=pos;
                    }
                }
            }

            if (maxpos != 0 && eic->intensity[maxpos] != 0 ) {
                Peak peak(eic,maxpos);
                eic->findPeakBounds(peak);
                eic->getPeakDetails(peak);
                this->addPeak(peak);
                filledInCount++;
            }
        }
    }

    //cerr << "fillInPeaks" << rtmin << " " << rtmax << " " << eics.size() << " " peaks.size() << endl;
    //    if (filledInCount > 0) { this->fillInPeaks(eics); }
}

void PeakGroup::reduce() { // make sure there is only one peak per sample

    map<mzSample*, Peak> maxPeaks;
    map<mzSample*, Peak>::iterator itr;

    if (peaks.size() < 2 ) return;

       //CURRENTLY NOT USED
//    float groupMeanRt=0;
//    float totalWeight=1;

//    for( unsigned int i=0; i < peaks.size(); i++)  { totalWeight +=  peaks[i].peakIntensity; }
//    for( unsigned int i=0; i < peaks.size(); i++)  { groupMeanRt += peaks[i].rt * peaks[i].peakIntensity/totalWeight;  }


    for( unsigned int i=0; i < peaks.size(); i++) {
        mzSample* c = peaks[i].getSample();
        //float rtdiff = abs(groupMeanRt-peaks[i].rt);

        /*
        //In each group, take peak that closest to the mean retention time of a group
        if ( maxPeaks.count(c) == 0 ||  rtdiff < abs( groupMeanRt - maxPeaks[c].rt) ) {
            maxPeaks[c].copyObj(peaks[i]);
        }
        */

        //new approach
//        //In each group, take the most intense peak
//        if (maxPeaks.find(c) == maxPeaks.end()) {
//            maxPeaks.insert(make_pair(c, peaks[i]));
//        } else if (maxPeaks[c].peakIntensity < peaks[i].peakIntensity){
//            maxPeaks[c] = peaks[i];
//        }

        //old approach
        if ( maxPeaks.find(c) == maxPeaks.end() || maxPeaks[c].peakIntensity < peaks[i].peakIntensity) {
            maxPeaks[c].copyObj(peaks[i]);
        }
    }

    peaks.clear();
    for( itr = maxPeaks.begin(); itr != maxPeaks.end(); ++itr){
        const Peak& peak = (*itr).second;
        addPeak(peak);
    }
    //	cerr << "\t\t\treduce() from " << startSize << " to " << peaks.size() << endl;
}

void PeakGroup::updateQuality() {
    maxQuality=0;
    goodPeakCount=0;
    for(unsigned int i=0; i< peaks.size(); i++) {
        if(peaks[i].quality > maxQuality) maxQuality = peaks[i].quality;
        if(peaks[i].quality > 0.5) goodPeakCount++;
    }
}

void PeakGroup::groupStatistics() { 
    float rtSum = 0;
    float mzSum = 0;
    maxIntensity = 0;

    blankMax =0;
    blankSampleCount=0;

    sampleMax=0;
    sampleCount=0;
    sampleMean=0;

    maxNoNoiseObs=0;
    minRt = 0;
    maxRt = 0;
    minMz = 0;
    maxMz = 0;

    maxPeakFracionalArea=0;
    maxQuality=0;
    goodPeakCount=0;
    maxSignalBaselineRatio=0;
    int nonZeroCount=0;

    for(unsigned int i=0; i< peaks.size(); i++) {
        if(peaks[i].pos != 0) { rtSum += peaks[i].rt; mzSum += peaks[i].baseMz; nonZeroCount++; }

        if(peaks[i].peakIntensity>maxIntensity) {
            maxIntensity = peaks[i].peakIntensity;
            meanMz=peaks[i].baseMz;
            meanRt=peaks[i].rt;
        }

        if(peaks[i].noNoiseObs>maxNoNoiseObs) maxNoNoiseObs = peaks[i].noNoiseObs;
        if(minRt == 0 || peaks[i].rtmin < minRt) minRt = peaks[i].rtmin;
        if(maxRt == 0 || peaks[i].rtmax > maxRt) maxRt = peaks[i].rtmax;
        if(minMz == 0 || peaks[i].mzmin < minMz) minMz = peaks[i].mzmin;
        if(maxMz == 0 || peaks[i].mzmax > maxMz) maxMz = peaks[i].mzmax;
        if(peaks[i].peakAreaFractional > maxPeakFracionalArea) maxPeakFracionalArea=peaks[i].peakAreaFractional;
        if(peaks[i].quality > maxQuality) maxQuality = peaks[i].quality;
        if(peaks[i].quality > 0.5) goodPeakCount++;
        if(peaks[i].signalBaselineRatio > maxSignalBaselineRatio) maxSignalBaselineRatio =  peaks[i].signalBaselineRatio;


        if(peaks[i].fromBlankSample) {
            blankSampleCount++;
            if(peaks[i].peakIntensity > blankMax) blankMax = peaks[i].peakIntensity;
        } else {
            sampleMean += peaks[i].peakIntensity;
            sampleCount++;
            if(peaks[i].peakIntensity > sampleMax) sampleMax = peaks[i].peakIntensity;
        }
    }

    if (sampleCount>0) sampleMean = sampleMean/sampleCount;

    if ( nonZeroCount ) {
        meanRt = rtSum/nonZeroCount;
        meanMz = mzSum/nonZeroCount;
    }

    groupOverlapMatrix();
}

void PeakGroup::groupOverlapMatrix() {

    for(unsigned int i=0; i< peaks.size(); i++) peaks[i].groupOverlapFrac=0;

    for(unsigned int i=0; i< peaks.size(); i++) {
        Peak& a = peaks[i];
        for(unsigned int j=i; j< peaks.size(); j++) {
            Peak& b = peaks[j];
            float overlap = checkOverlap(a.rtmin,a.rtmax,b.rtmin,b.rtmax); //check for overlap
            if (overlap > 0 ) { b.groupOverlapFrac += log(overlap); a.groupOverlapFrac += log(overlap); }

            /*
                    if ( overlap > 0.1 ) { 
						b.peakAreaFractional < 1 ? a.groupOverlapFrac += log(1-b.peakAreaFractional) : a.groupOverlapFrac += log(0.01);
						a.peakAreaFractional < 1 ? b.groupOverlapFrac += log(1-a.peakAreaFractional) : b.groupOverlapFrac += log(0.01);
					}
					*/
        }
    }
    //normalize
    for(unsigned int i=0; i< peaks.size(); i++)  peaks[i].groupOverlapFrac /= peaks.size();
}

void PeakGroup::summary() {
    cerr 	<< tagString << endl;
    cerr
            <<"\t" << "meanRt=" << meanRt << endl
            <<"\t" << "meanMz=" << meanMz << endl
            <<"\t" << "goodPeakCount=" << goodPeakCount << endl
            <<"\t" << "maxQuality=" <<  maxQuality << endl
            <<"\t" << "maxNoNoiseObs=" << maxNoNoiseObs << endl
            <<"\t" << "sampleCount=" << sampleCount << endl
            <<"\t" << "maxSignalBaselineRatio=" << maxSignalBaselineRatio << endl
            <<"\t" << "maxPeakFracionalArea=" << maxPeakFracionalArea << endl
            <<"\t" << "blankMean=" << blankMean << endl
            <<"\t" << "sampleMean=" << sampleMean << endl
            <<"\t" << "maxIntensity=" << maxIntensity << endl
            << endl;

    for (unsigned int i=0; i < peaks.size(); i++ ) {
        cerr << "\t\t" << "Q:" << peaks[i].quality<< " "
                << "pAf:" << peaks[i].peakAreaFractional<< " "
                << "noNf" << peaks[i].noNoiseFraction << " "
                << "noObs:" << peaks[i].noNoiseObs   << " "
                << "w:"<< peaks[i].width<< " "
                << "sn:" << peaks[i].signalBaselineRatio << " "
                << "ovp:" << peaks[i].groupOverlapFrac << endl;
    }

    for(unsigned int i=0; i < children.size(); i++ ) children[i].summary();
}

PeakGroup::PeakGroup(const PeakGroup& o)  { 
    copyObj(o);
}

PeakGroup& PeakGroup::operator=(const PeakGroup& o)  {
    copyObj(o);
    return *this;
}

bool PeakGroup::operator==(const PeakGroup& o)  {
    if ( this->meanMz == o.meanMz
         && this->meanRt == o.meanRt
         && this->minMz  == o.minMz
         && this->maxMz  == o.maxMz
         && this->minRt  == o.minRt
         && this->maxRt  == o.maxRt
    ) return true;
    return false;
}

bool PeakGroup::operator==(const PeakGroup* o)  {
    if ( this == o )  return true;
    return false;
}

Peak* PeakGroup::getPeak(mzSample* s ) {
    if ( s == nullptr ) return nullptr;
    for(unsigned int i=0; i < peaks.size(); i++ ) {
        if ( peaks[i].getSample() == s ) {
            return &peaks[i];
        }
    }
    return nullptr;
}


void PeakGroup::reorderSamples() {  
    std::sort(peaks.begin(), peaks.end(), Peak::compIntensity);
    for(unsigned int i=0; i < peaks.size(); i++ ) {
        mzSample* s = peaks[i].getSample();
        if ( s != nullptr ) s->setSampleOrder(i);
    }
}

string PeakGroup::getName() {
    string tag;
    if (compound) tag = compound->name;
    if (tagString.empty()) tag += " | " + tagString;
    if (srmId.empty()) tag +=  " | " + srmId;
    if (tag.empty()) tag = integer2string(groupId);
    return tag;
}

vector<Scan*> PeakGroup::getRepresentativeFullScans() {
    vector<Scan*>matchedscans;
    for(unsigned int i=0; i < peaks.size(); i++ ) {
        mzSample* sample = peaks[i].getSample();
        if ( sample == nullptr ) continue;
        Scan* scan = sample->getScan(peaks[i].scan);
        if (scan and scan->mslevel == 1) matchedscans.push_back(scan);
    }
    return matchedscans;
}

vector<Scan*> PeakGroup::getFragmentationEvents() {

    vector<Scan*> matchedscans;

    for(unsigned int i=0; i < peaks.size(); i++ ) {
        mzSample* sample = peaks[i].getSample();
        if (!sample) continue;

        for( unsigned int j=0; j < sample->scans.size(); j++ ) {
            Scan* scan = sample->scans[j];
            if (scan->mslevel <= 1) continue; //ms2 + scans only
            if (scan->rt < peaks[i].rtmin) continue;
            if (scan->rt > peaks[i].rtmax) break;
            if( scan->precursorMz >= minMz and scan->precursorMz <= maxMz) {
                matchedscans.push_back(scan);
            }
        }
    }
    return matchedscans;
}

void PeakGroup::findHighestPurityMS2Pattern(float prePpmTolr) {


    //build consensus ms2 specta
    vector<Scan*>ms2events = getFragmentationEvents();
    if (ms2events.size() < 1 ) return;

    Scan* best=ms2events.front();
    float highestPurity = best->getPrecursorPurity(prePpmTolr) * log10(best->nobs()+1);

    cerr << "findHighesPurityMS2Pattern: " << prePpmTolr << " ms2ev:" << ms2events.size() <<  " hP=" << highestPurity << endl;

    for( Scan* x: ms2events) {
        float xpurity = x->getPrecursorPurity(prePpmTolr) * log10(x->nobs()+1);
        if (xpurity > highestPurity) {
            highestPurity = xpurity;
            best = x;
        }
    }

    if (best) {
        fragmentationPattern = Fragment(best,0.01,1,1024);
        //for(Scan* s : ms2events) {  fragmentationPattern.addFragment(new Fragment(s,0,0.01,1024)); }
        //fragmentationPattern.consensus = nullptr;
        fragmentationPattern.sortByMz();;
        ms2EventCount = ms2events.size();
    }

}


void PeakGroup::computeFragPattern(float productPpmTolr)  {
    //build consensus ms2 specta
    vector<Scan*>ms2events = getFragmentationEvents();
    if (ms2events.size() == 0 ) return;
    sort(ms2events.begin(),ms2events.end(),Scan::compIntensity);

    Fragment f(ms2events[0],0.01,1,1024);
    for(Scan* s : ms2events) {  f.addFragment(new Fragment(s,0,0.01,1024)); }
    f.buildConsensus(productPpmTolr);
    f.consensus->sortByMz();
    fragmentationPattern = f.consensus;
    ms2EventCount = static_cast<int>(ms2events.size());
}

void PeakGroup::computePeaksSearchFragPattern(shared_ptr<PeaksSearchParameters> params) {

    //build consensus ms2 specta
//    vector<Scan*>ms2events = getFragmentationEvents();
//    if (ms2events.size() == 0 ) return;
//    sort(ms2events.begin(),ms2events.end(),Scan::compIntensity);

    //TODO: modify all of this
    //Fragment f(ms2events[0],0.01,1,1024);
//    for(Scan* s : ms2events) {  f.addFragment(new Fragment(s,0,0.01,1024)); }
//    f.buildConsensus(params->ms2PpmTolr);

//    f.consensus->sortByMz();
//    fragmentationPattern = f.consensus;

    vector<Scan*> ms2Scans = getFragmentationEvents();
    ms2EventCount = static_cast<int>(ms2Scans.size());

    if (ms2Scans.empty()) return;   //this can happen when a peakgroup is associated with the wrong meanMz.

    Fragment *f = nullptr;
    for (auto& scan : ms2Scans) {
        if (!f){
            f = new Fragment(scan,
                             params->scanFilterMinFracIntensity,
                             params->scanFilterMinSNRatio,
                             params->scanFilterMaxNumberOfFragments,
                             params->scanFilterBaseLinePercentile,
                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                             params->scanFilterPrecursorPurityPpm,
                             params->scanFilterMinIntensity);
        } else {
            Fragment *brother = new Fragment(scan,
                                             params->scanFilterMinFracIntensity,
                                             params->scanFilterMinSNRatio,
                                             params->scanFilterMaxNumberOfFragments,
                                             params->scanFilterBaseLinePercentile,
                                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                             params->scanFilterPrecursorPurityPpm,
                                             params->scanFilterMinIntensity);

            f->addFragment(brother);
        }
    }

    f->buildConsensus(params->consensusPpmTolr,
                      params->consensusIntensityAgglomerationType,
                      params->consensusIsIntensityAvgByObserved,
                      params->consensusIsNormalizeTo10K,
                      params->consensusMinNumMs2Scans,
                      params->consensusMinFractionMs2Scans
                      );

    f->consensus->sortByMz();

    this->fragmentationPattern = f->consensus;

    //cleanup
    if (f) delete(f);
}

void PeakGroup::computeDIFragPattern(shared_ptr<DirectInfusionSearchParameters> params){

    vector<Scan*> ms2Scans;
    for(unsigned int i=0; i < peaks.size(); i++ ) {
        mzSample* sample = peaks[i].getSample();
        if (!sample) continue;

        for (Scan *scan : sample->scans){
            if (scan->mslevel == 2 ){
                if (meanMz >= scan->getPrecMzMin() && meanMz <= scan->getPrecMzMax()){
                    ms2Scans.push_back(scan);
                }
            }
        }
    }

    if (ms2Scans.empty()) return;   //this can happen when a peakgroup is associated with the wrong meanMz.

    Fragment *f = nullptr;
    for (auto& scan : ms2Scans) {
        if (!f){
            f = new Fragment(scan,
                             params->scanFilterMinFracIntensity,
                             params->scanFilterMinSNRatio,
                             params->scanFilterMaxNumberOfFragments,
                             params->scanFilterBaseLinePercentile,
                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                             params->scanFilterPrecursorPurityPpm,
                             params->scanFilterMinIntensity);
        } else {
            Fragment *brother = new Fragment(scan,
                                             params->scanFilterMinFracIntensity,
                                             params->scanFilterMinSNRatio,
                                             params->scanFilterMaxNumberOfFragments,
                                             params->scanFilterBaseLinePercentile,
                                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                             params->scanFilterPrecursorPurityPpm,
                                             params->scanFilterMinIntensity);

            f->addFragment(brother);
        }
    }

    f->buildConsensus(params->consensusPpmTolr,
                      params->consensusIntensityAgglomerationType,
                      params->consensusIsIntensityAvgByObserved,
                      params->consensusIsNormalizeTo10K,
                      params->consensusMinNumMs2Scans,
                      params->consensusMinFractionMs2Scans
                      );

    f->consensus->sortByMz();

    this->fragmentationPattern = f->consensus;

    //cleanup
    if (f) delete(f);

}

Scan* PeakGroup::getAverageFragmentationScan(float productPpmTolr)  {
    //build consensus ms2 specta
    computeFragPattern(productPpmTolr);
    Scan* avgScan = new Scan(nullptr,0,0,0,0,0);

    for(unsigned int i=0; i<fragmentationPattern.mzs.size();i++) {
        avgScan->mz.push_back(fragmentationPattern.mzs[i]);
        avgScan->intensity.push_back(fragmentationPattern.intensity_array[i]);
    }
    return avgScan;
}

/**
 * Clustering function that organizes peak groups into clusters based on several independent metrics.
 *
 * EIC comparisons are always done on the highest-intensity sample.
 *
 * TODO: PeakGroup A has most intense sample as X, Should PeakGroup B use the same sample, or just some other (more intense) sample?
 * Should this function try to fall back to the most intense sample found in both A and B?
 *
 * TODO: extend this function to only consider peak groups with known m/z deltas from each other
 * for example, isotopic peaks, or differences between common adducts, etc
 *
 * TODO: more information about sample correlations
 *
 * TODO: option to do smoothing on correlations?
 *
 * @brief PeakGroup::clusterGroups
 * @param allgroups
 * @param samples
 * @param maxRtDiff
 * @param minSampleCorrelation
 * @param minPeakShapeCorrelation
 * @param ppm
 * @param mzDeltas
 */
void PeakGroup::clusterGroups(vector<PeakGroup> &allgroups, vector<mzSample*>samples, double maxRtDiff, double minSampleCorrelation, double minPeakShapeCorrelation, double ppm, vector<double> mzDeltas) {
    sort(allgroups.begin(),allgroups.end(), PeakGroup::compRt);
    int metaGroupId = 0;

    //clear cluster information
    for(unsigned int i=0; i<allgroups.size(); i++) allgroups[i].metaGroupId=0;
    map<int,PeakGroup*>parentGroups;

    for(unsigned int i=0; i<allgroups.size(); i++) {
        PeakGroup& grp1 = allgroups[i];

        if (grp1.metaGroupId == 0) {  //create new cluster
            grp1.metaGroupId=++metaGroupId;
            parentGroups[metaGroupId]=&grp1;
        }

        //cluster parent
        PeakGroup* parent = parentGroups[ metaGroupId ];

        mzSample* largestSample=nullptr;
        double maxIntensity=0;

        for(unsigned int i=0; i < grp1.peakCount(); i++ ) {
            mzSample* sample = grp1.peaks[i].getSample();
            if ( grp1.peaks[i].peakIntensity > maxIntensity ) largestSample=sample;
        }

        if (largestSample == nullptr ) continue;
        vector<float>peakIntensityA = grp1.getOrderedIntensityVector(samples,PeakGroup::AreaTop);

        for(unsigned int j=i+1; j<allgroups.size(); j++) {
            PeakGroup& grp2 = allgroups[j];
            if (grp2.metaGroupId > 0 ) continue;

            //retention time distance
            float rtdist  = abs(parent->meanRt-grp2.meanRt);
            if (rtdist > maxRtDiff*2 ) continue;

            //retention time overlap
            float rtoverlap = mzUtils::checkOverlap(grp1.minRt, grp1.maxRt, grp2.minRt, grp2.maxRt );
            if (rtoverlap < 0.1) continue;

            //peak intensity correlation
            vector<float>peakIntensityB = grp2.getOrderedIntensityVector(samples,PeakGroup::AreaTop);
            float cor = correlation(peakIntensityA,peakIntensityB);
            if (cor < minSampleCorrelation) continue;

            //peak shape correlation
            float cor2 = largestSample->correlation(grp1.meanMz,grp2.meanMz,ppm,grp1.minRt,grp1.maxRt);
            if (cor2 < minPeakShapeCorrelation) continue;

            //passed all the filters.. group grp1 and grp2 into a single metagroup
            //cerr << rtdist << " " << cor << " " << cor2 << endl;
            grp2.metaGroupId = grp1.metaGroupId;
        }
    }
}

Peak* PeakGroup::getHighestIntensityPeak() {
    Peak* highestIntensityPeak=0;
    for(int i=0; i < peaks.size(); i++ ) {
        if (!highestIntensityPeak or peaks[i].peakIntensity > highestIntensityPeak->peakIntensity) {
           highestIntensityPeak = &peaks[i];
        }
    }
    return highestIntensityPeak;
}


int PeakGroup::getChargeStateFromMS1(float ppm) {
    Peak* highestIntensityPeak=this->getHighestIntensityPeak();
    int scanum = highestIntensityPeak->scan;
    Scan* s = highestIntensityPeak->getSample()->getScan(scanum);

    int peakPos=0;
    if (s) peakPos = s->findHighestIntensityPos(highestIntensityPeak->peakMz,ppm);

    if (s and peakPos > 0 and peakPos < s->nobs() ) {
        vector<int> parentPeaks = s->assignCharges(ppm);
        return parentPeaks[peakPos];
    } else {
        return 0;
    }
}

bool PeakGroup::isMonoisotopic( float ppm) {
    Peak* highestIntensityPeak=this->getHighestIntensityPeak();
    if(! highestIntensityPeak) return false;

    int scanum = highestIntensityPeak->scan;
    Scan* s = highestIntensityPeak->getSample()->getScan(scanum);
    return s->isMonoisotopicPrecursor(highestIntensityPeak->peakMz,ppm, this->chargeState);
}

bool PeakGroup::isGroupGood() {
    return isGroupLabeled('g');
}

bool PeakGroup::isGroupBad() {
    return isGroupLabeled('b');
}

bool PeakGroup::isGroupLabeled(char label) {
    return find(labels.begin(), labels.end(), label) != labels.end();
}

void PeakGroup::markGroupGood() {
    if (isGroupGood()) return;
    labels.erase(remove(labels.begin(), labels.end(), 'b'), labels.end());
    labels.push_back('g');
}

void PeakGroup::markGroupBad() {
    if (isGroupBad()) return;
    labels.erase(remove(labels.begin(), labels.end(), 'g'), labels.end());
    labels.push_back('b');

}

void PeakGroup::toggleLabel(char label) {
    processLabel(label, true);
}

void PeakGroup::addLabel(char label) {
    processLabel(label, false);
}

void PeakGroup::processLabel(char label, bool isToggle) {
    if (label == '\0') {
        //reserved character: clear all labels

        labels.clear();

    } else if (find(labels.begin(), labels.end(), label) != labels.end()) {

        if (isToggle){
        //If the label already exists, remove it

            labels.erase(remove(labels.begin(), labels.end(), label), labels.end());

        } else {
        //If the label already exists, do nothing

           return;
        }
    } else if (label == 'g') {
        //reserved character: 'good', add mark and remove 'bad' label

        labels.erase(remove(labels.begin(), labels.end(), 'b'), labels.end());
        labels.push_back('g');

    } else if (label == 'b') {
        //reserved character: 'bad', add mark and remove 'good' label

        labels.erase(remove(labels.begin(), labels.end(), 'g'), labels.end());
        labels.push_back('b');


    } else {
        //non-reserved character label that is not already in vector: append to vector

        labels.push_back(label);
    }
}

string PeakGroup::getPeakGroupLabel() {
    return string(labels.begin(), labels.end());
}
