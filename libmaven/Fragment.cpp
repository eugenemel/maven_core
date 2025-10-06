#include "Fragment.h"
#include "mzSample.h"

#include "directinfusionprocessor.h"

//empty constructor
Fragment::Fragment() { 
	precursorMz = 0; 
	polarity=1; 
	scanNum=0; 
	rt=0;
   	collisionEnergy=0; 
    consensus=nullptr;
	precursorCharge=0; 
	isDecoy=false; 
	sortedBy=None; 
    group=nullptr;
    mergeCount=0;
    purity=0;
    clusterId=0;
	mergedScore=0;
    scanNumMap={};
    isNLSpectrum=false;
}


/**
 * @brief Fragment::Fragment
 * from a single MS2 scan, create an adjusted MS2 spectrum to compare to library spectra.
 *
 * @param scan
 * @param minFracIntensity: retain peaks with intensity at or above this proportion of max intensity peak
 * @param minSNRatio: retain peaks with S:N at or above this value (use @param baseLineLevel for noise level)
 * @param maxNumberOfFragments: maximum number of spectral peaks to include
 * @param baseLinePercentile: (expressed as a percentage) intensity percentile corresponding to S:N ratio of 1
 * @param isRetainFragmentsAbovePrecursorMz: if true, keep peaks above precursor m/z. If false, skip peaks
 * Issue 570: Note that if the precursorMz value is set to 0, all peaks will be excluded.
 * @param precursorPurityPpm: ppm value used to determine precursor purity (search MS1 scans)
 * @param minIntensity: minimum intensity value for a peak to be retained
 */
Fragment::Fragment(Scan* scan,
                   float minFracIntensity,
                   float minSNRatio,
                   int maxNumberOfFragments,
                   int baseLinePercentile,
                   bool isRetainFragmentsAbovePrecursorMz,
                   float precursorPurityPpm,
                   float minIntensity
                   ) {
    this->precursorMz = scan->precursorMz;
    this->collisionEnergy = scan->collisionEnergy;
    this->polarity = scan->getPolarity();
    this->sampleName = scan->getSampleName();
    this->scanNum = scan->scannum;
    this->precursorCharge = scan->precursorCharge;
    this->group = nullptr;
    this->consensus = nullptr;
	this->mergeCount=0;
	this->mergedScore=0;
	this->clusterId=0;
    this->scanNumMap={};
    scanNumMap.insert(make_pair(scan->getSampleName(), unordered_set<int>()));
    scanNumMap[scan->getSampleName()].insert(scan->scannum);

    //faster implementation when fewer filters specified

    bool isKeepAllFragments = maxNumberOfFragments < 0 || maxNumberOfFragments >= static_cast<int>(scan->nobs());

    if (minFracIntensity <= 0 && minSNRatio <= 0 && isKeepAllFragments && baseLinePercentile <= 0 && isRetainFragmentsAbovePrecursorMz){
        if (minIntensity <= 0) {
            this->mzs = scan->mz;
            this->intensity_array = scan->intensity;
        } else {
            for (unsigned int i = 0; i < scan->nobs(); i++) {
                if (scan->intensity[i] >= minIntensity) {
                    this->mzs.push_back(scan->mz[i]);
                    this->intensity_array.push_back(scan->intensity[i]);
                }
            }
        }

        this->sortedBy = SortType::Mz;
    } else {

        //<intensity, mz>
        vector<pair<float,float> >mzarray = scan->getTopPeaks(minFracIntensity, minSNRatio, baseLinePercentile, minIntensity);

        if (maxNumberOfFragments < 0) {
            maxNumberOfFragments = INT_MAX;
        }

        for (unsigned int j=0; j<mzarray.size() && j < static_cast<unsigned int>(maxNumberOfFragments); j++ ) {
            if (isRetainFragmentsAbovePrecursorMz || mzarray[j].second <= this->precursorMz) { //remove fragments higher than precursorMz
                this->mzs.push_back(mzarray[j].second);
                this->intensity_array.push_back(mzarray[j].first);
            }
        }

        this->sortedBy = SortType::None;
    }

    this->obscount = vector<int>( this->mzs.size(), 1);
    this->fragment_labels = vector<string>(this->mzs.size(), "");

    this->rt = scan->rt;
    this->purity = precursorPurityPpm > 0 ? scan->getPrecursorPurity(precursorPurityPpm) : 0.0f;  //this might be slow
    this->isNLSpectrum = false;
    this->sortByMz();
}

/**
 * @brief Fragment::Fragment
 * @param scan
 *
 * The data from the scan is directly copied into the Fragment, with no filtering or adjusting.
 */
Fragment::Fragment(Scan *scan){
    this->precursorMz = scan->precursorMz;
    this->collisionEnergy = scan->collisionEnergy;
    this->polarity = scan->getPolarity();
    this->sampleName = scan->getSampleName();
    this->scanNum = scan->scannum;
    this->precursorCharge = scan->precursorCharge;
    this->group = nullptr;
    this->consensus = nullptr;
    this->rt = scan->rt;
    this->purity = 0;
    this->mergeCount=0;
    this->mergedScore=0;
    this->clusterId=0;
    this->mzs = scan->mz;
    this->intensity_array = scan->intensity;
    this->fragment_labels = vector<string>(this->mzs.size(), "");
    this->sortedBy = SortType::Mz; // scans should always be encoded in increasing m/z.
    this->obscount = vector<int>( this->mzs.size(), 1); //used when creating consensus spectra.
    this->scanNumMap={};
    scanNumMap.insert(make_pair(scan->getSampleName(), unordered_set<int>()));
    scanNumMap[scan->getSampleName()].insert(scan->scannum);
    this->isNLSpectrum = false;
}

//delete
Fragment::~Fragment() {
    if(consensus) {
        delete(consensus);
        consensus = nullptr;
    }
    mzUtils::delete_all(brothers);
}

//make a copy of Fragment.
Fragment::Fragment( Fragment* other) { 
    this->precursorMz = other->precursorMz;
    this->polarity = other->polarity;
    this->mzs = other->mzs;
    this->intensity_array = other->intensity_array;
    this->fragment_labels = other->fragment_labels;
    this->obscount = other->obscount;
    this->consensus = other->consensus;
    this->scanNum = other->scanNum;
    this->sampleName = other->sampleName;
    this->collisionEnergy = other->collisionEnergy;
    this->precursorCharge= other->precursorCharge;
    this->sortedBy = other->sortedBy;
    this->group = other->group;
    this->mergeCount = other->mergeCount;
    this->purity = other->purity;
    this->tmtQuant = other->tmtQuant;
    this->clusterId= other->clusterId;
	this->mergedScore = other->mergedScore;
    this->scanNumMap = other->scanNumMap;
    this->consensusPositionToScanIntensities = other->consensusPositionToScanIntensities;
    this->isNLSpectrum = other->isNLSpectrum;
}

void Fragment::appendBrothers(Fragment* other) {
    //copy other's brothers
    for(unsigned int i=0; i < other->brothers.size(); i++ ) {
        this->brothers.push_back( other->brothers[i]);
    }

    //empty brothers link
    other->brothers.clear();

    //append other to brother link
    this->brothers.push_back(other);
}

void Fragment::printMzList() { 
    for(unsigned int i=0; i<mzs.size(); i++ ) { cerr << setprecision(3) << mzs[i] << " "; }
}

int Fragment::findClosestHighestIntensityPos(float _mz, float ppmTolr) {
    float tolr = _mz/1e6*ppmTolr;
    float mzmin = _mz - tolr;
    float mzmax = _mz + tolr+0.001;

    vector<float>::iterator itr = lower_bound(mzs.begin(), mzs.end(), mzmin-1);
    int lb = itr-mzs.begin();
    float highestIntensity=0; 
    for(unsigned int k=lb; k < mzs.size(); k++ ) {
        if (mzs[k] < mzmin) continue; 
        if (mzs[k] > mzmax) break;
        if (intensity_array[k] > highestIntensity) highestIntensity=intensity_array[k];
    }

    int bestPos=-1; float bestScore=0;
    for(unsigned int k=lb; k < mzs.size(); k++ ) {
        if (mzs[k] < mzmin) continue; 
        if (mzs[k] > mzmax) break;
        float deltaMz = (mzs[k]-_mz); 
        float alignScore = sqrt(intensity_array[k] / highestIntensity)-(deltaMz*deltaMz);
        //	cerr << _mz << "\t" << k << "\t" << deltaMz << " " << alignScore << endl;
        if (bestScore < alignScore) { bestScore=alignScore; bestPos=k; }
    }
    //if(bestPos>=0) cerr << "best=" << bestPos << endl;
    return bestPos;
}


void Fragment::printFragment(float productPpmToll,unsigned int limitTopX=10) {
    cerr << setprecision(10) << " preMz=" << precursorMz << " ce=" << this->collisionEnergy <<  " scan=" << this->scanNum << endl;
    cerr << " mzs: \t";
    for(unsigned int i=0; i<mzs.size() && i < limitTopX; i++ ) {
        cerr << setprecision(6) << mzs[i] << " ";
        if (annotations.count(i))  {
            string ionName = annotations[i];
            cerr  << "[" << ionName << "] "; 
        }
        if (!fragment_labels[i].empty()) {
            cerr << "[" << fragment_labels[i] << "]";
        }

        //cerr  << "[" << (int) intensity_array[i] << "] "; 
        // cerr << "[" << obscount[i] << "] ";
    }
    cerr << endl;

    for(unsigned int i=0; i<brothers.size(); i++) {
        double matchScore = this->compareToFragment(brothers[i],productPpmToll);
        cerr << setprecision(3) << " similarity=" << matchScore;

        cerr << " brother mzs=\t";
        for(unsigned int j=0; j<brothers[i]->mzs.size() && j< limitTopX; j++ ) {
            cerr << setprecision(3) << brothers[i]->mzs[j] << " ";
            //ycerr << "[" << (int) brothers[i]->intensity_array[j] << "] ";
            //cerr  << "[" << brothers[i]->obscount[j] << "] ";
        }
        cerr << endl;
    }

    if (this->consensus) { 
        cerr << "Cons mzs=\t";
        for(unsigned int j=0; j<this->consensus->mzs.size(); j++ )  cerr << setprecision(3) << this->consensus->mzs[j] << " "; 
        cerr << endl;
    }
}

void Fragment::printConsensusMS2(ostream& outstream, double minConsensusFraction) {

    if(this->consensus == NULL) return; 
    int consensusSize = this->brothers.size()+1;

    outstream << "S" << "\t" << consensus->scanNum << "\t" << consensus->scanNum << "\t" << consensus->precursorMz << endl;

    for(unsigned int i=0; i<consensus->mzs.size(); i++ ) { 
        float fracObserved = ((float) consensus->obscount[i])/consensusSize;
        if (fracObserved > minConsensusFraction )  {
            outstream << setprecision(5) << consensus->mzs[i] << "\t";
            outstream << consensus->intensity_array[i] << endl;
        }
    }
    outstream << endl;
}


void Fragment::printConsensusMGF(ostream& outstream, double minConsensusFraction) {

    if(this->consensus == NULL) return;
    int consensusSize = this->brothers.size()+1;
    int precursorCharge=this->precursorCharge;



    outstream << "BEGIN IONS" << endl;
    if (!sampleName.empty()) { outstream << "TITLE=" <<  sampleName << "." << consensus->scanNum << "." << consensus->scanNum << "." << precursorCharge << endl; }
    outstream << "PEPMASS=" << setprecision(8) << precursorMz << endl;
    outstream << "RTINSECONDS=" << setprecision(9) << rt*60 << "\n";
    outstream << "CHARGE=" << precursorCharge; if(polarity < 0) outstream << "-"; else outstream << "+"; outstream << endl;

    for(unsigned int i=0; i<consensus->mzs.size(); i++ ) {
        float fracObserved = ((float) consensus->obscount[i])/consensusSize;
        if (fracObserved > minConsensusFraction )  {
            outstream << setprecision(8) << consensus->mzs[i] << "\t";
            outstream << consensus->intensity_array[i] << endl;
        }
    }
    outstream << "END IONS" << endl;
}

float Fragment::consensusRt() { 

	if (brothers.size() == 0) 
		return this->rt;

    //compute retention time window
    StatisticsVector<float>retentionTimes; retentionTimes.push_back(this->rt);
    for(unsigned int i=0; i<brothers.size();i++ ) retentionTimes.push_back(brothers[i]->rt);
	return retentionTimes.mean();
}

float Fragment::consensusPurity() {

    if (brothers.size() == 0)
        return this->purity;

    //compute average purity
    StatisticsVector<float>v; v.push_back(this->purity);
    for(unsigned int i=0; i<brothers.size();i++ ) v.push_back(brothers[i]->purity);
    return v.mean();
}


void Fragment::printConsensusNIST(ostream& outstream, double minConsensusFraction, float productPpmToll, Compound* compound=0, Adduct* adduct=0) {

	//this->buildConsensus(productPpmToll);
    //if(this->consensus == NULL) return;

    string COMPOUNDNAME = "Unknown";
    if (compound)  COMPOUNDNAME = compound->id;

    //compute average MVH score
    float avgMVH=0;

    //compute retention time window
    StatisticsVector<float>retentionTimes; retentionTimes.push_back(this->rt);
    for(unsigned int i=0; i<brothers.size();i++ ) retentionTimes.push_back(brothers[i]->rt);

    float precursorMz =  this->precursorMz;

    //scan list
    string scanList = this->sampleName + "." + integer2string(this->scanNum);
    for(unsigned int i=0; i<brothers.size();i++ )  scanList += ";" + brothers[i]->sampleName + "." + integer2string(brothers[i]->scanNum);

    //int consensusSize = this->brothers.size()+1;
    int consensusSize = (this->mergeCount+1);

	float avgRt = this->consensusRt();

    if (compound) {
        if (adduct) {
            outstream << "Name: " << compound->id << " " << adduct->name << endl;
        } else {
            outstream << "Name: " << compound->id << endl;
        }

        outstream << "Id: " << compound->id << endl;
        if (adduct) {
            outstream << "ADDUCT: " << adduct->name << endl;
        }

        if (!compound->formula.empty()) outstream << "FORMULA: " << compound->formula << endl;
        if (compound->category.size())    outstream << "CATEGORY: " << compound->category.front() << endl;
        if (!compound->smileString.empty()) outstream << "SMILE: " << compound->smileString << endl;
        outstream << "LOGP: " << compound->logP << endl;

    } else {
            outstream << "Name: " << "unknown" << ppmround(this->precursorMz,1000) << "@" << ppmround(avgRt,10) << endl;
	}

    outstream << "PrecursorMZ: " << setprecision(10) << precursorMz << endl;
    outstream << "RT: " << setprecision(5) << avgRt << endl;
    outstream << "Comment: ";
    outstream << " Spec=Consensus";
    outstream << " Parent=" << setprecision(10) << precursorMz;
    if(this->collisionEnergy) outstream << " collisionEnergy=" << (int)this->collisionEnergy;
    if(avgMVH) outstream << " AvgMVH=" << avgMVH;
    outstream << " AvgRt=" << avgRt;

    //outstream << " Scanlist=" << scanList;
    outstream << " ConsensusSize=" << consensusSize << endl;
    outstream << "NumPeaks: " << mzs.size() << endl;


    for(unsigned int i=0; i< this->mzs.size(); i++ ) {
        float fracObserved = ((float) this->obscount[i])/consensusSize;

        if (fracObserved > minConsensusFraction )  {

            outstream << setprecision(7) << this->mzs[i] << " ";
            outstream << (int) this->intensity_array[i] << " ";

            string ionName = "?";
            if (annotations.count(i)) ionName = annotations[i];
            if (!fragment_labels[i].empty()) ionName = fragment_labels[i];
            outstream << "\"" << ionName << " " << this->obscount[i]<< "\"" << endl;
        }
    }
    outstream << endl;
}

void Fragment::printInclusionList(bool printHeader, ostream& outstream, string COMPOUNDNAME) {
    if(this->consensus == NULL) return; 
    bool printAvgRtTime=true;

    if (printHeader) { 
        outstream << "mz,ionizationmode,rtmin,rtmax,ce,charge,rt,comment\n";
    }

    //compute retention time window
    StatisticsVector<float>retentionTimes; retentionTimes.push_back(this->rt);
    for(unsigned int i=0; i<brothers.size();i++ ) retentionTimes.push_back(brothers[i]->rt);

    float medRt = retentionTimes.median();
    float stdRt	 = sqrt(retentionTimes.variance());

    if (stdRt > 1 ) stdRt = 1;
    if (stdRt <= 0 ) stdRt = 0.5;
    float minRt =  medRt - 3*stdRt;
    float maxRt =  medRt + 3*stdRt;

    //start otput

    outstream << setprecision(10) << consensus->precursorMz << ",";
    this->polarity > 0 ? outstream << "Positive," : outstream << "Negative,";
    outstream << minRt << ",";
    outstream << maxRt << ",";

    if (this->collisionEnergy ) {
        outstream << this->collisionEnergy << ","; 	//CE
    } else {
        outstream << 25 << ","; 	//default CE set to 25
    }
    outstream << "1,";

    if (printAvgRtTime) { 
        outstream << medRt << ",";
    }

    outstream << COMPOUNDNAME << ":" << setprecision(8) << consensus->precursorMz << " cid:" << consensus->collisionEnergy << endl;

}

/**
 * @brief Fragment::scoreMatch
 * @param other - observed MS/MS spectrum ('this' is the library MS/MS spectrum)
 * @param productPpmTolr
 * @return
 */
FragmentationMatchScore Fragment::scoreMatch(Fragment* other, float productPpmTolr) {
    FragmentationMatchScore s;
    if (mzs.size() < 2 or other->mzs.size() < 2) return s;

    //which one is smaller;
    Fragment* a = this;
    Fragment* b = other;

    s.ppmError = abs((a->precursorMz-b->precursorMz)/a->precursorMz*1e6);

    //use library precursorMz if it exists.  If it is not provided, use the experimental precursorMz
    double precursorMz = a->precursorMz > 0 ? a->precursorMz : b->precursorMz;

    float maxDeltaMz = (productPpmTolr * static_cast<float>(precursorMz))/ 1000000;

    /*
     * ranks[x] = y
     * x = index of frag peak in a
     * y = index of frag peak in b
     */
    s.ranks = findFragPairsGreedyMz(a, b, maxDeltaMz);

    for(int rank: s.ranks) { if(rank != -1) s.numMatches++; }

//    cerr << "Num Matches=" << s.numMatches << endl;

    //annotate?
    for(int i=0; i < s.ranks.size(); i++){
        if (s.ranks[i] != -1) {
            other->annotations[s.ranks[i]]=this->annotations[i];
        }
    }

    s.fractionMatched = s.numMatches / a->nobs();
    s.spearmanRankCorrelation = spearmanRankCorrelation(s.ranks);
    s.ticMatched = ticMatched(s.ranks);
    s.mzFragError =  mzErr(s.ranks,b);
    s.dotProduct = dotProduct(b);
    s.hypergeomScore  = SHP(s.numMatches,a->nobs(),b->nobs(),100000) + s.ticMatched; // ticMatch is tie breaker
    s.mvhScore = MVH(s.ranks,b);
    s.weightedDotProduct = mzWeightedDotProduct(s.ranks,b);
    s.matchedQuantiles = matchedRankVector(s.ranks,b);
    //s.dotProductShuffle = this->dotProductShuffle(b,2000);

    s.dotProductNorm = normCosineScore(this, b, s.ranks);
    s.dotProductNormNL = normNLCosineScore(this, b, productPpmTolr);

    //cerr << "scoreMatch:\n" << a->nobs() << "\t" << b->nobs() << "\t" << s.numMatches << " hyper=" << s.hypergeomScore << "\n";

    return s;
}

map<string, int> Fragment::getDiagnosticMatches(vector<string>& labels, vector<int>& ranks){

    map<string, int> diagnosticMatches = {};

    for (auto x : labels) {
        diagnosticMatches.insert(make_pair(x, 0));
    }

    if (annotations.size() != mzs.size()) return diagnosticMatches;

    for (unsigned int i = 0; i < ranks.size(); i++){

        if (annotations[i].empty()) continue;

        if (ranks[i] != -1) { //indicates a matched fragment

            for (auto label : labels){
                if (annotations[i].find(label) == 0) diagnosticMatches[label]++;
            }
        }
    }

    return diagnosticMatches;
}

double Fragment::compareToFragment(Fragment* other, float productPpmTolr) {
    if (mzs.size() < 2 or other->mzs.size() < 2) return 0; 	
    //which one is smaller;
    Fragment* a = this;
    Fragment* b =  other;
    //if (b->mzs.size() < a->mzs.size() ) { a=other; b=this; }
    vector<int>ranks = compareRanks(a,b,productPpmTolr);
    //return spearmanRankCorrelation(ranks);
    //return fractionMatched(ranks);
    return ticMatched(ranks);
}

/**
 * Note change to use findFragPairsGreedyMz()
 */
vector<int> Fragment::compareRanks(Fragment* a, Fragment* b, float productPpmTolr) {
    float maxDeltaMz = (productPpmTolr * static_cast<float>(a->precursorMz))/ 1000000;
    return findFragPairsGreedyMz(a, b, maxDeltaMz);
}

/**
 * Method to return all m/z matches between two MS/MS spectra (represented as Fragment objects).
 *
 * Requires that Fragments be sorted by m/z.  If they are not, they will be in this method.
 *
 * @brief Fragment::findMatches
 * @param a (MS/MS spectrum 1) (typically from library)
 * @param b (MS/MS spectrum 2) (typically from experimental data)
 * @param maxMzDiff (tolerance, translated into a maximum m/z difference)
 * @param isMatchToNLMassShift: Match to NL-adjusted m/zs (precursorMz).
 * If both an NL-adjusted m/z and regular m/z match, the lowest m/z diff is retained.
 * [default=false]
 *
 * @return ranks vector<int>, of length of a.
 *
 * rank[a_position] = b_position
 * when sorted by m/z
 *
 * Note that these pairs are only valid if the two MS/MS spectra remain sorted by Mz.
 */
vector<int> Fragment::findFragPairsGreedyMz(Fragment* library, Fragment* observed, float maxMzDiff, bool isMatchToNLMassShift) {

    //Sort spectra by m/z
    library->sortByMz();
    observed->sortByMz();

    map<unsigned int, float> NLMzs{};
    if (isMatchToNLMassShift) {
        for (unsigned int i = 0; i < library->mzs.size(); i++) {
            double libraryMz = static_cast<double>(library->mzs[i]);
            if (libraryMz < library->precursorMz) {
                float nLMz = static_cast<float>(library->precursorMz - libraryMz);
                NLMzs.insert(make_pair(i, nLMz));
            }
        }
    }

    //TODO: need to intermix NLMzs with library mzs.
    //TODO: use special class, mapping int to some special container class?

    //Identify all valid possible fragment pairs (based on tolerance),
    //record with mzDelta

    vector<pair<float,pair<unsigned int, unsigned int>>> fragPairsWithMzDeltas;

    for (unsigned int i = 0; i < library->mzs.size(); i++){

        float mzLibraryFrag = library->mzs.at(i);

        for (unsigned int j = 0; j < observed->mzs.size(); j++) {

            float mzObservedFrag = observed->mzs.at(j);

            if (mzLibraryFrag - mzObservedFrag > maxMzDiff) {
                //Out of tolerance - could possibly come into tolerance as j increases.
            } else if (mzObservedFrag - mzLibraryFrag > maxMzDiff) {
                //Out of tolerance - cannot possibly come into tolerance as j increases.
                break;
            } else {
                //In tolerance - record dissimilarity as candidate match.

                float mzDelta = abs(mzLibraryFrag - mzObservedFrag);

                pair<unsigned int, unsigned int> peakPair (i, j); //First position is reserved for a, second for b

                pair<float,pair<unsigned int, unsigned int>> fragPairWithMzDelta (mzDelta, peakPair);

                fragPairsWithMzDeltas.push_back(fragPairWithMzDelta);
            }

        }
    }

    //sort standard ions by m/z for matching.
    sort(fragPairsWithMzDeltas.begin(), fragPairsWithMzDeltas.end(),
         [ ](const pair < float,pair < int, int > > & lhs, const pair < float,pair < int, int > > & rhs){
           if (abs(static_cast<double>(lhs.first) - static_cast<double>(rhs.first)) < 1e-6) {
             if (lhs.second.first == rhs.second.first) {
               return lhs.second.second < rhs.second.second;
             } else {
               return lhs.second.first < rhs.second.first;
             }
           } else {
             return lhs.first < rhs.first;
           }
         });

    //Once a fragment has been claimed in a frag pair, it may not be involved in any other
    //frag pair.
    vector<pair<unsigned int, unsigned int>> matches;
    vector<int> ranks (library->mzs.size(),-1);
    set<unsigned int> claimedLibraryFrags;
    set<unsigned int> claimedObservedFrags;

//    cerr << "match: ref <--> obs" << endl;

    for (pair<float,pair<unsigned int, unsigned int>> fragPairWithMzDelta : fragPairsWithMzDeltas){
        pair<unsigned int, unsigned int> fragPair = fragPairWithMzDelta.second;
        unsigned int libraryFrag = fragPair.first;
        unsigned int observedFrag = fragPair.second;

        if (claimedLibraryFrags.count(libraryFrag) == 0 && claimedObservedFrags.count(observedFrag) == 0){

            matches.push_back(fragPair);
            ranks[libraryFrag] = static_cast<int>(observedFrag);

//            cerr << "match: " << a->mzs.at(a_frag) << " <--> " << b->mzs.at(b_frag) << endl;

            claimedLibraryFrags.insert(libraryFrag);
            claimedObservedFrags.insert(observedFrag);
        }
    }

    return ranks;
}


vector<int> Fragment::locatePositions( Fragment* a, Fragment* b, float productPpmTolr) {
    bool verbose=false;
    if (verbose) { cerr << "\t\t "; a->printMzList(); cerr << "\nvs \n"; b->printMzList();  }
    b->sortByMz();
    vector<int> ranks (a->mzs.size(),-1);	//missing value == -1
    for(unsigned int i=0; i<a->mzs.size(); i++ ) {
        int pos = b->findClosestHighestIntensityPos(a->mzs[i],productPpmTolr);
        ranks[i] = pos;
    }
    if (verbose) { cerr << " compareranks: "; for(unsigned int i=0; i < ranks.size(); i++ ) cerr << ranks[i] << " "; }
    return ranks;
}

void Fragment::mergeFragment(Fragment* brother, float productPpmTolr) {

	vector<int>ranks=compareRanks(brother,this,productPpmTolr);
    for(unsigned int j=0; j<ranks.size(); j++ ) {
            int   posA = ranks[j];
            float mzB = brother->mzs[j];
            float intB = brother->intensity_array[j];
            if (posA >= 0)  {
                this->intensity_array[ posA ] += intB;
                this->obscount[ posA ] += 1;
            } else if ( posA == -1 ) {
                this->mzs.push_back(mzB);
                this->intensity_array.push_back(intB);
                this->obscount.push_back(1);
            }
      }

	this->mergeCount++;
}

void Fragment::truncateTopN(int n) { 
	if (this->nobs() < n) return;

	this->sortByIntensity();
	this->mzs.resize(n);
	this->intensity_array.resize(n);
    this->fragment_labels.resize(n);
	this->obscount.resize(n);
}


/**
 * @brief Fragment::buildConsensus
 * @param productPpmTolr: m/z tolerance for associating peaks together
 * @param consensusIntensityAgglomerationType: method to agglomerate intensities
 * @param isIntensityAvgByObserved: intensity averaged by only scans where peak observed (instead of all possible scans)
 * @param isNormalizeIntensityArray: scale all intensities so that max intensity = 10000
 * @param minNumMs2ScansForConsensus: retain peaks found in at least this many scans
 * @param minFractionMs2ScansForConsensus: retain peaks found in at least this proportion of scans
 * @param  bool isRetainOriginalScanIntensities: retain map of <consensus position, vector of scan intensities>
 * @param mzToRemove: vector<float> of m/z values that are removed from the scan after all other processing steps
 * @param mzToRemoveTol: ppm tolerance to remove m/z values from scan (based on proximity in m/z)
 */
void Fragment::buildConsensus(float productPpmTolr,
                              ConsensusIntensityAgglomerationType consensusIntensityAgglomerationType,
                              bool isIntensityAvgByObserved,
                              bool isNormalizeIntensityArray,
                              int minNumMs2ScansForConsensus,
                              float minFractionMs2ScansForConsensus,
                              bool isRetainOriginalScanIntensities,
                              vector<float> mzToRemove,
                              float mzToRemoveTol) {

    if(this->consensus) {
        delete(this->consensus);
        this->consensus=nullptr;
    }

    Fragment* seed= this;

    //find brother with largest nobs
    for(Fragment* b: brothers) if (b->nobs() > seed->nobs()) seed = b;

    //Issue 213: if the seed was changed, reassign brothers appropriately.
    vector<Fragment*> brothers = this->brothers;
    if (seed != this) {
        brothers.erase(remove(brothers.begin(), brothers.end(), seed), brothers.end());
        brothers.push_back(this);
    }

    Fragment* Cons = new Fragment(seed);  //make a copy of self
    this->consensus = Cons;
    Cons->sortByMz();

    unsigned long N = 1 + brothers.size();

    //Issue 217
    map<int, vector<float>> posToIntensityMap{};
    if (isRetainOriginalScanIntensities || consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean){
        for (unsigned int i = 0; i < Cons->intensity_array.size(); i++){
            posToIntensityMap.insert(make_pair(i, vector<float>{Cons->intensity_array[i]}));
        }
    }

    //Issue 245: ensure that a precursor m/z is set for all Fragment* where the m/z is missing.
    //Will always be a problem for MS1 scans.
    double precursorMzIfMissing = 0.0;
    if (!seed->mzs.empty()) precursorMzIfMissing = static_cast<double>(seed->mzs[seed->mzs.size()-1]);

    for (auto brother : brothers) {
        if (!brother->mzs.empty()){
            double candidatePrecursorMz = static_cast<double>(brother->mzs[brother->mzs.size()-1]);
            if (candidatePrecursorMz > precursorMzIfMissing) {
                precursorMzIfMissing = candidatePrecursorMz;
            }
        }
    }

    if (seed->precursorMz <= 0.0) seed->precursorMz = precursorMzIfMissing;

    for (auto brother : brothers) {
        if (brother->precursorMz <= 0.0) brother->precursorMz = precursorMzIfMissing;
    }

    for(unsigned int i=0; i<brothers.size(); i++) {

        Fragment* brother = brothers[i];
        vector<int>ranks=compareRanks(brother,Cons,productPpmTolr);	//location

        for(unsigned int j=0; j<ranks.size(); j++ ) {
            int   posA = ranks[j];
            float mzB = brother->mzs[j];
            float intB = brother->intensity_array[j];
            if (posA >= 0)  {
                Cons->intensity_array[ posA ] += intB;
                Cons->obscount[ posA ] += 1;
                Cons->fragment_labels[posA] = brother->fragment_labels[j];

                //Issue 217
                if (isRetainOriginalScanIntensities || consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean) {
                    posToIntensityMap[posA].push_back(intB);
                }

            } else if ( posA == -1 ) {
                Cons->mzs.push_back(mzB);
                Cons->intensity_array.push_back(intB);
                Cons->obscount.push_back(1);
                Cons->fragment_labels.push_back(brother->fragment_labels[j]);

                //Issue 217
                if (isRetainOriginalScanIntensities || consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean){
                    posToIntensityMap.insert(make_pair(Cons->mzs.size()-1,vector<float>{intB}));
                }
            }
        }

        //Issue 217: restructure map based on sorting
        if (isRetainOriginalScanIntensities || consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean){

            vector<int> mzOrderInc = Cons->mzOrderInc();

            map<int, vector<float>> posToIntensityMapSorted{};

            for(unsigned int j = 0; j < mzOrderInc.size(); j++) {

                int sortedPos = j;
                int unsortedPos = mzOrderInc[j];

                auto entry = posToIntensityMap.find(unsortedPos);
                if (entry != end(posToIntensityMap)) {
                    auto const value = std::move(entry->second);
                    posToIntensityMapSorted.insert({sortedPos, std::move(value)});
                } else {
                    cerr << "Illegal mapping in Fragment::buildConsensusSpectrum() for consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean!" << endl;
                    cerr << "unsortedPos was not found in posToIntensityMap." << endl;
                    cerr << "Exiting program." << endl;
                    abort();
                }

            }

            posToIntensityMap = posToIntensityMapSorted;
        }

        Cons->sortedBy = SortType::None;
        Cons->sortByMz();

        map<string, unordered_set<int>> brotherMap = brother->scanNumMap;
        for (auto it = brotherMap.begin(); it != brotherMap.end(); ++it) {

            string sample = it->first;
            unordered_set<int> scans = it->second;

            if (Cons->scanNumMap.find(sample) == Cons->scanNumMap.end()) {
                Cons->scanNumMap.insert(make_pair(sample, unordered_set<int>()));
            }

            for (auto x : scans) {
                Cons->scanNumMap[sample].insert(x);
            }
        }
    }

    for (auto it = this->scanNumMap.begin(); it != scanNumMap.end(); ++it){
        string sample = it->first;
        unordered_set<int> scans = it->second;

        if (Cons->scanNumMap.find(sample) == Cons->scanNumMap.end()) {
            Cons->scanNumMap.insert(make_pair(sample, unordered_set<int>()));
        }

        for (auto x : scans) {
            Cons->scanNumMap[sample].insert(x);
        }
    }

    //compute retention time window
    Cons->rt  = this->consensusRt();

    //compute average purity
    Cons->purity = this->consensusPurity();

    //optionally filter values from scans based on number of observations
    if (minNumMs2ScansForConsensus > 0 || minFractionMs2ScansForConsensus > 0) {

        vector<float> filtered_mzs;
        vector<float> filtered_intensities;
        vector<string> filtered_fragment_labels;
        vector<int> filtered_obscount;

        //Issue 227
        vector<vector<float>> medianIntensities;

        for (unsigned int i = 0; i < Cons->mzs.size(); i++){
            float frac = static_cast<float>(Cons->obscount[i]) / static_cast<float>(N);
            if (Cons->obscount[i] >= minNumMs2ScansForConsensus && frac >= minFractionMs2ScansForConsensus) {
                filtered_mzs.push_back(Cons->mzs[i]);
                filtered_intensities.push_back(Cons->intensity_array[i]);
                filtered_fragment_labels.push_back(Cons->fragment_labels[i]);
                filtered_obscount.push_back(Cons->obscount[i]);

                //Issue 227
                if (consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean){
                    medianIntensities.push_back(posToIntensityMap[static_cast<int>(i)]);
                }
            }
        }

        Cons->mzs = filtered_mzs;
        Cons->intensity_array = filtered_intensities;
        Cons->fragment_labels = filtered_fragment_labels;
        Cons->obscount = filtered_obscount;

        //Issue 227
        if (isRetainOriginalScanIntensities || consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean){

            posToIntensityMap.clear();
            for (unsigned int i = 0; i < medianIntensities.size(); i++){
                posToIntensityMap.insert(make_pair(i, medianIntensities[i]));
            }

        }
    }

    //agglomerate intensity values
    if( brothers.size() >= 1) {
        for( unsigned int i=0; i < Cons->intensity_array.size(); i++){
            if (consensusIntensityAgglomerationType == Mean) {
                Cons->intensity_array[i] /= (isIntensityAvgByObserved ? Cons->obscount[i] : N);
            } else if (consensusIntensityAgglomerationType == Median) {
                Cons->intensity_array[i] = mzUtils::median(posToIntensityMap[static_cast<int>(i)]);
            } else if (consensusIntensityAgglomerationType == Sum) {
                Cons->intensity_array[i] = accumulate(posToIntensityMap[static_cast<int>(i)].begin(), posToIntensityMap[static_cast<int>(i)].end(), 0.0f);
            } else if (consensusIntensityAgglomerationType == Max) {
                auto max_iterator = std::max_element(posToIntensityMap[static_cast<int>(i)].begin(), posToIntensityMap[static_cast<int>(i)].end());
                Cons->intensity_array[i] = *max_iterator;
            }
        }
	}
	
    //normalize intensities (if specified)
    if (isNormalizeIntensityArray) {
        Cons->sortByIntensity();
        if (Cons->intensity_array.size() > 1) {
            float maxValue = Cons->intensity_array[0];
            for(unsigned int i=0; i<Cons->intensity_array.size(); i++)  {
                Cons->intensity_array[i] = Cons->intensity_array[i]/maxValue*10000;
            }
        }
    }

    //Issue 797
    if (!mzToRemove.empty()) {

        // Must be sorted by m/z, otherwise std::lower_bound() returns nonsense
        Cons->sortByMz();

        set<int> indexesToRemove{};

        for (float mz : mzToRemove) {

            const float mzLowerBound = mz - mz*mzToRemoveTol/1e6;
            const float mzUpperBound = mz + mz*mzToRemoveTol/1e6;

            auto it = std::lower_bound(Cons->mzs.begin(), Cons->mzs.end(), mzLowerBound);

            while(it != Cons->mzs.end()) {
                const float mzI = *it;
                size_t pos = std::distance(Cons->mzs.begin(), it);
                if (mzI <= mzUpperBound) {
                    indexesToRemove.insert(pos);
                } else {
                    break;
                }
                it++;
            }
        }

        // update via removing indexes
        Cons->mzs = mzUtils::removeIndexes(Cons->mzs, indexesToRemove);
        Cons->intensity_array = mzUtils::removeIndexes(Cons->intensity_array, indexesToRemove);
        Cons->fragment_labels = mzUtils::removeIndexes(Cons->fragment_labels, indexesToRemove);
        Cons->obscount = mzUtils::removeIndexes(Cons->obscount, indexesToRemove);
        posToIntensityMap = mzUtils::removeIndexesFromMap(posToIntensityMap, indexesToRemove);
    }

    //Issue 468, 532
    if (isRetainOriginalScanIntensities) {
        Cons->consensusPositionToScanIntensities = posToIntensityMap;
    } else {
        Cons->consensusPositionToScanIntensities.clear();
    }

    this->consensus = Cons; 
}

vector<unsigned int> Fragment::intensityOrderDesc() {
    unsigned int nobs = intensity_array.size();
    vector<pair<float,int> > _pairsarray(nobs);
    vector<unsigned int>position(nobs);
    for(unsigned int pos=0; pos < nobs; pos++ ) {
        _pairsarray[pos] = make_pair(intensity_array[pos],pos);
    }

    //reverse sort first key [ ie intensity ]
    sort(_pairsarray.rbegin(), _pairsarray.rend());

    //return positions in order from highest to lowest intenisty
    for(unsigned int i=0; i < _pairsarray.size(); i++) { position[i] = _pairsarray[i].second; }
    return position;
}


vector<int> Fragment::mzOrderInc() {
    unsigned int nobs = mzs.size();
    vector<pair<float,int> > _pairsarray(nobs);
    vector<int>position(nobs);
    for(unsigned int pos=0; pos < nobs; pos++ ) {
        _pairsarray[pos] = make_pair(mzs[pos],pos);
    }

    //forward sort first key
    sort(_pairsarray.begin(), _pairsarray.end());

    //return positions in order of increasing m/z
    for(unsigned int i=0; i < _pairsarray.size(); i++) { position[i] = _pairsarray[i].second; }
    return position;
}


void Fragment::sortByIntensity() { 
    if(sortedBy == Fragment::SortType::Intensity) return; //speedup already sorted

    vector<unsigned int>order = intensityOrderDesc();
    vector<float> a(mzs.size());
    vector<float> b(intensity_array.size());
    vector<int> c(obscount.size());
    vector<string> fragLabels(fragment_labels.size());
    vector<vector<float>> medianIntensities(order.size());

    map<int,string>d;
    for(unsigned int i=0; i<order.size(); i++) {
        b[i] = intensity_array[order[i]];
        a[i] = mzs[order[i]];
        fragLabels[i] = fragment_labels[order[i]];
        if(order[i] < obscount.size()) c[i] = obscount[order[i]];
        if(annotations.count(order[i])>0) d[i] = annotations[order[i]];

        if (!consensusPositionToScanIntensities.empty()){
            medianIntensities.at(i) = consensusPositionToScanIntensities[order[i]];
        }

    }

    mzs = a;
    intensity_array = b;
    obscount = c;
    annotations=d;
    fragment_labels = fragLabels;

    if (!medianIntensities.empty()) {
        consensusPositionToScanIntensities.clear();
        for (unsigned int i = 0; i < medianIntensities.size(); i++){
            consensusPositionToScanIntensities.insert(make_pair(i, medianIntensities[i]));
        }
    }

    sortedBy = Fragment::SortType::Intensity;
}	

void Fragment::sortByMz() {

    if(sortedBy == Fragment::SortType::Mz) return; //speedup already sorted

    vector<int>order = mzOrderInc();
    vector<float> a(mzs.size());
    vector<float> b(intensity_array.size());
    vector<int> c(obscount.size());
    vector<string> fragLabels(fragment_labels.size());
    vector<vector<float>> medianIntensities(order.size());

    map<int,string>d;

    for(unsigned int i=0; i<order.size(); i++) {
        b[i] = intensity_array[order[i]];
        a[i] = mzs[order[i]];
        fragLabels[i] = fragment_labels[order[i]];
        if(order[i] < obscount.size()) c[i] = obscount[order[i]];
        if (annotations.count(order[i])>0) d[i] = annotations[order[i]];

        if (!consensusPositionToScanIntensities.empty()){
            medianIntensities.at(i) = consensusPositionToScanIntensities[order[i]];
        }

    }

    mzs = a;
    intensity_array = b;
    obscount = c;
    annotations=d;
    fragment_labels = fragLabels;

    if (!medianIntensities.empty()) {
        consensusPositionToScanIntensities.clear();
        for (unsigned int i = 0; i < medianIntensities.size(); i++){
            consensusPositionToScanIntensities.insert(make_pair(i, medianIntensities[i]));
        }
    }

    sortedBy = Fragment::SortType::Mz;
}


double Fragment::spearmanRankCorrelation(const vector<int>& X) {
    double d2=0; 
    int N = min((int)nobs(), (int) X.size()); //max elements in the second vector
    int n=0;
    for(int i=0; i<N;i++ ) {	
        //cerr << i << "\t" << n << "\t" << X[i] << endl;
        if (X[i] != -1 ) { //mising values set to average distance
            d2 += (n-X[i])*(n-X[i]);
			n++;
        } else  {
            d2 += (i-N)*(i-N);
		}
        //cerr << "n=" << n << " cor=" << d2 << endl;
        //}
    }
    if(n>1) { 
        //cerr << "n=" << n << "\t" << d2 << endl;
		//double p = 1.00-(6.0*d2)/(n*((n*n)-1));
		double p = 1.00-(6.0*d2)/(n*((n*n)-1));
		return p;
	} else {
		return 0;
	}
}

double Fragment::fractionMatched(const vector<int>& X) {
    if (X.size() == 0) return 0;
    int matchCount=0;
    for(unsigned int i=0; i<X.size();i++ ) {	if (X[i] != -1 ) matchCount++; }
    //if (verbose) { cerr << "\t\t fractionMatched:" << matchCount << endl; }
    return ((double)matchCount / X.size());
}

double Fragment::mzErr(const vector<int>& X, Fragment* other) {
    if (X.size() == 0) return 0;
	int n=0;
    double ERR=0;
    for(unsigned int i=0; i<X.size(); i++ ) if (X[i] != -1) { n++; ERR += POW2( mzs[i] - other->mzs[ X[i] ]);  }
	if (n) { return sqrt(ERR/n); }
	return 0;
}

double Fragment::totalIntensity() {
    double TIC=0;
    for(unsigned int i=0; i<this->nobs(); i++) TIC += this->intensity_array[i];
    return TIC;
}

vector<float> Fragment::asDenseVector(float mzmin, float mzmax, int nbins=2000) { 
	vector<float>v(nbins,0);
	double mzrange = mzmax-mzmin;
	for(int i=0; i < mzs.size(); i++) {
		if(mzs[i]<mzmin or mzs[i]>mzmax) continue;
		int bin = (int) (mzs[i]-mzmin)/mzrange*nbins;
		if(bin>0 and bin<nbins) v[bin] += this->intensity_array[i];
	}
	return v;
}

/**
 * @deprecated
 **/
void Fragment::normalizeIntensity(vector<float>&x, int binSize=100)  { 
    //condition
    for(int i=0; i<x.size(); i+=binSize) {
	float maxI = 0;
    	for(int j=i; j<i+binSize and j<x.size(); j++) if(x[j]>maxI) maxI=x[j];
	if(maxI <=0) continue;
    	for(int j=i; j<i+binSize and j<x.size(); j++) x[j] = x[j] /maxI;
    }
}

double Fragment::dotProductShuffle(Fragment* other, int nbins=1000) {
   
    double thisTIC = totalIntensity();
    double otherTIC = other->totalIntensity();

    if(thisTIC == 0 or otherTIC == 0) return 0;
    vector<float> va = this->asDenseVector(0,2000,nbins);
    vector<float> vb = other->asDenseVector(0,2000,nbins);
    //return mzUtils::crossCorrelationZ(va,vb,0.45);

	float obs_corr = mzUtils::correlation(va,vb);
	float exp_corr = 0;
	float exp_n    = 20;

	//expectation
	for(int i=0; i<exp_n; i++ ) {
		shuffle(vb);
		exp_corr += mzUtils::correlation(va,vb);
	}
	exp_corr /= exp_n;
	float score = (obs_corr - exp_corr)/0.05;
	if (score > 100) score = 100;
	return(score);
}

double Fragment::dotProduct(Fragment* other) {
    double thisTIC = totalIntensity();
    double otherTIC = other->totalIntensity();

    if(thisTIC == 0 or otherTIC == 0) return 0;
    vector<float> va = this->asDenseVector(100,2000,2000);
    vector<float> vb = other->asDenseVector(100,2000,2000);
    return mzUtils::correlation(va,vb);

    /*
    double dotP=0;
    for(unsigned int i=0; i<X.size(); i++ )
        if (X[i] != -1)  dotP += this->intensity_array[i] * other->intensity_array[X[i]];

    //return dotP/sqrt(thisTIC*thisTIC*otherTIC*otherTIC);   //DIP
    return sqrt(dotP)/sqrt(thisTIC*otherTIC);               //SIM
    */
}

double Fragment::mzWeightedDotProduct(const vector<int>& X, Fragment* other) {
    if (X.size() == 0) return 0;
    double thisTIC = 0;
    double otherTIC = 0;

    for(unsigned int i=0; i<this->nobs();  i++) thisTIC +=  mzs[i] * intensity_array[i];
    for(unsigned int j=0; j<other->nobs(); j++) otherTIC += other->mzs[j] * other->intensity_array[j];

    if(thisTIC == 0 or otherTIC == 0) return 0;

    double dotP=0;
    for(unsigned int i=0; i<X.size(); i++ ) {
        int j = X[i];
        if (j != -1)  dotP += mzs[i]*intensity_array[i] * other->mzs[j]*other->intensity_array[j];
    }

    //return dotP/sqrt(thisTIC*thisTIC*otherTIC*otherTIC);   //DIP
    return sqrt(dotP)/sqrt(thisTIC*otherTIC);               //SIM
}



double Fragment::ticMatched(const vector<int>& X) {
    if (X.size() == 0) return 0;
    double TIC = this->totalIntensity();
    double matchedTIC=0;

    for(unsigned int i=0; i<X.size(); i++ ) { 
        if (X[i] != -1) matchedTIC += intensity_array[i];
    }
    /*
       for(int i=0; i<other->nobs();i++ ) TIC += other->intensity_array[i]; 
       for(int i=0; i<X.size();     i++ ) if (X[i] != -1) matchedTIC += other->intensity_array[ X[i] ];
       */
    if (TIC>0) return matchedTIC/TIC; else return 0;
}

double Fragment::logNchooseK(int N,int k) {
    if (N == k || k == 0) return 0;
    if (N == k) return -1;  //not defined..

    double x = k / (double) N;
    return N*x*log(1/x)+(1-x)*log(1/(1-x));
}

double Fragment::SHP(int k, int m, int n, int N=100000) {   //k=matched, m=len1, n=len2
    if (k==0) return 0;
    //if (k > min(m,n)) { return -1; } //error
    if (k > min(m,n))  k = min(m,n);
    double A=logNchooseK(m,k);
    double B=logNchooseK(N-m,n-k);
    double C=logNchooseK(N,n);
    //cerr << " A,B,C " << "\t" << A << "\t" << B << "\t" << C << endl;
    return -(A+B-C);
}

vector<double> Fragment::matchedRankVector(const vector<int>& X, Fragment* other) {
    int n = other->nobs();
    vector<double> Qcuts  = {0.2, 0.5, 0.8, 1.0};
    vector<double> Counts( Qcuts.size(),0);

    if (X.size() == 0 or n == 0 ) return(Counts);

    for(unsigned int i=0; i<X.size(); i++ ) {
	int j = X[i];
        if (j == -1)  continue;
    	for(int qi=0; qi < Qcuts.size(); qi++ ) {
	   if (j < Qcuts[qi]*n ) { Counts[qi]++; break; }
	}
    }

    for(unsigned int i=0; i < Counts.size(); i++ ) Counts[i] /= n;
    return(Counts);
}

double Fragment::MVH(const vector<int>& X, Fragment* other) {
    //other is experimental spectra
    int N = 100000;
    if (X.size() == 0) return 0;
    //int m = this->nobs();
    int n = other->nobs();

    int Ak =   0; int Am=0.2 * n;
    int Bk =   0; int Bm=0.5 * n;
    int Ck =   0; int Cm=0.8 * n;
    int Dk =   0;

    for(unsigned int i=0; i<X.size(); i++ ) {
        int j = X[i];
        if (j == -1)  Dk++;
        else if (j < 0.2*n)   Ak++; //class A matched
        else if (j < 0.5*n)   Bk++; //class B matched
        else if (j < 0.8*n)   Ck++; //class C matched
        else Dk++;                  //class D matched
    }

    if (Ak>Am) Ak=Am;
    if (Bk>Bm) Bk=Bm;
    if (Ck>Cm) Ck=Cm;


    double A=logNchooseK(Am,Ak) + 0.1*logNchooseK(Bm,Bk) + 0.001*logNchooseK(Cm,Ck);
    double B=logNchooseK(N-Am-Bm-Cm,n-Ak-Bk-Ck);
    double C=logNchooseK(N,n);
    //printf("MVH %d %d  %d %d  %d %d  %d -> %e, %e, %e\n", Ak,Am, Bk,Bm, Ck,Cm, Dk,  A,B,C);
    return -(A+B-C);
}

bool Fragment::hasMz(float mzValue, float ppmTolr) {
    for(unsigned int i=0; i < nobs(); i++)  if (mzUtils::ppmDist(mzs[i],mzValue) < ppmTolr) return true;
    return false;
}

bool Fragment::hasNLS(float NLS, float ppmTolr) {
    float mzValue = precursorMz-NLS;
    for(unsigned int i=0; i < nobs(); i++)  if (mzUtils::ppmDist(mzs[i],mzValue) < ppmTolr) return true;
    return false;
}

/**
  * @deprecated
  * Never used, if this should be used, need to revisit the code, see
  * Fragment::convertToNLFragment() for considerations that should be handled here
**/
void Fragment::addNeutralLosses() {
    int N = mzs.size();
    for(int i=0; i < N; i++) {
        float nLMass = mzs[i]-precursorMz;
        float nLIntensity = intensity_array[i];
        if(nLMass < -1) {
            mzs.push_back(nLMass);
            intensity_array.push_back(nLIntensity);
            fragment_labels.push_back(fragment_labels[i]+ " NL");
            obscount.push_back(1);
        }
    }
    sortByIntensity();
}

//Issue 752: for NL scoring
void Fragment::convertToNLSpectrum() {
    //already an NL spectrum - nothing to do
    if (this->isNLSpectrum) return;

    vector<float> updatedMzs;
    vector<float> updatedIntensities;
    vector<string> updatedLabels;

    this->sortByMz();

    for(unsigned int i = 0; i < nobs(); i++) {
        float originalMass = mzs[i];

        //Do not retain any masses above the precursor m/z
        //Add 1 Da for tolerance wiggle room
        if (originalMass > static_cast<float>(precursorMz+1.0)) break;

        //precursor mass is preserved.
        float nLmass = originalMass;
        if (abs(originalMass-static_cast<float>(precursorMz)) > 0.1f){
            nLmass = static_cast<float>(precursorMz)-originalMass;
        }

        updatedMzs.push_back(nLmass);
        updatedIntensities.push_back(intensity_array[i]);

        string originalLabel = fragment_labels[i];
        string updatedLabel = "";
        if (!originalLabel.empty()) {
            string updatedLabel = "NL " + originalLabel;
        }
        updatedLabels.push_back(updatedLabel);
    }

    this->mzs = updatedMzs;
    this->intensity_array = updatedIntensities;
    this->fragment_labels = updatedLabels;

    this->sortedBy = SortType::None;
    this->sortByMz();

    this->isNLSpectrum = true;
}

double Fragment::normNLCosineScore(Fragment *library, Fragment *observed, float productPpmTolr){

    Fragment *observedNL = new Fragment(observed);
    observedNL->convertToNLSpectrum();

    vector<int> ranks = Fragment::compareRanks(library, observedNL, productPpmTolr);

    double normCosineNLScore = Fragment::normCosineScore(library, observedNL, ranks);

    delete(observedNL);

    return normCosineNLScore;
}

//Warning: do not use this with direct infusion code. only LC-MS/MS analysis.
int Fragment::getNumDiagnosticFragmentsMatched(string fragLblStartsWith, vector<string> labels, vector<int> ranks) {

    int numDiagnosticFragmentsMatched = 0;

    if (labels.size() != ranks.size()) return 0;

    for(int i=0; i < ranks.size(); i++){
        if (ranks[i] != -1) {
            if (labels[i].find(fragLblStartsWith) == 0) {
                numDiagnosticFragmentsMatched++;
            }
        }
    }

    return numDiagnosticFragmentsMatched;
}

/**
 * @brief Fragment::normalizeIntensityArray
 * Irreversible transformation, intensity values are converted to normalized value
 * Max intensity value takes on value of @param normValue. All others scaled commensurately.
 *
 * @param normValue
 */
void Fragment::normalizeIntensityArray(float normValue){

    sortByIntensity();
    if (!intensity_array.empty()) {
        float maxValue = intensity_array[0];
        for(unsigned int i=0; i<intensity_array.size(); i++)  {
            intensity_array[i] = intensity_array[i]/maxValue*normValue;
        }
    }

}

/**
 * @brief Fragment::agglomerateMzs
 * Irreversible transformation, where m/z values that are too close to each other
 * (based on @param minMzDelta) are identified and combined together.
 * The highest intensity associated with each respective (m/z, I) spectral peak is retained.
 * Note that this is a transitive agglomeration process, so if minMzDelta is too high,
 * many spectral peaks might get pulled in together.
 *
 * if @isMinMzDeltaPpm is true, @param minMzDelta is a ppm measurement.
 * If false, it is a Dalton measurement.
 *
 * @param minMzDelta
 */
void Fragment::agglomerateMzs(float minMzDelta, bool isMinMzDeltaPpm){

    if (minMzDelta <= 0) return;

    sortByMz();

    vector<vector<unsigned int>> agglomeratedMzs{};
    bool isHasMzAgglomeration = false;

    if (mzs.size() > 1) {

        //i == 0 case
        vector<unsigned int> currentGroup{0};

        for (unsigned int i = 1; i < mzs.size(); i++) {

            bool isAgglomerateMzs = false;
            if (isMinMzDeltaPpm) {
                isAgglomerateMzs = mzUtils::ppmDist(mzs.at(i), mzs.at(i-1)) < minMzDelta;
            } else {
                float mzDelta = mzs.at(i) - mzs.at(i-1);
                isAgglomerateMzs = mzDelta < minMzDelta;
            }

            if (!isAgglomerateMzs) {
                agglomeratedMzs.push_back(currentGroup);
                currentGroup = vector<unsigned int>{};
            } else {
                isHasMzAgglomeration = true;
            }
            currentGroup.push_back(i);
        }

        agglomeratedMzs.push_back(currentGroup);
    } else {
        for (unsigned int i = 0; i < mzs.size(); i++) {
            agglomeratedMzs.push_back(vector<unsigned int>{i});
        }
    }

    if (!isHasMzAgglomeration) return; //no need to agglomerate if no m/zs are too close to each other

    vector<float> updatedMzs{};
    vector<float> updatedIntensities{};
    vector<string> updatedFragmentLabels{};
    vector<int> updatedObsCount{};
    vector<vector<float>> updatedMedianIntensities{};

    for (unsigned int i = 0; i < agglomeratedMzs.size(); i++) {
        vector<unsigned int> aggMzs = agglomeratedMzs[i];

        unsigned int preferredPos = 0;

        float updatedMz = 0.0f;
        //float updatedIntensity = 0.0f;

        if (aggMzs.size() > 1) {
            float maxIntensity = -1.0f;
            for (unsigned int j = 0; j < aggMzs.size(); j++) {

                unsigned int posKey = aggMzs[j];

                updatedMz += mzs[posKey];

                float positionRepresentativeIntensity = 0.0f;
                if (consensusPositionToScanIntensities.find(static_cast<int>(posKey)) != consensusPositionToScanIntensities.end()) {
                    positionRepresentativeIntensity = *max_element(consensusPositionToScanIntensities[static_cast<int>(posKey)].begin(), consensusPositionToScanIntensities[static_cast<int>(posKey)].end());
                } else {
                    positionRepresentativeIntensity = intensity_array[posKey];
                }

                if (positionRepresentativeIntensity > maxIntensity) {
                    preferredPos = posKey;
                    maxIntensity = positionRepresentativeIntensity;
                }
            }

            updatedMz /= aggMzs.size();

        } else {

            updatedMz = mzs[aggMzs[0]];
            preferredPos = aggMzs[0];
        }

        if (!consensusPositionToScanIntensities.empty()) {
            updatedMedianIntensities.push_back(consensusPositionToScanIntensities[static_cast<int>(preferredPos)]);
        }

        updatedMzs.push_back(updatedMz);
        updatedIntensities.push_back(intensity_array[preferredPos]);
        updatedFragmentLabels.push_back(fragment_labels[preferredPos]);
        updatedObsCount.push_back(obscount[preferredPos]);
    }

    if (!updatedMedianIntensities.empty()) {
        consensusPositionToScanIntensities.clear();
        for (unsigned int i = 0; i < updatedMedianIntensities.size(); i++){
            consensusPositionToScanIntensities.insert(make_pair(i, updatedMedianIntensities[i]));
        }
    }

    mzs = updatedMzs;
    intensity_array = updatedIntensities;
    fragment_labels = updatedFragmentLabels;
    obscount = updatedObsCount;
}


/**
 * @brief Fragment::filterByMinIntensity
 * Irreversible transformation where all intensity values less than @param minIntensity
 * are dropped, and all relevant data structures are adjusted accordingly.
 *
 * @param minIntensity
 */
void Fragment::filterByMinIntensity(float minIntensity) {

    vector<float> filtered_mzs;
    vector<float> filtered_intensities;
    vector<string> filtered_fragment_labels;
    vector<int> filtered_obscount;
    vector<vector<float>> medianIntensities;

    for (unsigned int i= 0; i < nobs(); i++) {

        if (intensity_array[i] >= minIntensity) {

            filtered_mzs.push_back(mzs[i]);
            filtered_intensities.push_back(intensity_array[i]);
            filtered_fragment_labels.push_back(fragment_labels[i]);
            filtered_obscount.push_back(obscount[i]);

            if (!consensusPositionToScanIntensities.empty()) {
                medianIntensities.push_back(consensusPositionToScanIntensities[static_cast<int>(i)]);
            }

        }
    }

    if (!medianIntensities.empty()) {
        consensusPositionToScanIntensities.clear();
        for (unsigned int i = 0; i < medianIntensities.size(); i++){
            consensusPositionToScanIntensities.insert(make_pair(i, medianIntensities[i]));
        }
    }

    mzs = filtered_mzs;
    intensity_array = filtered_intensities;
    fragment_labels = filtered_fragment_labels;
    obscount = filtered_obscount;

}

/**
 * @brief Fragment::removeCoIsolations
 * Remove any m/z values that are within a 1 Da window centered around precursorMz.
 * Technically, this removes the original precursor along with any co-isolating species within tolerance.
 *
 * @param precursorMz
 * @param ppmTolr
 */
void Fragment::removeCoIsolations(float precursorMz, float ppmTolr){

    vector<float> filtered_mzs;
    vector<float> filtered_intensities;
    vector<string> filtered_fragment_labels;
    vector<int> filtered_obscount;
    vector<vector<float>> medianIntensities;

    float minPrecMz = precursorMz - 0.5f;
    float maxPrecMz = precursorMz + 0.5f;

    for (unsigned int i= 0; i < nobs(); i++) {

        if (mzs[i] < minPrecMz || mzs[i] > maxPrecMz || mzUtils::ppmDist(precursorMz, mzs[i]) <= ppmTolr) {

            filtered_mzs.push_back(mzs[i]);
            filtered_intensities.push_back(intensity_array[i]);
            filtered_fragment_labels.push_back(fragment_labels[i]);
            filtered_obscount.push_back(obscount[i]);

            if (!consensusPositionToScanIntensities.empty()) {
                medianIntensities.push_back(consensusPositionToScanIntensities[static_cast<int>(i)]);
            }

        }
    }

    if (!medianIntensities.empty()) {
        consensusPositionToScanIntensities.clear();
        for (unsigned int i = 0; i < medianIntensities.size(); i++){
            consensusPositionToScanIntensities.insert(make_pair(i, medianIntensities[i]));
        }
    }

    mzs = filtered_mzs;
    intensity_array = filtered_intensities;
    fragment_labels = filtered_fragment_labels;
    obscount = filtered_obscount;
}

double Fragment::normCosineScore(Fragment* a, Fragment* b, vector<int> ranks) {
    if (!a || !b) return -1.0;

    vector<float> aNorm = mzUtils::sumOfSquaresNorm(a->intensity_array);
    vector<float> bNorm = mzUtils::sumOfSquaresNorm(b->intensity_array);

    double score = 0.0;
    for (unsigned int i = 0; i < ranks.size(); i++) {
        int observedIndex = ranks[i];
        if (observedIndex != -1) {
            score += static_cast<double>(aNorm[i]) * static_cast<double>(bNorm[static_cast<unsigned int>(observedIndex)]);
        }
    }

    if (score > 1.0) {
        score = 1.0;
    } else if (score < 0) {
        score = 0;
    }
    return score;
}

/**
 * @brief Fragment::matchedPeakCosineScore
 * @param a
 * @param b
 * @param ranks
 *
 * Discard any unmatched peaks from a or b and compute a regular cosine score
 * considering only the matched peaks.
 *
 * Note that an unmatched peak is identical to matching to a peak with intensity 0.
 *
 * One might generalize this function by setting a minimum intensity threshold for
 * a matched peak intensity to be real, with a default value of 0.
 * @return
 */
double Fragment::matchedPeakCosineScore(Fragment* a, Fragment* b, vector<int> ranks) {
    if (!a || !b) return -1.0;

    vector<float> aMatchedIntensities{};
    vector<float> bMatchedIntensities{};

    for (unsigned int i = 0; i < ranks.size(); i++) {
        if (ranks[i] != -1) {
            aMatchedIntensities.push_back(a->intensity_array[i]);
            bMatchedIntensities.push_back(b->intensity_array[static_cast<unsigned long>(ranks[i])]);
        }
    }

    // Need at least two peaks for this score to make any sense.
    if (aMatchedIntensities.size() < 2) return 0;

    vector<float> aNorm = mzUtils::sumOfSquaresNorm(aMatchedIntensities);
    vector<float> bNorm = mzUtils::sumOfSquaresNorm(bMatchedIntensities);

    double score = 0.0;
    for (unsigned int i = 0; i < aNorm.size(); i++) {
        score += static_cast<double>(aNorm[i]) * static_cast<double>(bNorm[i]);
    }

    if (score > 1.0) {
        score = 1.0;
    } else if (score < 0) {
        score = 0;
    }
    return score;
}

/**
 * @brief encodeMsMsSpectrum
 * Convert masses and intensities from *fragment into string of form {{masses}{intensities}}
 * e.g.{{mass1,mass2,mass3},{intensity1,intensity2,intensity3}}
 * e.g.{{57.070015,58.0655441,60.0447273},{0.0146532,0.018739,0.0063873}}
 *
 * @param fragment: should have mzs and intensities
 * @param numDigits: number of digits after decimal point to include in string encoding.
 * Any fragments where intensity would have string representation of 0 is automatically
 * excluded from the encoded spectrum.
 *
 * @return encodedString
 */
string Fragment::encodeSpectrum(int numDigits, string type, bool isNormalizeToMaxIntensity) {

    sortedBy = Fragment::SortType::None;
    sortByMz();

    vector<float> mzs = this->mzs;
    auto c = intensity_array;
    vector<float> normalized_intensities(c.size());

    if (isNormalizeToMaxIntensity) {

        float m = *std::max_element(c.begin(), c.end());
        std::transform(c.begin(), c.end(), normalized_intensities.begin(), [&m](float intensity){return intensity /= m;});
        float minIntensity = static_cast<float>(pow(10, -numDigits));

        if (minIntensity > 0) {

            vector<float> mzCleaned;
            vector<float> intensityCleaned;

            for (unsigned int i = 0; i < mzs.size(); i++) {

                if (normalized_intensities[i] >= minIntensity) {
                    mzCleaned.push_back(mzs[i]);
                    intensityCleaned.push_back(normalized_intensities[i]);
                }
            }

            mzs = mzCleaned;
            normalized_intensities = intensityCleaned;
        }
    } else {
        normalized_intensities = c;
    }

    stringstream encodedSpectrum;

    encodedSpectrum << std::fixed << setprecision(numDigits);

    if (type == "encoded") {

        stringstream mzEncoding;
        stringstream intensityEncoding;

        mzEncoding << std::fixed << setprecision(numDigits);

        if (isNormalizeToMaxIntensity) {
            intensityEncoding << std::fixed << setprecision(numDigits);
        } else {
            intensityEncoding << std::fixed << setprecision(0);
        }

        for (unsigned int i = 0; i < mzs.size(); i++) {
            if (i > 0) {
                mzEncoding <<  ",";
                intensityEncoding << ",";
            }
            mzEncoding << mzs[i];
            intensityEncoding  << normalized_intensities[i];
        }

        encodedSpectrum << "{{" << mzEncoding.str() << "},{" << intensityEncoding.str() << "}}";

    } else {
        //fall back to tab-delimited list
        for (unsigned int i = 0; i < mzs.size(); i++) {
            if (isNormalizeToMaxIntensity) {
                encodedSpectrum << mzs[i] << "\t" << normalized_intensities[i] << "\n";
            } else {
                encodedSpectrum << std::fixed << setprecision(numDigits)
                                << mzs[i] << "\t"
                                << std::fixed << setprecision(0)
                                << normalized_intensities[i] << "\n";
            }
        }
    }

    return encodedSpectrum.str();
}

//Issue 585
void FragmentationMatchScore::addLabelSpecificMatches(string compoundLabel, bool debug) {

    bool isDiagnosticFragment = false;
    bool isAcylChainFragment = false;
    bool isSn1Fragment = false;
    bool isSn2Fragment = false;
    bool isSn3Fragment = false;
    bool isSn4Fragment = false;
    bool isOxidationFragment = false;

    vector<string> labelElements{};
    mzUtils::split(compoundLabel, '/', labelElements);

    //check key characters to determine fragment type

    for (string frag : labelElements) {
        if (debug) {
            cout << "frag: '" << frag;
        }

        if (frag.size() > 0) {
            char c = frag[0];
            if (debug) {
                cout << "', char: '" << c << "':" << endl;
            }

            if (c == '*'){
                isDiagnosticFragment = true;
            } else if (c == '@') {
                isSn1Fragment = true;
                isAcylChainFragment = true;
            } else if (c == '$'){
                isSn2Fragment = true;
                isAcylChainFragment = true;
            } else if (c == '!') {
                isSn3Fragment = true;
                isAcylChainFragment = true;
            } else if (c == '^') {
                isSn4Fragment = true;
                isAcylChainFragment = true;
            } else if (c == '&') {
                isOxidationFragment = true;
            }
        }
    }

    if (debug) {
        cout << "\tisDiagnosticFragment? " << (isDiagnosticFragment ? "yes" : "no") << "\n"
             << "\tisAcylChainFragment? " << (isAcylChainFragment ? "yes" : "no") << "\n"
             << "\tisSn1Fragment? " << (isSn1Fragment ? "yes" : "no") << "\n"
             << "\tisSn2Fragment? " << (isSn2Fragment ? "yes" : "no") << "\n"
             << "\tisSn3Fragment? " << (isSn3Fragment ? "yes" : "no") << "\n"
             << "\tisSn4Fragment? " << (isSn4Fragment ? "yes" : "no") << "\n"
             << "\tisOxidationFragment? " << (isOxidationFragment ? "yes" : "no") << "\n";
    }

    if (isDiagnosticFragment) numDiagnosticMatches++;
    if (isAcylChainFragment) numAcylChainMatches++;
    if (isSn1Fragment) numSn1Matches++;
    if (isSn2Fragment) numSn2Matches++;
    if (isSn3Fragment) numSn3Matches++;
    if (isSn4Fragment) numSn4Matches++;
    if (isOxidationFragment) numOxidations++;
}

Fragment* Fragment::createFromScans(vector<Scan*>& scans, shared_ptr<PeaksSearchParameters> params, bool debug) {

    Fragment *parent = nullptr;
    vector<Fragment*> fragments{};

    for (Scan* scan : scans) {

        if (debug) {
            cout << "[Fragment::createFromScans()]: Scan #" << scan->scannum << ": " << scan->nobs() << " peaks." << endl;

            cout << "[Fragment::createFromScans()]: Fragment Params: "
                 << "minFracIntensity=" << params->scanFilterMinFracIntensity << ", "
                 << "minSNRatio=" << params->scanFilterMinSNRatio << ", "
                 << "maxNumberOfFragments=" << params->scanFilterMaxNumberOfFragments << ", "
                 << "baseLinePercentile=" << params->scanFilterBaseLinePercentile << ", "
                 << "isRetainFragmentsAbovePrecursorMz=" << params->scanFilterIsRetainFragmentsAbovePrecursorMz << ", "
                 << "precursorPurityPpm=" << params->scanFilterPrecursorPurityPpm << ", "
                 << "minIntensity=" << params->scanFilterMinIntensity
                 << endl;
        }

        Fragment *frag = new Fragment(
            scan,
            params->scanFilterMinFracIntensity,
            params->scanFilterMinSNRatio,
            params->scanFilterMaxNumberOfFragments,
            params->scanFilterBaseLinePercentile,
            params->scanFilterIsRetainFragmentsAbovePrecursorMz,
            params->scanFilterPrecursorPurityPpm,
            params->scanFilterMinIntensity
            );

        if (debug) {
            cout << "[Fragment::createFromScans()]: Fragment: " << frag->nobs() << " peaks." << endl;
        }

        fragments.push_back(frag);
    }

    if (fragments.empty()) {
        parent = new Fragment();
        if (debug) cout << "[Fragment::createFromScans()]: No valid scans detected!." << endl;
    } else if (fragments.size() == 1) {
        parent = fragments.at(0);

        //Here, a 'consensus' built out of a single fragment just passes back the parent fragment data.
        parent->consensus = new Fragment();
        parent->consensus->mzs = parent->mzs;
        parent->consensus->intensity_array = parent->intensity_array;

        if (debug) cout << "[Fragment::createFromScans()]: A single-scan consensus spectrum was created with " << parent->nobs() << " peaks." << endl;
    } else {
        parent = fragments.at(0);
        for (unsigned int i = 1; i < fragments.size(); i++) {
            Fragment *frag = fragments[i];
            parent->addFragment(frag);
        }
        if (debug) {
            cout << "[Fragment::createFromScans()]: Successfully defined parent Fragment*, and added "
                 << (fragments.size()-1)
                 << " child Fragment*s."
                 <<  endl;
        }
        parent->buildConsensus(
            params->consensusPpmTolr,
            params->consensusIntensityAgglomerationType,
            params->consensusIsIntensityAvgByObserved,
            params->consensusIsNormalizeTo10K,
            params->consensusMinNumMs2Scans,
            params->consensusMinFractionMs2Scans,
            params->consensusIsRetainOriginalScanIntensities,
            mzUtils::decodeMzRemovedStr(params->consensusMs2MzRemovedStr),
            params->consensusMs2MzRemovedTol);

        if (debug) {
            cout << "[Fragment::createFromScans()]: Successfully built a true consensus spectrum." << endl;
        }

    }

    if (debug) {
        cout << "[Fragment::createFromScans()]: Final Fragment  # peaks=" << parent->consensus->nobs() << endl;
    }

    return parent;
}
