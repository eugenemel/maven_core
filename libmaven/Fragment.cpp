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
    this->sampleName = scan->sample->sampleName;
    this->scanNum = scan->scannum;
    this->precursorCharge = scan->precursorCharge;
    this->group = nullptr;
    this->consensus = nullptr;
	this->mergeCount=0;
	this->mergedScore=0;
	this->clusterId=0;
    this->scanNumMap={};
    scanNumMap.insert(make_pair(scan->sample, unordered_set<int>()));
    scanNumMap[scan->sample].insert(scan->scannum);

    //faster implementation when fewer filters specified
    if (minSNRatio <= 0 && (maxNumberOfFragments < 0 || maxNumberOfFragments >= scan->nobs()) && baseLinePercentile <= 0 && isRetainFragmentsAbovePrecursorMz){
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

        for (unsigned int j=0; j<mzarray.size() && j < maxNumberOfFragments; j++ ) {
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
    this->sampleName = scan->sample->sampleName;
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
    scanNumMap.insert(make_pair(scan->sample, unordered_set<int>()));
    scanNumMap[scan->sample].insert(scan->scannum);
}

//delete
Fragment::~Fragment() {
    mzUtils::delete_all(brothers);
    if(consensus) delete(consensus);
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

    //TODO: assess numDiagnosticFragments matched
    //TODO: respect intensity thresholds?

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
 * @return ranks vector<int>, of length of a.
 *
 * rank[a_position] = b_position
 * when sorted by m/z
 *
 * Note that these pairs are only valid if the two MS/MS spectra remain sorted by Mz.
 */
vector<int> Fragment::findFragPairsGreedyMz(Fragment* a, Fragment* b, float maxMzDiff) {

    //Sort spectra by m/z
    a->sortByMz();
    b->sortByMz();

    //Identify all valid possible fragment pairs (based on tolerance),
    //record with mzDelta

    vector<pair<float,pair<unsigned int, unsigned int>>> fragPairsWithMzDeltas;

    for (unsigned int i = 0; i < a->mzs.size(); i++){

        float mz_a = a->mzs.at(i);

        for (unsigned int j = 0; j < b->mzs.size(); j++) {

            float mz_b = b->mzs.at(j);

            if (mz_a - mz_b > maxMzDiff) {
                //Out of tolerance - could possibly come into tolerance as j increases.
            } else if (mz_b - mz_a > maxMzDiff) {
                //Out of tolerance - cannot possibly come into tolerance as j increases.
                break;
            } else {
                //In tolerance - record dissimilarity as candidate match.

                float mzDelta = abs(mz_a - mz_b);

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
    vector<int> ranks (a->mzs.size(),-1);
    set<unsigned int> claimedAFrags;
    set<unsigned int> claimedBFrags;

//    cerr << "match: ref <--> obs" << endl;

    for (pair<float,pair<unsigned int, unsigned int>> fragPairWithMzDelta : fragPairsWithMzDeltas){
        pair<unsigned int, unsigned int> fragPair = fragPairWithMzDelta.second;
        unsigned int a_frag = fragPair.first;
        unsigned int b_frag = fragPair.second;

        if (claimedAFrags.count(a_frag) == 0 && claimedBFrags.count(b_frag) == 0){

            matches.push_back(fragPair);
            ranks[a_frag] = static_cast<int>(b_frag);

//            cerr << "match: " << a->mzs.at(a_frag) << " <--> " << b->mzs.at(b_frag) << endl;

            claimedAFrags.insert(a_frag);
            claimedBFrags.insert(b_frag);
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
 */
void Fragment::buildConsensus(float productPpmTolr,
                              ConsensusIntensityAgglomerationType consensusIntensityAgglomerationType,
                              bool isIntensityAvgByObserved,
                              bool isNormalizeIntensityArray,
                              int minNumMs2ScansForConsensus,
                              float minFractionMs2ScansForConsensus) {

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
    if (consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean){
        for (unsigned int i = 0; i < Cons->intensity_array.size(); i++){
            posToIntensityMap.insert(make_pair(i, vector<float>{Cons->intensity_array[i]}));
        }
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
                if (consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean) posToIntensityMap[posA].push_back(intB);

            } else if ( posA == -1 ) {
                Cons->mzs.push_back(mzB);
                Cons->intensity_array.push_back(intB);
                Cons->obscount.push_back(1);
                Cons->fragment_labels.push_back(brother->fragment_labels[j]);

                //Issue 217
                if (consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean) posToIntensityMap.insert(make_pair(Cons->mzs.size()-1,vector<float>{intB}));
            }
        }

        //Issue 217: restructure map based on sorting
        if (consensusIntensityAgglomerationType != ConsensusIntensityAgglomerationType::Mean){

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

        map<mzSample*, unordered_set<int>> brotherMap = brother->scanNumMap;
        for (auto it = brotherMap.begin(); it != brotherMap.end(); ++it) {

            mzSample* sample = it->first;
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
        mzSample* sample = it->first;
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

    //compute avererage purity
    Cons->purity = this->consensusPurity();

    //optionally filter values from scans based on number of observations
    if (minNumMs2ScansForConsensus > 0 || minFractionMs2ScansForConsensus > 0) {

        vector<float> filtered_mzs;
        vector<float> filtered_intensities;
        vector<string> filtered_fragment_labels;
        vector<int> filtered_obscount;

        for (unsigned int i = 0; i < Cons->mzs.size(); i++){
            float frac = static_cast<float>(Cons->obscount[i]) / static_cast<float>(N);
            if (Cons->obscount[i] >= minNumMs2ScansForConsensus && frac >= minFractionMs2ScansForConsensus) {
                filtered_mzs.push_back(Cons->mzs[i]);
                filtered_intensities.push_back(Cons->intensity_array[i]);
                filtered_fragment_labels.push_back(Cons->fragment_labels[i]);
                filtered_obscount.push_back(Cons->obscount[i]);
            }
        }

        Cons->mzs = filtered_mzs;
        Cons->intensity_array = filtered_intensities;
        Cons->fragment_labels = filtered_fragment_labels;
        Cons->obscount = filtered_obscount;
    }

    //agglomerate intensity values
	if( brothers.size() >= 1) {
        if (consensusIntensityAgglomerationType == Mean) {
            for( unsigned int i=0; i < Cons->intensity_array.size(); i++){
                Cons->intensity_array[i] /= (isIntensityAvgByObserved ? Cons->obscount[i] : N);
            }
        } else if (consensusIntensityAgglomerationType == Median) {
            for( unsigned int i=0; i < Cons->intensity_array.size(); i++){

                vector<float> intensityValues = posToIntensityMap[static_cast<int>(i)];

                //Issue 227 debugging
                if (Cons->mzs[i] > 124 && Cons->mzs[i] < 124.1) {
                    cerr << "i=" << i << ", mz=" << Cons->mzs[i] << ", intensity vec=" << intensityValues.size() << endl;
                }

                //Issue 227
                if (isIntensityAvgByObserved) {
                    intensityValues.erase(remove_if(
                                              intensityValues.begin(), intensityValues.end(),
                                          [](const float& intensity){
                                            return intensity <= 0;
                                        }), intensityValues.end()
                                );

                }
                Cons->intensity_array[i] = mzUtils::median(intensityValues);
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
    if(sortedBy == Fragment::SortType::Intensity) return; //sppedup already sorted

    vector<unsigned int>order = intensityOrderDesc();
    vector<float> a(mzs.size());
    vector<float> b(intensity_array.size());
    vector<int> c(obscount.size());
    vector<string> fragLabels(fragment_labels.size());

    map<int,string>d;
    for(unsigned int i=0; i<order.size(); i++) {
        b[i] = intensity_array[order[i]];
        a[i] = mzs[order[i]];
        fragLabels[i] = fragment_labels[order[i]];
        if(order[i] < obscount.size()) c[i] = obscount[order[i]];
        if(annotations.count(order[i])>0) d[i] = annotations[order[i]];

    };
    mzs = a;
    intensity_array = b;
    obscount = c;
    annotations=d;
    fragment_labels = fragLabels;
    sortedBy = Fragment::SortType::Intensity;
}	

void Fragment::sortByMz() {

    if(sortedBy == Fragment::SortType::Mz) return; //speedup already sorted

    vector<int>order = mzOrderInc();
    vector<float> a(mzs.size());
    vector<float> b(intensity_array.size());
    vector<int> c(obscount.size());
    vector<string> fragLabels(fragment_labels.size());

    map<int,string>d;

    for(unsigned int i=0; i<order.size(); i++) {
        b[i] = intensity_array[order[i]];
        a[i] = mzs[order[i]];
        fragLabels[i] = fragment_labels[order[i]];
        if(order[i] < obscount.size()) c[i] = obscount[order[i]];
        if (annotations.count(order[i])>0) d[i] = annotations[order[i]];

    };

    mzs = a;
    intensity_array = b;
    obscount = c;
    annotations=d;
    fragment_labels = fragLabels;
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

