#include "Fragment.h"
#include "mzSample.h"

//empty constructor
Fragment::Fragment() { 
	precursorMz = 0; 
	polarity=1; 
	scanNum=0; 
	rt=0;
   	collisionEnergy=0; 
	consensus=NULL; 
	precursorCharge=0; 
	isDecoy=false; 
	sortedBy=None; 
	group=NULL;
    mergeCount=0;
    purity=0;
    consensus=NULL;
    group=NULL;
    clusterId=0;
	mergedScore=0;
}


//build fragment based on MS2 scan
Fragment::Fragment(Scan* scan, float minFractionalIntensity, float minSigNoiseRatio,unsigned int maxFragmentSize) {
    this->precursorMz = scan->precursorMz;
    this->collisionEnergy = scan->collisionEnergy;
    this->polarity = scan->getPolarity();
    this->sampleName = scan->sample->sampleName;
    this->scanNum = scan->scannum;
    this->precursorCharge = scan->precursorCharge;
    this->sortedBy= Fragment::SortType::Mz;
	this->mergeCount=0;
	this->mergedScore=0;
	this->clusterId=0;

    //TODO: make this configurable? as a scoring parameter?
    int baseLineLevel=5; //lowest 5% of data are considered to be baseline

    //don't worry about baseline.. keeping all points
	/*
    if(scan->nobs()<maxFragmentSize) { 
        minSigNoiseRatio=0; 
        minFractionalIntensity=0;
    }
	*/

    vector<pair<float,float> >mzarray = scan->getTopPeaks(minFractionalIntensity,minSigNoiseRatio,baseLineLevel);

    for(unsigned int j=0; j<mzarray.size() && j < maxFragmentSize; j++ ) {
		if (mzarray[j].second < this->precursorMz-1 ) { //remove fragments higher than precursorMz
			this->mzs.push_back(mzarray[j].second);
			this->intensity_array.push_back(mzarray[j].first);
		}
    }
    this->obscount = vector<int>( this->mzs.size(), 1);
    this->group = NULL;
    this->consensus =NULL;
    this->rt = scan->rt;
    this->purity = scan->getPrecursorPurity(10.00);  //this might be slow
    this->sortByMz();
}


//delete
Fragment::~Fragment() {
    mzUtils::delete_all(brothers);
    if(consensus != NULL) delete(consensus);
}

//make a copy of Fragment.
Fragment::Fragment( Fragment* other) { 
    this->precursorMz = other->precursorMz;
    this->polarity = other->polarity;
    this->mzs = other->mzs;
    this->intensity_array = other->intensity_array;
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

//WARNING: caller Compound::scoreCompoundHit calls Fragment::sortByIntensity() just before Fragment::scoreMatch().
//This seems dangerous.  If this sorting is required for scoring, it should be done inside this method,
//not just before calling.

FragmentationMatchScore Fragment::scoreMatch(Fragment* other, float productPpmTolr) {
    FragmentationMatchScore s;
    if (mzs.size() < 2 or other->mzs.size() < 2) return s;

    //which one is smaller;
    Fragment* a = this;
    Fragment* b = other;

    s.ppmError = abs((a->precursorMz-b->precursorMz)/a->precursorMz*1e6);

    float maxDeltaMz = (productPpmTolr * static_cast<float>(a->precursorMz))/ 1000000;

    /*
     * ranks[x] = y
     * x = index of frag peak in a
     * y = index of frag peak in b
     */
    vector<int> ranks = findFragPairsGreedyMz(a, b, maxDeltaMz);
    //vector<int>ranks = compareRanks(a,b,productPpmTolr);
    //vector<int>ranks = locatePositions(a,b,productPpmTolr);

    for(int rank: ranks) { if(rank != -1) s.numMatches++; }

    //annotate?
    for(int i=0; i<ranks.size();i++) other->annotations[ranks[i]]=this->annotations[i];

    s.fractionMatched = s.numMatches / a->nobs();
    s.spearmanRankCorrelation = spearmanRankCorrelation(ranks);
    s.ticMatched = ticMatched(ranks);
    s.mzFragError =  mzErr(ranks,b);
    s.dotProduct = dotProduct(b);
    s.hypergeomScore  = SHP(s.numMatches,a->nobs(),b->nobs(),100000) + s.ticMatched; // ticMatch is tie breaker
    s.mvhScore = MVH(ranks,b);
    s.weightedDotProduct = mzWeightedDotProduct(ranks,b);
    s.matchedQuantiles = matchedRankVector(ranks,b);
    //s.dotProductShuffle = this->dotProductShuffle(b,2000);

    //cerr << "scoreMatch:\n" << a->nobs() << "\t" << b->nobs() << "\t" << s.numMatches << " hyper=" << s.hypergeomScore << "\n";

    return s;
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

vector<int> Fragment::compareRanks(Fragment* a, Fragment* b, float productPpmTolr) {
    bool verbose=false;
    vector<int> ranks (a->mzs.size(),-1);	//missing value == -1
    for(unsigned int i=0; i<a->mzs.size(); i++ ) {
        for( unsigned int j=0; j < b->mzs.size(); j++ ) {
            if (mzUtils::ppmDist(a->mzs[i],b->mzs[j])<productPpmTolr) { ranks[i] = j; break; }   //this needs optimization..
        }
    }
    if (verbose) {
        cerr << " compareranks: " << a->sampleName << endl;
        for(unsigned int i=0; i < ranks.size(); i++ ) {
            float mz2=0;
            float ints2=0;

            if(ranks[i]>=0) { mz2=b->mzs[ ranks[i] ]; ints2= b->intensity_array[ ranks[i] ]; }
                std::cerr << ranks[i] << "," << a->mzs[i] << "\t" << mz2 << "\t\t" <<  a->intensity_array[i] << "\t" << ints2 << endl;
            }

        }

        return ranks;
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

    vector<pair<float,pair<uint, uint>>> fragPairsWithMzDeltas;

    for (uint i = 0; i < a->mzs.size(); i++){

        float mz_a = a->mzs.at(i);

        for (uint j = 0; j < b->mzs.size(); j++) {

            float mz_b = b->mzs.at(j);

            if (mz_a - mz_b > maxMzDiff) {
                //Out of tolerance - could possibly come into tolerance as j increases.
            } else if (mz_b - mz_a > maxMzDiff) {
                //Out of tolerance - cannot possibly come into tolerance as j increases.
                break;
            } else {
                //In tolerance - record dissimilarity as candidate match.

                float mzDelta = abs(mz_a - mz_b);

                pair<uint,uint> peakPair (i, j); //First position is reserved for a, second for b

                pair<float,pair<uint, uint>> fragPairWithMzDelta (mzDelta, peakPair);

                fragPairsWithMzDeltas.push_back(fragPairWithMzDelta);
            }

        }
    }

    //sort pairs in increasing order by mzDelta
    std::sort(fragPairsWithMzDeltas.begin(), fragPairsWithMzDeltas.end(),
              [ ](const pair<float,pair<uint, uint>>& lhs, const pair<float,pair<uint, uint>>& rhs){
        if (lhs.first < rhs.first) {
            return -1;
        } else if (lhs.first > rhs.first) {
            return 1;
        } else {
            if (lhs.second.first < rhs.second.first) {
                return -1;
            } else if (lhs.second.first > rhs.second.first){
                return 1;
            } else if (lhs.second.second < rhs.second.second){
                return -1;
            } else if (lhs.second.second > rhs.second.second){
                return 1;
            } else {
                return 0;
            }
        }
    });

    //Once a fragment has been claimed in a frag pair, it may not be involved in any other
    //frag pair.
    vector<pair<uint, uint>> matches;
    vector<int> ranks (a->mzs.size(),-1);
    set<uint> claimedAFrags;
    set<uint> claimedBFrags;

    for (pair<float,pair<uint,uint>> fragPairWithMzDelta : fragPairsWithMzDeltas){
        pair<uint, uint> fragPair = fragPairWithMzDelta.second;
        uint a_frag = fragPair.first;
        uint b_frag = fragPair.second;

        if (claimedAFrags.count(a_frag) == 0 && claimedBFrags.count(b_frag) == 0){

            matches.push_back(fragPair);
            ranks[a_frag] = static_cast<int>(b_frag);

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
	this->obscount.resize(n);
}


void Fragment::buildConsensus(float productPpmTolr) {
    if(this->consensus != NULL) {  delete(this->consensus); this->consensus=NULL; }

	//find brother with largest nobs
	Fragment* seed= this;
    for(Fragment* b: brothers) if (b->nobs() > seed->nobs()) seed = b;

    Fragment* Cons = new Fragment(seed);  //make a copy of self
    this->consensus = Cons;
    Cons->sortByMz();

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
            } else if ( posA == -1 ) {
                Cons->mzs.push_back(mzB);
                Cons->intensity_array.push_back(intB);
                Cons->obscount.push_back(1);
            }
        }
        Cons->sortByMz();
    }

    //compute retention time window
    Cons->rt  = this->consensusRt();

    //compute avererage purity
    Cons->purity = this->consensusPurity();

    //average values 
	if( brothers.size() >= 1) {
		int N = 1+brothers.size();
		for(unsigned int i=0; i<Cons->intensity_array.size(); i++) Cons->intensity_array[i] /= N;
	}
	
    Cons->sortByIntensity();
	if (Cons->intensity_array.size() > 1) {
		float maxValue = Cons->intensity_array[0];
		for(unsigned int i=0; i<Cons->intensity_array.size(); i++)  {
		   	Cons->intensity_array[i] = Cons->intensity_array[i]/maxValue*10000;
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

    //return positions in order from highest to lowest intenisty
    for(unsigned int i=0; i < _pairsarray.size(); i++) { position[i] = _pairsarray[i].second; }
    return position;
}


void Fragment::sortByIntensity() { 
    if(sortedBy == Fragment::SortType::Intensity) return; //sppedup already sorted

    vector<unsigned int>order = intensityOrderDesc();
    vector<float> a(mzs.size());
    vector<float> b(intensity_array.size());
    vector<int> c(obscount.size());
    map<int,string>d;
    for(unsigned int i=0; i<order.size(); i++) {
        b[i] = intensity_array[order[i]];
        a[i] = mzs[order[i]];
        if(order[i] < obscount.size()) c[i] = obscount[order[i]];
        if(annotations.count(order[i])>0) d[i] = annotations[order[i]];

    };
    mzs = a;
    intensity_array = b;
    obscount = c;
    annotations=d;
    sortedBy = Fragment::SortType::Intensity;
}	

void Fragment::sortByMz() {

    if(sortedBy == Fragment::SortType::Mz) return; //sppedup already sorted

    vector<int>order = mzOrderInc();
    vector<float> a(mzs.size());
    vector<float> b(intensity_array.size());
    vector<int> c(obscount.size());
    map<int,string>d;

    for(unsigned int i=0; i<order.size(); i++) {
        b[i] = intensity_array[order[i]];
        a[i] = mzs[order[i]];
        if(order[i] < obscount.size()) c[i] = obscount[order[i]];
        if (annotations.count(order[i])>0) d[i] = annotations[order[i]];

    };

    mzs = a;
    intensity_array = b;
    obscount = c;
    annotations=d;
    sortedBy = Fragment::SortType::Mz;
}

void Fragment::buildConsensusAvg() { 
    map<float,double> mz_intensity_map;
    map<float,double> mz_bin_map;
    map<float,int> mz_count;

    vector<Fragment*> fragmentList = brothers;
    fragmentList.push_back(this);

    for(unsigned int i=0; i<fragmentList.size(); i++) {
        Fragment* brother = fragmentList[i];
        for(unsigned int j=0; j < brother->mzs.size(); j++) {
            float bin = (round(brother->mzs[j]+0.5)*10)/10;
            mz_intensity_map[bin] += ((double) brother->intensity_array[j]);
            mz_bin_map[bin] += ((double)(brother->intensity_array[j])*(brother->mzs[j]));
            mz_count[bin]++;
        }
    }
    map<float,double>::iterator itr;
    for(itr = mz_intensity_map.begin(); itr != mz_intensity_map.end(); ++itr ) {
        float bin = (*itr).first;
        double totalIntensity=(*itr).second;
        double avgMz =  mz_bin_map[bin] / totalIntensity;
        cerr << "\t" << setprecision(3) << avgMz << " " << totalIntensity/mz_count[bin] << endl;
        //avgScan->mz.push_back((float)avgMz);
        //avgScan->intensity.push_back((float) totalIntensity / mz_count[bin]);
    }
}


double Fragment::spearmanRankCorrelation(const vector<int>& X) {
    double d2=0; 
    int N = X.size();
    int n=0;
    for(int i=0; i<N;i++ ) {	
        if (X[i] == -1 ) { //mising values set to average distance
            d2 += 2*i;
            n++;
        } else  {
            n++;
            d2 += (i-X[i])*(i-X[i]);
            //cerr << "n=" << n << " cor=" << d2 << endl;
        }
    }
    //double p = 1.00-(6.0*d2)/(n*((n*n)-1));
    return 1.00-(6.0*d2)/(n*((n*n)-1));
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
    double ERR=1000;
    for(unsigned int i=0; i<X.size(); i++ ) if (X[i] != -1) ERR = POW2( mzs[i] - other->mzs[ X[i] ]);
    return sqrt(ERR);
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
            obscount.push_back(1);
        }
    }
    sortByIntensity();
}

