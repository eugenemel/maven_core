#include "directinfusionprocessor.h"
#include "lipidsummarizationutils.h"

#include <chrono>

using namespace std;
using namespace mzUtils;

shared_ptr<DirectInfusionSearchSet> DirectInfusionProcessor::getSearchSet(mzSample* sample,
                                                                              const vector<Compound*>& compounds,
                                                                              const vector<Adduct*>& adducts,
                                                                              shared_ptr<DirectInfusionSearchParameters> params,
                                                                              bool debug) {

    shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet = shared_ptr<DirectInfusionSearchSet>(new DirectInfusionSearchSet());

    for (Scan* scan : sample->scans){
        if (scan->mslevel == 2){
            int mapKey = static_cast<int>(round(scan->precursorMz+0.001f)); //round to nearest int

            if (directInfusionSearchSet->mzRangesByMapKey.find(mapKey) == directInfusionSearchSet->mzRangesByMapKey.end()) {
                float precMzMin = scan->getPrecMzMin();
                float precMzMax = scan->getPrecMzMax();

                directInfusionSearchSet->mzRangesByMapKey.insert(make_pair(mapKey, make_pair(precMzMin, precMzMax)));
            }

            directInfusionSearchSet->mapKeys.insert(mapKey);
        }
    }

    typedef map<int, pair<float, float>>::iterator mzRangeIterator;

    MassCalculator massCalc;

    if (debug) cerr << "Organizing database into map for fast lookup..." << endl;

    for (Compound *compound : compounds) {
        for (Adduct *adduct : adducts) {

            if (SIGN(adduct->charge) != SIGN(compound->charge)) {
                continue;
            }

            float compoundMz = compound->precursorMz;

            if (params->ms1IsRequireAdductPrecursorMatch){

                if (compound->adductString != adduct->name){
                    continue;
                }

                //TODO: this code works for compounds that do not have an 'adductString' set, but the adduct in the name.
                //delete this eventually.
//                if(compound->name.length() < adduct->name.length() ||
//                   compound->name.compare (compound->name.length() - adduct->name.length(), adduct->name.length(), adduct->name) != 0){
//                    continue;
//                }
            } else {
                compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));
            }

            //determine which map key to associate this compound, adduct with

            for (mzRangeIterator it = directInfusionSearchSet->mzRangesByMapKey.begin(); it != directInfusionSearchSet->mzRangesByMapKey.end(); ++it) {
                int mapKey = it->first;
                pair<float, float> mzRange = it->second;

                if (compoundMz > mzRange.first && compoundMz < mzRange.second) {

                    if (directInfusionSearchSet->compoundsByMapKey.find(mapKey) == directInfusionSearchSet->compoundsByMapKey.end()){
                        directInfusionSearchSet->compoundsByMapKey.insert(make_pair(mapKey, vector<pair<Compound*, Adduct*>>()));
                    }

                    directInfusionSearchSet->compoundsByMapKey[mapKey].push_back(make_pair(compound, adduct));
                    break;
                }
            }
        }
    }

    return directInfusionSearchSet;

}

vector<Ms3Compound*> DirectInfusionProcessor::getMs3CompoundSet(const vector<Compound*>& compounds,
                                                                bool debug){
    vector<Ms3Compound*> ms3Compounds(compounds.size());

    for (unsigned int i = 0; i < compounds.size(); i++) {
        ms3Compounds[i] = new Ms3Compound(compounds[i]); //WARNING: this delete this at some point to avoid memory leak
    }

    if (debug) cout << "Created database of " << ms3Compounds.size() << " Ms3Compounds." << endl;

    return ms3Compounds;
}


vector<vector<tuple<double, double, Scan*>>> DirectInfusionProcessor::organizeMs3ScansByPrecursor(
        vector<tuple<double, double, Scan*>> allMs3Scans,
        double ms3AnalysisMs1PrecursorPpmTolr,
        double ms3PrecursorPpmTolr,
        bool debug){

    sort(allMs3Scans.begin(), allMs3Scans.end(), [](const tuple<double, double, Scan*>& lhs, const tuple<double, double, Scan*>& rhs){

        double lhsMs1PrecMz = get<0>(lhs);
        double rhsMs1PrecMz = get<0>(rhs);

        double lhsMs2PrecMz = get<1>(lhs);
        double rhsMs2PrecMz = get<1>(rhs);

        int lhsScanNum = get<2>(lhs)->scannum;
        int rhsScanNum = get<2>(rhs)->scannum;

        if (abs(lhsMs1PrecMz - rhsMs1PrecMz) < 1e-7) {
            if (abs(lhsMs2PrecMz - rhsMs2PrecMz) < 1e-7) {
                return lhsScanNum < rhsScanNum;
            } else {
                return lhsMs2PrecMz < rhsMs2PrecMz;
            }
        } else {
            return lhsMs1PrecMz < rhsMs1PrecMz;
        }
    });

    vector<vector<tuple<double, double, Scan*>>> ms3ScanGroups;
    vector<tuple<double,double,Scan*>> scanGroup;

    pair<double,double> lastPrecMz = make_pair(0.0, 0.0);
    unsigned int numProcessedPairs = 0;

    for (unsigned int i = 0; i < allMs3Scans.size(); i++) {

        scanGroup = vector<tuple<double,double,Scan*>>();

        tuple<double, double, Scan*> ithScanPair = allMs3Scans[i];

        lastPrecMz = make_pair(get<0>(ithScanPair), get<1>(ithScanPair));

        scanGroup.push_back(ithScanPair);

        for (unsigned int j = i+1; j < allMs3Scans.size(); j++) {

            tuple<double, double, Scan*> jthScanPair = allMs3Scans[j];

            bool isMatchingMs1PrecursorTolr = mzUtils::ppmDist(get<0>(ithScanPair), get<0>(jthScanPair)) <= ms3AnalysisMs1PrecursorPpmTolr;
            bool isMatchingMs2PrecursorTolr = mzUtils::ppmDist(get<1>(ithScanPair), get<1>(jthScanPair)) <= ms3PrecursorPpmTolr;

            if (isMatchingMs1PrecursorTolr && isMatchingMs2PrecursorTolr) {
                scanGroup.push_back(jthScanPair);
                lastPrecMz = make_pair(get<0>(jthScanPair), get<1>(jthScanPair));

                if (j == allMs3Scans.size()-1) {
                    i = static_cast<unsigned int>(allMs3Scans.size()); // avoid outer loop
                    break;
                }
            } else {

                i = j;
                ms3ScanGroups.push_back(scanGroup);

                numProcessedPairs += scanGroup.size();

                if (debug) cout << "i=" << i << ", numProcessedPairs= " << numProcessedPairs << endl;

                i--; // necessary b/c outer for loop will increment i
                break;
            }

//            if (mzUtils::ppmDist(jthScanPair.first, ithScanPair.first) > ms3PrecursorPpmTolr) {
//                i = j;
//                ms3ScanGroups.push_back(scanGroup);

//                numProcessedPairs += scanGroup.size();

//                if (debug) cout << "i=" << i << ", numProcessedPairs= " << numProcessedPairs << endl;

//                i--; // necessary b/c outer for loop will increment i
//                break;
//            } else {
//                scanGroup.push_back(jthScanPair);
//                lastPrecMz = jthScanPair.first;

//                if (j == allMs3Scans.size()-1) {
//                    i = static_cast<unsigned int>(allMs3Scans.size()); //avoid outer loop
//                    break;
//                }
//            }
        }
    }

    if (!scanGroup.empty()){
        ms3ScanGroups.push_back(scanGroup);
        numProcessedPairs += scanGroup.size();
        if (debug) cout << "i=" << allMs3Scans.size() << ", numProcessedPairs="<< numProcessedPairs << endl;
    }

    //debugging
    if (debug) {

        unsigned int spCounter = 0;
        unsigned int grpCounter = 0;
        for (auto sg : ms3ScanGroups) {

            grpCounter++;
            cout<< "group #" << grpCounter << ": scans=";
            for (auto sp : sg) {
                cout << get<2>(sp)->scannum << ", ";
            }
            cout << endl;

            spCounter += sg.size();
        }

        cout << "ms3 scans: # all=" << allMs3Scans.size() << ", # grouped=" << spCounter << endl;
    }

    return ms3ScanGroups;
}

vector<Ms3SingleSampleMatch*> DirectInfusionProcessor::processSingleMs3Sample(mzSample* sample,
                                                                                  const vector<Ms3Compound*>& ms3Compounds,
                                                                                  shared_ptr<DirectInfusionSearchParameters> params,
                                                                                  bool debug){

    //initialize output
    vector<Ms3SingleSampleMatch*> output;

    map<int, vector<Scan>> ms3ScansByMzPrecursor{};

    //          <ms1Pre, ms2Pre>
    vector<tuple<double, double, Scan*>> allMs3Scans;

    vector<Scan*> validMs1Scans;

    for (Scan* scan : sample->scans) {

        if (scan->mslevel == 3 &&
                (params->scanFilterMs3MinRt <= -1.0f || scan->rt >= params->scanFilterMs3MinRt) &&
                (params->scanFilterMs3MaxRt <= -1.0f || scan->rt <= params->scanFilterMs3MaxRt)) {

            tuple<double, double, Scan*> ms3ScanData = tuple<double, double, Scan*>(scan->ms1PrecursorForMs3, scan->precursorMz, scan);
            allMs3Scans.push_back(ms3ScanData);

        } else if (scan->mslevel == 1 &&
                              scan->filterString.find(params->ms1ScanFilter) != string::npos &&
                              (params->scanFilterMs1MinRt <= -1.0f || scan->rt >= params->scanFilterMs1MinRt) &&
                              (params->scanFilterMs1MaxRt <= -1.0f || scan->rt <= params->scanFilterMs1MaxRt)) {

            validMs1Scans.push_back(scan);
            if (debug) cout << "Added valid MS1 scan: " << scan->scannum << " " << scan->filterString << endl;
        }
    }

    if (debug) cout << "Computing consensus MS1 scan from " << validMs1Scans.size() << " MS1 scans..." << endl;

    Fragment *ms1Fragment = nullptr;
    for (auto & scan: validMs1Scans) {
        if (!ms1Fragment) {
            ms1Fragment = new Fragment(scan,
                                       params->scanFilterMinFracIntensity,
                                       params->scanFilterMinSNRatio,
                                       params->scanFilterMaxNumberOfFragments,
                                       params->scanFilterBaseLinePercentile,
                                       params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                       params->scanFilterPrecursorPurityPpm,
                                       params->scanFilterMinIntensity);
        } else {
            Fragment *ms1Brother = new Fragment(scan,
                                             params->scanFilterMinFracIntensity,
                                             params->scanFilterMinSNRatio,
                                             params->scanFilterMaxNumberOfFragments,
                                             params->scanFilterBaseLinePercentile,
                                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                             params->scanFilterPrecursorPurityPpm,
                                             params->scanFilterMinIntensity);

            ms1Fragment->addFragment(ms1Brother);
        }
    }

    if (ms1Fragment){
        ms1Fragment->buildConsensus(params->consensusMs1PpmTolr,
                                    params->consensusIntensityAgglomerationType,
                                    params->consensusIsIntensityAvgByObserved,
                                    params->consensusIsNormalizeTo10K,
                                    params->consensusMinNumMs1Scans,
                                    params->consensusMinFractionMs1Scans
                                    );

        ms1Fragment->consensus->sortByMz();
    }

    if (debug) cout << "Finished computing consensus MS1 scan." << endl;
    if (debug){
        cout << "ms1Fragment->buildConsensus() parameters:" << endl;
        cout << "\tconsensusMs1PpmTolr: " << to_string(params->consensusMs1PpmTolr) << endl;
        cout << "\tconsensusMinNumMs1Scans: " << to_string(params->consensusMinNumMs1Scans) << endl;
        cout << "\tconsensusMinFractionMs1Scans: " <<to_string(params->consensusMinFractionMs1Scans) << endl;
    }

    if (debug) {
        cout << "scans in consensus MS1 scan:" << endl;
        for (auto it = ms1Fragment->consensus->scanNumMap.begin(); it != ms1Fragment->consensus->scanNumMap.end(); ++it){
            cout << it->first->sampleName << ": " << endl;
            for (auto x : it->second) {
                cout << x << " ";
            }
            cout << endl;
        }
    }

    //                  ms1Pre, ms2Pre
    vector<vector<tuple<double, double, Scan*>>> ms3ScanGroups = DirectInfusionProcessor::organizeMs3ScansByPrecursor(
                allMs3Scans,
                static_cast<double>(params->ms3AnalysisMs1PrecursorPpmTolr),            //refers to ms1 target m/z
                static_cast<double>(params->ms3PrecursorPpmTolr),   //refers to ms2 target m/z
                debug);

    vector<tuple<double, double, Fragment*>> consensusMs3Spectra(ms3ScanGroups.size());

    for (unsigned int i = 0; i < ms3ScanGroups.size(); i++) {

        auto pairVector = ms3ScanGroups[i];
        Fragment *f = nullptr;

        double avgMs1PrecMz = 0;
        double avgMs2PrecMz = 0;

        for (auto pair : pairVector) {
            avgMs1PrecMz += get<0>(pair);
            avgMs2PrecMz += get<1>(pair);
            if (!f) {
                f = new Fragment(get<2>(pair),
                                 params->scanFilterMinIntensity,
                                 params->scanFilterMinSNRatio,
                                 params->scanFilterMaxNumberOfFragments,
                                 params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                 params->scanFilterPrecursorPurityPpm,
                                 params->scanFilterMinIntensity);
            } else {
                Fragment *brother = new Fragment(get<2>(pair),
                                                 params->scanFilterMinIntensity,
                                                 params->scanFilterMinSNRatio,
                                                 params->scanFilterMaxNumberOfFragments,
                                                 params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                                 params->scanFilterPrecursorPurityPpm,
                                                 params->scanFilterMinIntensity);
                f->addFragment(brother);
            }
        }

        avgMs1PrecMz /= pairVector.size();
        avgMs2PrecMz /= pairVector.size();

        f->buildConsensus(params->consensusMs3PpmTolr,
                          params->consensusIntensityAgglomerationType,
                          params->consensusIsIntensityAvgByObserved,
                          params->consensusIsNormalizeTo10K,
                          params->consensusMinNumMs3Scans,
                          params->consensusMinFractionMs3Scans
                          );

        f->consensus->sortByMz();
        consensusMs3Spectra[i] = tuple<double, double, Fragment*>(avgMs1PrecMz, avgMs2PrecMz, f);
    }

    sort(consensusMs3Spectra.begin(), consensusMs3Spectra.end(), [](const tuple<double, double, Fragment*>& lhs, const tuple<double, double, Fragment*>& rhs){
        double lhsMs1PrecMz = get<0>(lhs);
        double rhsMs1PrecMz = get<0>(rhs);

        double lhsMs2PrecMz = get<1>(lhs);
        double rhsMs2PrecMz = get<1>(rhs);

        int lhsScanNum = get<2>(lhs)->scanNum;
        int rhsScanNum = get<2>(rhs)->scanNum;

        if (abs(lhsMs1PrecMz - rhsMs1PrecMz) < 1e-7) {
            if (abs(lhsMs2PrecMz - rhsMs2PrecMz) < 1e-7) {
                return lhsScanNum < rhsScanNum;
            } else {
                return lhsMs2PrecMz < rhsMs2PrecMz;
            }
        } else {
            return lhsMs1PrecMz < rhsMs1PrecMz;
        }
    });

    unsigned int compoundCounter = 0;

    for (auto ms3Compound : ms3Compounds) {

        int numMs3Matches = 0;
        map<int, pair<Fragment*, vector<int>>> matchData{};
        string matchInfoDebugString("");

        double ms1PrecMz = ms3Compound->baseCompound->precursorMz;

        double ms1MinMz = ms1PrecMz - ms1PrecMz * params->ms3AnalysisMs1PrecursorPpmTolr/1000000.0;
        double ms1MaxMz = ms1PrecMz + ms1PrecMz * params->ms3AnalysisMs1PrecursorPpmTolr/1000000.0;

        auto lb = lower_bound(consensusMs3Spectra.begin(), consensusMs3Spectra.end(), ms1MinMz, [](const tuple<double, double, Fragment*>& lhs, const double& rhs){
            return get<0>(lhs) < rhs;
        });

        for (auto it = ms3Compound->ms3_fragment_mzs.begin(); it != ms3Compound->ms3_fragment_mzs.end(); ++it){

            double ms2PrecMz = mzUtils::intKeyToMz(it->first);

            for (unsigned int pos = lb - consensusMs3Spectra.begin(); pos < consensusMs3Spectra.size(); pos++) {

                tuple<double, double, Fragment*> data = consensusMs3Spectra[pos];

                double targetMs1PrecMz = get<0>(data);
                double targetMs2PrecMz = get<1>(data);

                bool isMatchingMs1PrecursorTolr = mzUtils::ppmDist(targetMs1PrecMz, ms1PrecMz) <= params->ms3AnalysisMs1PrecursorPpmTolr;
                bool isMatchingMs2PrecursorTolr = mzUtils::ppmDist(targetMs2PrecMz, ms2PrecMz) <= params->ms3PrecursorPpmTolr;

                if (isMatchingMs1PrecursorTolr && isMatchingMs2PrecursorTolr) {

                    Fragment t;
                    t.precursorMz = ms2PrecMz;
                    t.mzs = it->second;
                    t.intensity_array = ms3Compound->ms3_fragment_intensity[it->first];
                    t.fragment_labels = ms3Compound->ms3_fragment_labels[it->first];

                    float maxDeltaMz = (params->ms3PpmTolr * static_cast<float>(t.precursorMz))/ 1000000;
                    vector<int> ranks = Fragment::findFragPairsGreedyMz(&t, get<2>(data)->consensus, maxDeltaMz);

                    bool isHasMatch = false;
                    for (unsigned long i = 0; i < ranks.size(); i++) {

                        int y = ranks[i];

                        if (y != -1) {
                            numMs3Matches++;
                            isHasMatch = true;
                            if (debug) {
                                matchInfoDebugString = matchInfoDebugString + "\t"
                                        + t.fragment_labels[i] + " " + to_string(t.mzs[i])
                                        + " <==> "
                                        + to_string(get<2>(data)->consensus->mzs[y]) + " (intensity=" + to_string(get<2>(data)->consensus->intensity_array[y]) + ")\n";
                            }
                        }
                    }

                    if (isHasMatch) {
                        //                        precMz,           <consensus Fragment*, ranks>
                        matchData.insert(make_pair(it->first, make_pair(get<2>(data), ranks)));
                    }
                } else if (targetMs1PrecMz > ms1MaxMz) {
                    //all subsequent consensus ms3 scans will not be associated with proper ms1 target m/z
                    break;
                }

            } // END ms3Compound m/z map

        } //END for (auto it = ms3Compound->ms3_fragment_mzs.begin();

        //Issue 240: ms1 precursor
        float observedMs1Intensity = 0.0f;

        if (params->ms1IsFindPrecursorIon && ms1Fragment && ms1Fragment->consensus) {
            double precMz = ms3Compound->baseCompound->precursorMz;

            double minMz = precMz - precMz*params->ms1PpmTolr/1e6;
            double maxMz = precMz + precMz*params->ms1PpmTolr/1e6;

            auto lb = lower_bound(ms1Fragment->consensus->mzs.begin(), ms1Fragment->consensus->mzs.end(), minMz);

            auto pos = lb - ms1Fragment->consensus->mzs.begin();

            for (unsigned int i = pos; i < ms1Fragment->consensus->mzs.size(); i++) {
                if (ms1Fragment->consensus->mzs[i] <= maxMz) {
                    if (ms1Fragment->consensus->intensity_array[i] > observedMs1Intensity) {
                        observedMs1Intensity = ms1Fragment->consensus->intensity_array[i];
                    }
                } else {
                    break;
                }
            }
        }

        bool isPassesMs1PrecursorRequirements = !params->ms1IsFindPrecursorIon || (observedMs1Intensity > 0.0f && observedMs1Intensity >= params->ms1MinIntensity);

        if (numMs3Matches >= params->ms3MinNumMatches && isPassesMs1PrecursorRequirements) {

            Ms3SingleSampleMatch *ms3SingleSampleMatch = new Ms3SingleSampleMatch;
            ms3SingleSampleMatch->ms3Compound = ms3Compound;
            ms3SingleSampleMatch->sample = sample;
            ms3SingleSampleMatch->numMs3Matches = numMs3Matches;
            ms3SingleSampleMatch->matchData = matchData;
            ms3SingleSampleMatch->observedMs1Intensity = observedMs1Intensity;

            output.push_back(ms3SingleSampleMatch);
            if (debug) cout << ms3Compound->baseCompound->name << " " << ms3Compound->baseCompound->adductString << ": " << numMs3Matches << " matches; observedMs1Intensity=" << ms3SingleSampleMatch->observedMs1Intensity << endl;
            if (debug) cout << matchInfoDebugString;
        }

        //Issue 244
        if (debug) cout << "Finished comparing compound #" << (compoundCounter+1)
                        << " (" << ms3Compound->baseCompound->name << " " << ms3Compound->baseCompound->adductString
                        << "). number of Ms3SingleSampleMatch matches so far: " << output.size() << endl;

        compoundCounter++;

    } // END for (auto ms3Compound : ms3Compounds)

    if (debug){
        long numOutputRows = 0;
        for (auto match : output) {
          numOutputRows += match->numMs3Matches;
        }
        cout << "Identified " << output.size() << " Ms3SingleSampleMatches, with a total of " << numOutputRows << " fragments matched." << endl;
    }

    return output;
}

map<int, DirectInfusionAnnotation*> DirectInfusionProcessor::processSingleSample(mzSample* sample,
                                                                              shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
                                                                              shared_ptr<DirectInfusionSearchParameters> params,
                                                                              bool debug) {

    MassCalculator massCalc;
    map<int, DirectInfusionAnnotation*> annotations = {};

    if (debug) cerr << "Started DirectInfusionProcessor::processSingleSample()" << endl;

    //Organize all scans by common precursor m/z

    map<int, vector<Scan*>> ms2ScansByBlockNumber = {};
    vector<Scan*> validMs1Scans;

    for (Scan* scan : sample->scans){

        if (scan->mslevel == 2 &&
                (params->scanFilterMs2MinRt <= -1.0f || scan->rt >= params->scanFilterMs2MinRt) &&
                (params->scanFilterMs2MaxRt <= -1.0f || scan->rt <= params->scanFilterMs2MaxRt)){

            int mapKey = static_cast<int>(round(scan->precursorMz+0.001f)); //round to nearest int

            if (ms2ScansByBlockNumber.find(mapKey) == ms2ScansByBlockNumber.end()) {
                ms2ScansByBlockNumber.insert(make_pair(mapKey, vector<Scan*>()));
            }

            ms2ScansByBlockNumber[mapKey].push_back(scan);
        }
        if (scan->mslevel == 1 &&
                scan->filterString.find(params->ms1ScanFilter) != string::npos &&
                (params->scanFilterMs1MinRt <= -1.0f || scan->rt >= params->scanFilterMs1MinRt) &&
                (params->scanFilterMs1MaxRt <= -1.0f || scan->rt <= params->scanFilterMs1MaxRt)) {
            validMs1Scans.push_back(scan);
        }
    }
    if (debug) cerr << "Performing search over map keys..." << endl;

    //For MS1 quant
    if (debug) cerr << "Computing consensus MS1 scan..." << endl;

    Fragment *ms1Fragment = nullptr;
    for (auto & scan: validMs1Scans) {
        if (!ms1Fragment) {
            ms1Fragment = new Fragment(scan,
                                       params->scanFilterMinFracIntensity,
                                       params->scanFilterMinSNRatio,
                                       params->scanFilterMaxNumberOfFragments,
                                       params->scanFilterBaseLinePercentile,
                                       params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                       params->scanFilterPrecursorPurityPpm,
                                       params->scanFilterMinIntensity);
        } else {
            Fragment *ms1Brother = new Fragment(scan,
                                             params->scanFilterMinFracIntensity,
                                             params->scanFilterMinSNRatio,
                                             params->scanFilterMaxNumberOfFragments,
                                             params->scanFilterBaseLinePercentile,
                                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                             params->scanFilterPrecursorPurityPpm,
                                             params->scanFilterMinIntensity);

            ms1Fragment->addFragment(ms1Brother);
        }
    }

    if (ms1Fragment){
        ms1Fragment->buildConsensus(params->consensusMs1PpmTolr,
                                    params->consensusIntensityAgglomerationType,
                                    params->consensusIsIntensityAvgByObserved,
                                    params->consensusIsNormalizeTo10K,
                                    params->consensusMinNumMs1Scans,
                                    params->consensusMinFractionMs1Scans
                                    );

        ms1Fragment->consensus->sortByMz();
    }

    if (debug) cerr << "Finished computing consensus MS1 scan." << endl;

    for (auto mapKey : directInfusionSearchSet->mapKeys){

        DirectInfusionAnnotation* directInfusionAnnotation = processBlock(mapKey,
                                                                          directInfusionSearchSet->mzRangesByMapKey[mapKey],
                                                                          sample,
                                                                          ms2ScansByBlockNumber[mapKey],
                                                                          ms1Fragment,
                                                                          directInfusionSearchSet->compoundsByMapKey[mapKey],
                                                                          params,
                                                                          debug);

        if (directInfusionAnnotation) annotations.insert(make_pair(mapKey, directInfusionAnnotation));

    }

    //Issue 232: prevent memory leak
    if (ms1Fragment) delete ms1Fragment;
    ms1Fragment = nullptr;

    return annotations;

}

DirectInfusionAnnotation* DirectInfusionProcessor::processBlock(int blockNum,
                                       const pair<float, float>& mzRange,
                                       mzSample* sample,
                                       const vector<Scan*>& ms2Scans,
                                       const Fragment *ms1Fragment,
                                       const vector<pair<Compound*, Adduct*>> library,
                                       const shared_ptr<DirectInfusionSearchParameters> params,
                                       const bool debug){

    //need MS2 scans and compounds to identify matches
    if (ms2Scans.empty()) return nullptr;
    if (library.empty()) return nullptr;

    //build search spectrum
    Fragment *f = nullptr;
    Scan* representativeScan = nullptr;
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
            representativeScan = scan;
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

    vector<shared_ptr<DirectInfusionMatchData>> libraryMatches;

    //Compare to library
    for (auto libraryEntry : library){

        unique_ptr<DirectInfusionMatchAssessment> matchAssessment = assessMatch(f, ms1Fragment, libraryEntry, params, debug);
        FragmentationMatchScore s = matchAssessment->fragmentationMatchScore;
        float fragmentMaxObservedIntensity = matchAssessment->fragmentMaxObservedIntensity;
        float observedMs1Intensity = matchAssessment->observedMs1Intensity;

        if (s.numMatches >= params->ms2MinNumMatches &&
                s.numDiagnosticMatches >= params->ms2MinNumDiagnosticMatches &&
                params->isDiagnosticFragmentMapAgreement(matchAssessment->diagnosticFragmentMatchMap)) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());

            directInfusionMatchData->compound = libraryEntry.first;
            directInfusionMatchData->adduct = libraryEntry.second;
            directInfusionMatchData->fragmentationMatchScore = s;
            directInfusionMatchData->fragmentMaxObservedIntensity = fragmentMaxObservedIntensity;
            directInfusionMatchData->observedMs1Intensity = observedMs1Intensity;

            libraryMatches.push_back(directInfusionMatchData);
        }
    }

    //agglomerate (if necessary), and return if valid matches exist.
    if (!libraryMatches.empty()){

        //output
        DirectInfusionAnnotation *directInfusionAnnotation = new DirectInfusionAnnotation();
        directInfusionAnnotation->precMzMin = mzRange.first;
        directInfusionAnnotation->precMzMax = mzRange.second;
        directInfusionAnnotation->sample = sample;
        directInfusionAnnotation->scan = representativeScan;
        directInfusionAnnotation->fragmentationPattern = f;

        //determine fragment match maps, and mutate compounds, if needed.
        //includes agglomerating compounds into SummarizedCompounds.
        unique_ptr<DirectInfusionMatchInformation> matchInfo = DirectInfusionProcessor::getMatchInformation(
                    libraryMatches,
                    f->consensus,
                    params,
                    debug);

        //TODO: relocate match info
        vector<shared_ptr<DirectInfusionMatchData>> processedMatchData(matchInfo->matchDataToFrags.size());

        unsigned int i = 0;
        for (auto it = matchInfo->matchDataToFrags.begin(); it != matchInfo->matchDataToFrags.end(); ++it){
            processedMatchData[i] = it->first;
            i++;
        }

        directInfusionAnnotation->compounds = processedMatchData;

        return directInfusionAnnotation;
    }

    return nullptr;
}

unique_ptr<DirectInfusionMatchAssessment> DirectInfusionProcessor::assessMatch(const Fragment *f, //ms2 fragment
                                                             const Fragment *ms1Fragment,
                                                             const pair<Compound*, Adduct*>& libraryEntry,
                                                             const shared_ptr<DirectInfusionSearchParameters> params,
                                                             const bool debug){
    //Initialize output
    unique_ptr<DirectInfusionMatchAssessment> directInfusionMatchAssessment = unique_ptr<DirectInfusionMatchAssessment>(new DirectInfusionMatchAssessment());

    Compound* compound = libraryEntry.first;
    Adduct *adduct = libraryEntry.second;

    //=============================================== //
    //START COMPARE MS1
    //=============================================== //

    float observedMs1Intensity = 0.0f;

    if (ms1Fragment && ms1Fragment->consensus) {
        double precMz = compound->precursorMz;
        if (!params->ms1IsRequireAdductPrecursorMatch) {

            //Compute this way instead of using compound->precursorMz to allow for possibility of matching compound to unexpected adduct
            MassCalculator massCalc;
            float compoundMz = adduct->computeAdductMass(massCalc.computeNeutralMass(compound->getFormula()));
            precMz = adduct->computeAdductMass(compoundMz);

        }

        double minMz = precMz - precMz*params->ms1PpmTolr/1e6;
        double maxMz = precMz + precMz*params->ms1PpmTolr/1e6;

        auto lb = lower_bound(ms1Fragment->consensus->mzs.begin(), ms1Fragment->consensus->mzs.end(), minMz);

        auto pos = lb - ms1Fragment->consensus->mzs.begin();

        for (unsigned int i = pos; i < ms1Fragment->consensus->mzs.size(); i++) {
            if (ms1Fragment->consensus->mzs[i] <= maxMz) {
                if (ms1Fragment->consensus->intensity_array[i] > observedMs1Intensity) {
                    observedMs1Intensity = ms1Fragment->consensus->intensity_array[i];
                }
            } else {
                break;
            }
        }
    }

    bool isPassesMs1PrecursorRequirements = !params->ms1IsFindPrecursorIon || (observedMs1Intensity > 0.0f && observedMs1Intensity >= params->ms1MinIntensity);

    if (!isPassesMs1PrecursorRequirements) return directInfusionMatchAssessment; // will return with no matching fragments, 0 for every score

    //=============================================== //
    //END COMPARE MS1
    //=============================================== //

    //=============================================== //
    //START COMPARE MS2
    //=============================================== //

    Fragment t;
    t.precursorMz = compound->precursorMz;
    t.mzs = compound->fragment_mzs;
    t.intensity_array = compound->fragment_intensity;
    t.fragment_labels = compound->fragment_labels;

    float maxDeltaMz = (params->ms2PpmTolr * static_cast<float>(t.precursorMz))/ 1000000;
    directInfusionMatchAssessment->fragmentationMatchScore.ranks = Fragment::findFragPairsGreedyMz(&t, f->consensus, maxDeltaMz);

    bool isHasLabels = compound->fragment_labels.size() == directInfusionMatchAssessment->fragmentationMatchScore.ranks.size();

    float fragmentMaxObservedIntensity = 0;

    map<string, int> diagnosticMatchesMap = {};
    for (auto it = params->ms2MinNumDiagnosticMatchesMap.begin(); it != params->ms2MinNumDiagnosticMatchesMap.end(); ++it){
        diagnosticMatchesMap.insert(make_pair(it->first, 0));
    }

    for (unsigned long i=0; i < directInfusionMatchAssessment->fragmentationMatchScore.ranks.size(); i++) {

        int y = directInfusionMatchAssessment->fragmentationMatchScore.ranks[i];

        if (y != -1 && f->consensus->intensity_array[y] >= params->ms2MinIntensity) {

            float fragmentObservedIntensity = f->consensus->intensity_array[y];

            if (fragmentObservedIntensity > fragmentMaxObservedIntensity) {
                fragmentMaxObservedIntensity = fragmentObservedIntensity;
            }

            directInfusionMatchAssessment->fragmentationMatchScore.numMatches++;

            if (!isHasLabels) continue;

            if (compound->fragment_labels[i].find("*") == 0) {
                directInfusionMatchAssessment->fragmentationMatchScore.numDiagnosticMatches++;
            }

            for (auto it = params->ms2MinNumDiagnosticMatchesMap.begin(); it != params->ms2MinNumDiagnosticMatchesMap.end(); ++it){
                string diagnosticFragLabel = it->first;
                if (compound->fragment_labels[i].find(diagnosticFragLabel) == 0) {
                    diagnosticMatchesMap[diagnosticFragLabel]++;
                }
            }

        }
    }

    directInfusionMatchAssessment->diagnosticFragmentMatchMap = diagnosticMatchesMap;
    directInfusionMatchAssessment->fragmentMaxObservedIntensity = fragmentMaxObservedIntensity;
    directInfusionMatchAssessment->observedMs1Intensity = observedMs1Intensity;

    //=============================================== //
    //END COMPARE MS2
    //=============================================== //

    return directInfusionMatchAssessment;
}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::getFragmentMatchMaps(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

       if (debug) cerr << "DirectInfusionProcessor::getFragmentMatchMaps()" << endl;

       unique_ptr<DirectInfusionMatchInformation> matchInfo = unique_ptr<DirectInfusionMatchInformation>(new DirectInfusionMatchInformation());

       for (auto directInfusionMatchData : allCandidates) {

           Compound *compound = directInfusionMatchData->compound;
           FragmentationMatchScore fragmentationMatchScore = directInfusionMatchData->fragmentationMatchScore;

           vector<int> compoundFrags(static_cast<unsigned int>(fragmentationMatchScore.numMatches));

           unsigned int matchCounter = 0;
           for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

               int observedIndex = fragmentationMatchScore.ranks[i];

               //Issue 209: peaks may be unmatched based on intensity as well as ranks[] position
               if (observedIndex == -1 || observedSpectrum->intensity_array[observedIndex] < params->ms2MinIntensity) continue;

               if (debug) cerr << "allCandidates[" << i << "]: " << compound->name << "|" << compound->adductString << " observedIndex=" << observedIndex << endl;

               int fragInt = mzToIntKey(compound->fragment_mzs[i]);

               compoundFrags[matchCounter] = fragInt;
               matchCounter++;

               pair<int, shared_ptr<DirectInfusionMatchData>> key = make_pair(fragInt, directInfusionMatchData);

               if (matchInfo->fragToMatchData.find(fragInt) == matchInfo->fragToMatchData.end()) {
                   matchInfo->fragToMatchData.insert(make_pair(fragInt, unordered_set<shared_ptr<DirectInfusionMatchData>>()));
               }

               matchInfo->fragToMatchData[fragInt].insert(directInfusionMatchData);

   //            if (debug) cerr << "allCandidates [end] i=" << i << ", ranks=" << fragmentationMatchScore.ranks.size() << endl;
           }

           //Issue 270: need these to be sorted
           sort(compoundFrags.begin(), compoundFrags.end());
           matchInfo->matchDataToFrags.insert(make_pair(directInfusionMatchData, compoundFrags));

       }

       //Issue 270: helpful to keep this here
       for (auto it = matchInfo->matchDataToFrags.begin(); it != matchInfo->matchDataToFrags.end(); ++it){
           vector<int> fragList = it->second;
           if (matchInfo->fragListToCompounds.find(fragList) == matchInfo->fragListToCompounds.end()) {
               matchInfo->fragListToCompounds.insert(make_pair(fragList, vector<shared_ptr<DirectInfusionMatchData>>()));
           }
           matchInfo->fragListToCompounds[fragList].push_back(it->first);
       }

       return matchInfo;

}

/**
 * @brief DirectInfusionProcessor::summarizeFragmentGroups
 * @param matchInfo
 * @param observedSpectrum
 * @param params
 * @param debug
 * @return
 *
 * A "fragment group" is the set of m/zs matched.
 *
 * The set of all fragments that match to exactly the same set of m/zs can be condensed into a summarized compound.
 *
 * Depending on the type of summarization, this may involve naming the fragment according to summed composition,
 * chain length summary, or in the most general, a collection of all compounds names with no adjustments.
 *
 */
unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::summarizeFragmentGroups(
        unique_ptr<DirectInfusionMatchInformation> matchInfo,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug) {

    if (debug) cout << "directInfusionProcessor::summarizeFragmentGroups()" << endl;

    map<vector<int>, vector<shared_ptr<DirectInfusionMatchData>>> summarizedFragListToCompounds{};
    map<shared_ptr<DirectInfusionMatchData>, vector<int>> summarizedMatchDataToFrags{};
    map<int, unordered_set<shared_ptr<DirectInfusionMatchData>>> summarizedFragToMatchData{};

    for (auto it = matchInfo->fragListToCompounds.begin(); it != matchInfo->fragListToCompounds.end(); ++it) {

        vector<shared_ptr<DirectInfusionMatchData>> compoundList = it->second;

        //Initialize output
        shared_ptr<DirectInfusionMatchData> summarizedMatchData;

        //Case: if only one match, do not try to do any summarization.
        if (compoundList.size() == 1) {
            summarizedMatchData = compoundList[0];
            summarizedMatchDataToFrags.insert(make_pair(summarizedMatchData, it->first));
            continue;
        }

        //Issue 267: consistent name order
        //Issue 272: avoid rounding errors by always computing average in the same order
        sort(compoundList.begin(), compoundList.end(), [](const shared_ptr<DirectInfusionMatchData>& lhs, const shared_ptr<DirectInfusionMatchData>& rhs){
           return lhs->compound->name < rhs->compound->name;
        });

        /*
         *  Determine how to summarize identical fragment matched compounds (acyl chain, composition, or general)
         */

        unordered_set<string> compositionLevel{};
        unordered_set<string> acylChainLevel{};

        bool isMissingCompositionLevel = false;
        bool isMissingAcylChainLevel = false;

        Adduct *adduct = nullptr;
        vector<Compound*> compounds;
        float observedMs1Intensity = 0.0f;

        for (auto matchData : compoundList) {

            compounds.push_back(matchData->compound);
            adduct = matchData->adduct;
            observedMs1Intensity += matchData->observedMs1Intensity;

            if (matchData->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()) != matchData->compound->metaDataMap.end()){
                compositionLevel.insert(matchData->compound->metaDataMap[LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()]);
            } else {
                isMissingCompositionLevel = true;
            }

            if (matchData->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()) != matchData->compound->metaDataMap.end()){
                acylChainLevel.insert(matchData->compound->metaDataMap[LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()]);
            } else {
                isMissingAcylChainLevel = true;
            }
        }
        observedMs1Intensity /= compoundList.size();

        bool isUseAcylChainOrCompositionSummarization = params->spectralCompositionAlgorithm == SpectralCompositionAlgorithm::AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION &&
                ((!isMissingAcylChainLevel && acylChainLevel.size() == 1) || (!isMissingCompositionLevel && compositionLevel.size() == 1));

        if (isUseAcylChainOrCompositionSummarization) {

            /*
             * try to summarize to acyl chain level, then composition level, then finally fall back to identical fragment matches,
             * if these summarizations are not possible.
             */

            SummarizedCompound *summarizedCompound = nullptr;

            if (!isMissingAcylChainLevel && acylChainLevel.size() == 1) {

                summarizedCompound = new SummarizedCompound(*acylChainLevel.begin(), compounds);

            } else if (!isMissingCompositionLevel && compositionLevel.size() == 1) {

                summarizedCompound = new SummarizedCompound(*compositionLevel.begin(), compounds);

            } else {
                cerr << "Problem in finding appropriate level of summarization for fragment matches!" << endl;
                abort();
            }

            //This is necessary, as the summarized form may actually map to multiple fragments (e.g., multiple TG(60:7)).
            string summarizedId("{" + summarizedCompound->name + summarizedCompound->adductString + "}={");
            for (unsigned int i = 0; i < compoundList.size(); i++) {
                if (i > 0) {
                    summarizedId += ";";
                }
                summarizedId += compoundList[i]->compound->name;
            }
            summarizedId += "}";

            summarizedCompound->adductString = compounds.at(0)->adductString;
            summarizedCompound->formula = compounds.at(0)->getFormula();
            summarizedCompound->precursorMz = compounds.at(0)->precursorMz;
            summarizedCompound->setExactMass(compounds.at(0)->getExactMass());
            summarizedCompound->charge = compounds.at(0)->charge;
            summarizedCompound->id = summarizedId;

            summarizedCompound->computeSummarizedData();

            summarizedMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());
            summarizedMatchData->compound = summarizedCompound;
            summarizedMatchData->adduct = adduct;
            summarizedMatchData->fragmentationMatchScore = summarizedCompound->scoreCompoundHit(observedSpectrum, params->ms2PpmTolr, false);

            summarizedMatchData->observedMs1Intensity = observedMs1Intensity;

        } else { //fall back to general summarization based on identical fragments

            if (debug) cout << "Creating general summarization based on identical fragments." << endl;

            string summarizedName("{");

            string adductString("");
            string formulaString("");
            double precursorMz = 0;
            float exactMass = 0;
            int charge = 0;

            map<string, int> adductStringMap{};
            map<string, int> formulaStringMap{};

            vector<Compound*> compoundPtrs(compoundList.size());

            for (unsigned int i = 0; i < compoundList.size(); i++) {
                shared_ptr<DirectInfusionMatchData> compound = compoundList[i];

                compoundPtrs[i] = compound->compound;

                if (i > 0) {
                    summarizedName += ";";
                }

                summarizedName = summarizedName + compound->compound->name + "|" + compound->adduct->name;

                precursorMz += compound->compound->precursorMz;
                exactMass += compound->compound->getExactMass();
                charge += compound->compound->charge;

                if (adductStringMap.find(compound->adduct->name) == adductStringMap.end()){
                    adductStringMap.insert(make_pair(compound->adduct->name, 0));
                }
                adductStringMap[compound->adduct->name]++;

                if (formulaStringMap.find(compound->compound->formula) == formulaStringMap.end()) {
                    formulaStringMap.insert(make_pair(compound->compound->formula, 0));
                }
                formulaStringMap[compound->compound->formula]++;
            }

            precursorMz /= compoundList.size();
            exactMass /= compoundList.size();
            charge = static_cast<int>(charge/compoundList.size());

            vector<pair<string, int>> adductNameCounts{};
            for (auto it = adductStringMap.begin(); it != adductStringMap.end(); ++it){
                adductNameCounts.push_back(make_pair(it->first, it->second));
            }

            vector<pair<string, int>> formulaNameCounts{};
            for (auto it = formulaStringMap.begin(); it != formulaStringMap.end(); ++it){
                formulaNameCounts.push_back(make_pair(it->first, it->second));
            }

            sort(adductNameCounts.begin(), adductNameCounts.end(), [](const pair<string, int>& lhs, const pair<string, int>& rhs){
                if (lhs.second == rhs.second) {
                    return lhs.first.compare(rhs.first);
                } else {
                    return lhs.second - rhs.second;
                }
            });

            sort(formulaNameCounts.begin(), formulaNameCounts.end(), [](const pair<string, int>& lhs, const pair<string, int>& rhs){
                if (lhs.second == rhs.second) {
                    return lhs.first.compare(rhs.first);
                } else {
                    return lhs.second - rhs.second;
                }
            });

            adductString = adductNameCounts[0].first;
            formulaString = formulaNameCounts[0].first;

            Adduct *adduct = nullptr;
            for (auto match : compoundList){
                if (match->adduct->name == adductString) {
                    adduct = match->adduct;
                    break;
                }
            }
            summarizedName += "}";

            SummarizedCompound *summarizedCompound = new SummarizedCompound(summarizedName, compoundPtrs);

            summarizedCompound->adductString = adductString;
            summarizedCompound->formula = formulaString;
            summarizedCompound->precursorMz = precursorMz;
            summarizedCompound->setExactMass(exactMass);
            summarizedCompound->charge = charge;
            summarizedCompound->id = summarizedName + adductString;

            summarizedCompound->computeSummarizedData();

            summarizedMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());
            summarizedMatchData->compound = summarizedCompound;
            summarizedMatchData->adduct = adduct;
            summarizedMatchData->fragmentationMatchScore = summarizedCompound->scoreCompoundHit(observedSpectrum, params->ms2PpmTolr, false);

            summarizedMatchData->observedMs1Intensity = observedMs1Intensity;
        }

        summarizedMatchDataToFrags.insert(make_pair(summarizedMatchData, it->first));
    }

    for (auto it = summarizedMatchDataToFrags.begin(); it != summarizedMatchDataToFrags.end(); ++it) {

        shared_ptr<DirectInfusionMatchData> matchData = it->first;
        vector<int> fragList = it->second;

        for (auto frag : fragList) {
            if (summarizedFragToMatchData.find(frag) == summarizedFragToMatchData.end()) {
                summarizedFragToMatchData.insert(make_pair(frag, unordered_set<shared_ptr<DirectInfusionMatchData>>()));
            }
            summarizedFragToMatchData[frag].insert(it->first);
        }

        //after summarization, every fragList corresponds to only one compound
        vector<shared_ptr<DirectInfusionMatchData>> matchDataVector = vector<shared_ptr<DirectInfusionMatchData>>(1);
        matchDataVector[0] = matchData;
        summarizedFragListToCompounds.insert(make_pair(fragList, matchDataVector));
    }

    matchInfo->fragListToCompounds = summarizedFragListToCompounds;
    matchInfo->fragToMatchData = summarizedFragToMatchData;
    matchInfo->matchDataToFrags = summarizedMatchDataToFrags;

    return matchInfo;

}

//Issue 270
unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::reduceBySimpleParsimony(
        unique_ptr<DirectInfusionMatchInformation> matchInfo,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    vector<vector<int>> fragmentGroups(matchInfo->fragListToCompounds.size());
    unsigned int counter = 0;
    for (auto it = matchInfo->fragListToCompounds.begin(); it != matchInfo->fragListToCompounds.end(); ++it) {
        fragmentGroups[counter] = it->first;
        counter++;
    }

    map<vector<int>, vector<shared_ptr<DirectInfusionMatchData>>> reducedFragListToCompounds{};

    map<shared_ptr<DirectInfusionMatchData>, vector<int>> reducedMatchDataToFrags{};
    map<int, unordered_set<shared_ptr<DirectInfusionMatchData>>> reducedFragToMatchData{};

    vector<vector<int>> fragmentGroupsReducedyByParsimony = mzUtils::simpleParsimonyReducer(fragmentGroups);

    for (vector<int> group : fragmentGroupsReducedyByParsimony) {

        vector<shared_ptr<DirectInfusionMatchData>> compounds = matchInfo->fragListToCompounds[group];

        reducedFragListToCompounds.insert(make_pair(group, compounds));

        for (shared_ptr<DirectInfusionMatchData> compound : compounds) {

            reducedMatchDataToFrags.insert(make_pair(compound, group));

            for (int frag : group) {

                if (reducedFragToMatchData.find(frag) == reducedFragToMatchData.end()) {
                    reducedFragToMatchData.insert(make_pair(frag, unordered_set<shared_ptr<DirectInfusionMatchData>>{}));
                }

                reducedFragToMatchData[frag].insert(compound);

            }
        }

    }

    matchInfo->fragListToCompounds = reducedFragListToCompounds;
    matchInfo->matchDataToFrags = reducedMatchDataToFrags;
    matchInfo->fragToMatchData = reducedFragToMatchData;

    return matchInfo;
}


unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::reduceByUniqueMatches(
            unique_ptr<DirectInfusionMatchInformation> matchInfo,
            Fragment* observedSpectrum,
            shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    map<vector<int>, vector<shared_ptr<DirectInfusionMatchData>>> reducedFragListToCompounds{};

    map<shared_ptr<DirectInfusionMatchData>, vector<int>> reducedMatchDataToFrags{};
    map<int, unordered_set<shared_ptr<DirectInfusionMatchData>>> reducedFragToMatchData = {};


    for (auto it = matchInfo->matchDataToFrags.begin(); it != matchInfo->matchDataToFrags.end(); ++it){

        shared_ptr<DirectInfusionMatchData> matchData = it->first;
        vector<int> matchDataFrags = it->second;

        if (matchData->numUniqueFragments >= params->ms2MinNumUniqueMatches) {

            reducedMatchDataToFrags.insert(make_pair(matchData, matchDataFrags));

            if (reducedFragListToCompounds.find(matchDataFrags) == reducedFragListToCompounds.end()) {
                reducedFragListToCompounds.insert(make_pair(matchDataFrags, vector<shared_ptr<DirectInfusionMatchData>>()));
            }
            reducedFragListToCompounds[matchDataFrags].push_back(matchData);

            for (auto frag : matchDataFrags) {
                if (reducedFragToMatchData.find(frag) == reducedFragToMatchData.end()) {
                    reducedFragToMatchData.insert(make_pair(frag, unordered_set<shared_ptr<DirectInfusionMatchData>>()));
                }
                reducedFragToMatchData[frag].insert(matchData);
            }
        }
    }

    matchInfo->fragListToCompounds = reducedFragListToCompounds;
    matchInfo->matchDataToFrags = reducedMatchDataToFrags;
    matchInfo->fragToMatchData = reducedFragToMatchData;

    return matchInfo;
}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::getMatchInformation(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    if (debug) cerr << "DirectInfusionProcessor::getMatchInformation()" << endl;

    unique_ptr<DirectInfusionMatchInformation> matchInfo = getFragmentMatchMaps(allCandidates, observedSpectrum, params, debug);

    if (params->isReduceBySimpleParsimony) {
        matchInfo = reduceBySimpleParsimony(move(matchInfo), params, debug);
    }

    //Issue 273: all summarization approaches based on identical fragment groups
    if (params->spectralCompositionAlgorithm != SpectralCompositionAlgorithm::ALL_CANDIDATES){
        matchInfo = summarizeFragmentGroups(move(matchInfo), observedSpectrum, params, debug);
    }

    addBlockSpecificMatchInfo(matchInfo.get(), observedSpectrum, params, debug);

    if (params->ms2MinNumUniqueMatches > 0) {
        matchInfo = reduceByUniqueMatches(move(matchInfo), observedSpectrum, params, debug);
    }

    if (debug) {
        cerr << "Fragments --> Compounds: (" << matchInfo->matchDataToFrags.size() << " passing compounds)" << endl;

        for (auto iterator = matchInfo->fragToMatchData.begin(); iterator != matchInfo->fragToMatchData.end(); ++iterator) {
            int frag = iterator->first;
            unordered_set<shared_ptr<DirectInfusionMatchData>> compounds = iterator->second;
            cerr<< "frag= " << intKeyToMz(frag) << " m/z : ";
            for (auto matchData : compounds) {
                cerr << matchData->compound->name << "|" << matchData->compound->adductString << " ";
            }
            cerr << endl;
        }

        cerr << "Compounds --> Fragments: (" << matchInfo->fragToMatchData.size() << " matched fragments)" << endl;

        for (auto iterator = matchInfo->matchDataToFrags.begin(); iterator != matchInfo->matchDataToFrags.end(); ++iterator) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = iterator->first;
            vector<int> frags = iterator->second;

            cerr << "Compound= " << directInfusionMatchData->compound->name << "|" << directInfusionMatchData->compound->adductString << ": ";
            for (auto frag : frags){
                cerr << intKeyToMz(frag) << " ";
            }
            cerr << endl;
        }
    }

    return matchInfo;
}


void DirectInfusionProcessor::addBlockSpecificMatchInfo(
        DirectInfusionMatchInformation *matchInfo,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    for (auto it = matchInfo->matchDataToFrags.begin(); it != matchInfo->matchDataToFrags.end(); ++it){
        shared_ptr<DirectInfusionMatchData> compound = it->first;
        compound->isFragmentUnique = vector<bool>(compound->compound->fragment_mzs.size(), false);
    }

    for (auto it = matchInfo->fragToMatchData.begin(); it != matchInfo->fragToMatchData.end(); ++it){

        int fragMzKey = it->first;
        unordered_set<shared_ptr<DirectInfusionMatchData>> compounds = it->second;

        if (compounds.size() == 1) { // unique fragment

            shared_ptr<DirectInfusionMatchData> matchData = *compounds.begin();
            double fragMzVal = intKeyToMz(fragMzKey) - 0.00001; // possible rounding error in int -> float conversion

            auto lb_it = lower_bound(matchData->compound->fragment_mzs.begin(), matchData->compound->fragment_mzs.end(), fragMzVal);
            long lb = lb_it - matchData->compound->fragment_mzs.begin();

            matchData->isFragmentUnique[lb] = true;
            matchData->numUniqueFragments++;
        }

    }
}

void DirectInfusionGroupAnnotation::clean() {
    for (map<mzSample*, DirectInfusionAnnotation*>::iterator it = annotationBySample.begin(); it != annotationBySample.end(); ++it) {

        if (it->second->fragmentationPattern) delete(it->second->fragmentationPattern);
        for (auto matchData : it->second->compounds) {
             //SummarizedCompounds are created transiently by directinfusionprocessor, Compounds are retrieved from DB.compounds
//            if (SummarizedCompound* sc = dynamic_cast<SummarizedCompound*>(matchData->compound)){
//                delete(sc);
//            }
        }
        if (it->second) delete(it->second);
    }
    annotationBySample.clear();
}

DirectInfusionGroupAnnotation* DirectInfusionGroupAnnotation::createByAverageProportions(vector<DirectInfusionAnnotation*> crossSampleAnnotations, shared_ptr<DirectInfusionSearchParameters> params, bool debug) {

    DirectInfusionGroupAnnotation *directInfusionGroupAnnotation = new DirectInfusionGroupAnnotation();

    directInfusionGroupAnnotation->precMzMin = crossSampleAnnotations.at(0)->precMzMin;
    directInfusionGroupAnnotation->precMzMax = crossSampleAnnotations.at(0)->precMzMax;

    if (debug) {
        cerr << "=========================================" << endl;
        cerr << "Merging peak groups in precMzRange = [" << directInfusionGroupAnnotation->precMzMin << " - " <<  directInfusionGroupAnnotation->precMzMax << "]" << endl;
    }

    Fragment *f = nullptr;

//    map<shared_ptr<DirectInfusionMatchData>, double, DirectInfusionMatchDataCompare> proportionSums = {};
//    map<shared_ptr<DirectInfusionMatchData>, FragmentationMatchScore, DirectInfusionMatchDataCompare> bestFragMatch = {};

    map<shared_ptr<DirectInfusionMatchData>, double, DirectInfusionMatchDataCompareByNames> proportionSums = {};
    map<shared_ptr<DirectInfusionMatchData>, FragmentationMatchScore, DirectInfusionMatchDataCompareByNames> bestFragMatch = {};

    unsigned int compoundInSampleMatchCounter = 0;

    /**
     * If a sample contains no compounds, it should be excluded from calculation for cross-sample adjusted proportions.
     *
     * Individual proportions of compounds within a single sample all sum to 1, except when no compounds are found in the sample,
     * in which case the individual proportions all sum to 0 (and should be excluded from re-calculating cross-sample contributions)
     */
    unsigned int numContributingSamples = 0;

    for (auto directInfusionAnnotation : crossSampleAnnotations){
        directInfusionGroupAnnotation->annotationBySample.insert(
                    make_pair(directInfusionAnnotation->sample,
                              directInfusionAnnotation)
                    );

        //Issue 218
        if (!f){
            f = new Fragment(directInfusionAnnotation->scan,
                             params->scanFilterMinFracIntensity,
                             params->scanFilterMinSNRatio,
                             params->scanFilterMaxNumberOfFragments,
                             params->scanFilterBaseLinePercentile,
                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                             params->scanFilterPrecursorPurityPpm,
                             params->scanFilterMinIntensity);
        } else {
            Fragment *brother = new Fragment(directInfusionAnnotation->scan,
                                             params->scanFilterMinFracIntensity,
                                             params->scanFilterMinSNRatio,
                                             params->scanFilterMaxNumberOfFragments,
                                             params->scanFilterBaseLinePercentile,
                                             params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                             params->scanFilterPrecursorPurityPpm,
                                             params->scanFilterMinIntensity);

            f->addFragment(brother);
        }

        //Issue 218
        if (directInfusionAnnotation->fragmentationPattern) {
            for (auto fragment : directInfusionAnnotation->fragmentationPattern->brothers) {
                if (fragment) {
                    Fragment *brother = new Fragment(fragment);
                    f->addFragment(brother);
                }
            }
        }

        if (debug) {
            cerr << "sample=" << directInfusionAnnotation->sample->sampleName
                 << ": " << directInfusionAnnotation->compounds.size() << " compounds."
                 << endl;
        }

        if (directInfusionAnnotation->compounds.size() > 0) {
            numContributingSamples++;
        }

        for (auto matchData : directInfusionAnnotation->compounds){

            compoundInSampleMatchCounter++;

            double runningSum = matchData->proportion;

            if (proportionSums.find(matchData) != proportionSums.end()){
                runningSum += proportionSums.at(matchData);
            } else {
                proportionSums.insert(make_pair(matchData, 0.0));
            }

            proportionSums.at(matchData) = runningSum;

            if (debug) {
                cerr << "sample=" << directInfusionAnnotation->sample->sampleName
                     << "(" << matchData->compound->name
                     << ", " << matchData->adduct->name
                     << ", proportion=" << matchData->proportion
                     << "): runningSum=" << runningSum << endl;
            }

            FragmentationMatchScore bestMatch = matchData->fragmentationMatchScore;
            if (bestFragMatch.find(matchData) != bestFragMatch.end()) {
                FragmentationMatchScore previousBestMatch = bestFragMatch.at(matchData);

                //TODO: how to decide on best match?
                if (bestMatch.hypergeomScore >= previousBestMatch.hypergeomScore){
                    bestFragMatch.at(matchData) =  bestMatch;
                }

            } else {
                bestFragMatch.insert(make_pair(matchData, bestMatch));
            }

        }

    }

    if (debug) {
        cerr << "Identified " << proportionSums.size() << " unique and " << compoundInSampleMatchCounter << " total compound-adduct pairs in all samples." << endl;
    }

    //Issue 218
    f->buildConsensus(params->consensusPpmTolr,
                      params->consensusIntensityAgglomerationType,
                      params->consensusIsIntensityAvgByObserved,
                      params->consensusIsNormalizeTo10K,
                      params->consensusMinNumMs2Scans,
                      params->consensusMinFractionMs2Scans
                      );
    f->consensus->sortByMz();

    directInfusionGroupAnnotation->fragmentationPattern = f;

    directInfusionGroupAnnotation->compounds.resize(proportionSums.size());

    double numSamples = static_cast<double>(directInfusionGroupAnnotation->annotationBySample.size());

    unsigned int annotationMatchIndex = 0;

    for (auto matchDataPair : proportionSums) {
       shared_ptr<DirectInfusionMatchData> groupMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());

       shared_ptr<DirectInfusionMatchData> matchData = matchDataPair.first;

       groupMatchData->compound = matchData->compound;
       groupMatchData->adduct = matchData->adduct;
       groupMatchData->proportion = matchDataPair.second / numContributingSamples;
       groupMatchData->fragmentationMatchScore = bestFragMatch.at(matchData);
       groupMatchData->fragmentMaxObservedIntensity = matchData->fragmentMaxObservedIntensity;
       groupMatchData->observedMs1Intensity = matchData->observedMs1Intensity;

       directInfusionGroupAnnotation->compounds.at(annotationMatchIndex) = groupMatchData;

       if (debug) {
           cerr << "Compound: " << groupMatchData->compound->name
                << ", Adduct: "<< groupMatchData->adduct->name
                << ", numMatches: " << groupMatchData->fragmentationMatchScore.numMatches
                << ", Proportion: " << groupMatchData->proportion
                << ", Fragment Max Observed Intensity: " << groupMatchData->fragmentMaxObservedIntensity
                << endl;
       }

       annotationMatchIndex++;
    }

    if (debug) {
        cerr << "Determined cross-sample proportions for " << annotationMatchIndex << " compounds." << endl;
        cerr << "=========================================" << endl;
    }

    return directInfusionGroupAnnotation;
}

