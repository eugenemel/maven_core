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

            bool isWroteCompoundToMs2ScansRange = false;
            for (mzRangeIterator it = directInfusionSearchSet->mzRangesByMapKey.begin(); it != directInfusionSearchSet->mzRangesByMapKey.end(); ++it) {
                int mapKey = it->first;
                pair<float, float> mzRange = it->second;

                if (compoundMz > mzRange.first && compoundMz < mzRange.second) {

                    if (directInfusionSearchSet->compoundsByMapKey.find(mapKey) == directInfusionSearchSet->compoundsByMapKey.end()){
                        directInfusionSearchSet->compoundsByMapKey.insert(make_pair(mapKey, vector<pair<Compound*, Adduct*>>()));
                    }

                    isWroteCompoundToMs2ScansRange = true;
                    directInfusionSearchSet->compoundsByMapKey[mapKey].push_back(make_pair(compound, adduct));
                    break;
                }
            }

            if (!isWroteCompoundToMs2ScansRange) {
                if (directInfusionSearchSet->compoundsByMapKey.find(DirectInfusionSearchSet::getNoMs2ScansMapKey()) == directInfusionSearchSet->compoundsByMapKey.end()) {
                    directInfusionSearchSet->compoundsByMapKey.insert(make_pair(DirectInfusionSearchSet::getNoMs2ScansMapKey(), vector<pair<Compound*, Adduct*>>()));
                    directInfusionSearchSet->mapKeys.insert(DirectInfusionSearchSet::getNoMs2ScansMapKey());
                }
                directInfusionSearchSet->compoundsByMapKey[DirectInfusionSearchSet::getNoMs2ScansMapKey()].push_back(make_pair(compound, adduct));
            }
        }
    }

    if (directInfusionSearchSet->compoundsByMapKey.find(DirectInfusionSearchSet::getNoMs2ScansMapKey()) != directInfusionSearchSet->compoundsByMapKey.end()) {

        //Currently, support any possible precursor m/z for compounds missing MS2 scans,
        //however this could also be limited by the MS1 scans in the experiment
        directInfusionSearchSet->mzRangesByMapKey.insert(make_pair(DirectInfusionSearchSet::getNoMs2ScansMapKey(), make_pair(0.0f, FLT_MAX)));
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

    //restructure as map for more efficient search
    map<pair<double, double>, vector<Scan*>> ms3ScanGroupMap{};

    for (auto ms3ScanGroup : ms3ScanGroups) {

        double ms1PrecMz = 0;
        double ms2PrecMz = 0;

        vector<Scan*> scans(ms3ScanGroup.size());
        for (unsigned int i = 0; i < ms3ScanGroup.size(); i++) {

            if (i == 0) {
                ms1PrecMz = get<0>(ms3ScanGroup[i]);
                ms2PrecMz = get<1>(ms3ScanGroup[i]);
            }

            scans[i] = get<2>(ms3ScanGroup[i]);

        }

        ms3ScanGroupMap.insert(make_pair(make_pair(ms1PrecMz, ms2PrecMz), scans));

        if (debug) cout << "(" << ms1PrecMz << ", " << ms2PrecMz << "): " << scans.size() << " scans." << endl;

    }

    unsigned int compoundCounter = 0;

    for (auto ms3Compound : ms3Compounds) {

        map<int, vector<float>> scanIntensitiesByMs3Mz{};
        map<pair<int, int>, vector<float>> scanIntensitiesByMs1Ms2Ms3Mzs{};

        string matchInfoDebugString("");

        double ms1PrecMz = ms3Compound->baseCompound->precursorMz;

        double ms1MinMz = ms1PrecMz - ms1PrecMz * params->ms3AnalysisMs1PrecursorPpmTolr/1000000.0;
        double ms1MaxMz = ms1PrecMz + ms1PrecMz * params->ms3AnalysisMs1PrecursorPpmTolr/1000000.0;

        for (auto it = ms3Compound->ms3_fragment_mzs.begin(); it != ms3Compound->ms3_fragment_mzs.end(); ++it){

            int ms2MzKey = it->first;

            double ms2PrecMz = mzUtils::intKeyToMz(ms2MzKey);

            for (unsigned int i = 0; i < it->second.size(); i++) {

                float ms3_mz = it->second[i];

                double ms3_mz_min = ms3_mz - params->ms3MatchTolrInDa;
                double ms3_mz_max = ms3_mz + params->ms3MatchTolrInDa;

                long ms3MzKey = mzUtils::mzToIntKey(ms3_mz);

                for (auto it2 = ms3ScanGroupMap.begin(); it2 != ms3ScanGroupMap.end(); ++it2) {

                    double targetMs1PrecMz = it2->first.first;
                    double targetMs2PrecMz = it2->first.second;

                    bool isMatchingMs1PrecursorTolr = mzUtils::ppmDist(targetMs1PrecMz, ms1PrecMz) <= params->ms3AnalysisMs1PrecursorPpmTolr;
                    bool isMatchingMs2PrecursorTolr = mzUtils::ppmDist(targetMs2PrecMz, ms2PrecMz) <= params->ms3PrecursorPpmTolr;

                    if (isMatchingMs1PrecursorTolr && isMatchingMs2PrecursorTolr) {

                        vector<Scan*> scans = it2->second;

                        vector<float> ms3Intensities{};

                        for (auto scan : scans) {

                          auto lb_ms3 = lower_bound(scan->mz.begin(), scan->mz.end(), ms3_mz_min);

                          float ms3_intensity = 0.0f;
                          float deltaMz = 99999;
                          bool isFoundMatch = false;

                          for (unsigned int ms3_pos = lb_ms3 - scan->mz.begin(); ms3_pos < scan->mz.size(); ms3_pos++) {

                            if (scan->mz[ms3_pos] > ms3_mz_max) {
                              break;
                            }

                            if (params->ms3IntensityType == Ms3IntensityType::ALL_MATCHES) {
                              ms3_intensity += scan->intensity[ms3_pos];
                              isFoundMatch = true;
                            } else if (params->ms3IntensityType == Ms3IntensityType::MAX_INTENSITY) {
                              if (scan->intensity[ms3_pos] > ms3_intensity) {
                                ms3_intensity = scan->intensity[ms3_pos];
                                isFoundMatch = true;
                              }
                            } else if (params->ms3IntensityType == Ms3IntensityType::CLOSEST_MZ) {
                              if (abs(scan->mz[ms3_pos] - ms3_mz) < deltaMz) {
                                deltaMz = abs(scan->mz[ms3_pos] - ms3_mz);
                                ms3_intensity = scan->mz[ms3_pos];
                                isFoundMatch = true;
                              }
                            }

                          } // END scans->mz

                          if (ms3_intensity >= params->ms3MinIntensity && isFoundMatch) {
                              ms3Intensities.push_back(ms3_intensity);
                          }

                        } // END scans

                        float ms3IntensityFraction = static_cast<float>(ms3Intensities.size())/static_cast<float>(scans.size());

                        if (ms3IntensityFraction >= params->ms3MinFractionScans && static_cast<int>(ms3Intensities.size()) >= params->ms3MinNumScans) {
                            for (auto ms3_intensity : ms3Intensities) {
                                if (scanIntensitiesByMs3Mz.find(ms3MzKey) == scanIntensitiesByMs3Mz.end()) {
                                    scanIntensitiesByMs3Mz.insert(make_pair(ms3MzKey, vector<float>()));
                                }
                                scanIntensitiesByMs3Mz[ms3MzKey].push_back(ms3_intensity);

                                pair<int, int> mzKey(ms2MzKey, i);
                                if (scanIntensitiesByMs1Ms2Ms3Mzs.find(mzKey) == scanIntensitiesByMs1Ms2Ms3Mzs.end()) {
                                    scanIntensitiesByMs1Ms2Ms3Mzs.insert(make_pair(mzKey, vector<float>()));
                                }
                                scanIntensitiesByMs1Ms2Ms3Mzs[mzKey].push_back(ms3_intensity);
                            }
                        }

                    } //END matching precursor m/z

                } // END ms3ScanGroupMap

            } // END ms3Compound->ms3_fragment_mzs vector<float>

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

        if (scanIntensitiesByMs1Ms2Ms3Mzs.size() >= params->ms3MinNumMatches && isPassesMs1PrecursorRequirements) {

            map<pair<int, int>, float> intensityByMs1Ms2Ms3Mzs{};

            for (auto it = scanIntensitiesByMs1Ms2Ms3Mzs.begin(); it != scanIntensitiesByMs1Ms2Ms3Mzs.end(); ++it) {
                sort(it->second.begin(), it->second.end());

                float ms3MzIntensity = 0.0f;

                if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
                    ms3MzIntensity = accumulate(it->second.begin(), it->second.end(), 0.0f) / it->second.size();
                } else if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
                    ms3MzIntensity = median(it->second);
                }

                intensityByMs1Ms2Ms3Mzs.insert(make_pair(it->first, ms3MzIntensity));

            }

            // Issue 226
            map<int, float> sumMs3IntensityByMs2Mz{};
            map<int, unordered_set<int>> matchingCoordsByMs2{};

            for (auto it = intensityByMs1Ms2Ms3Mzs.begin(); it != intensityByMs1Ms2Ms3Mzs.end(); ++it) {
                int ms2MzKey = it->first.first;

                if (sumMs3IntensityByMs2Mz.find(ms2MzKey) == sumMs3IntensityByMs2Mz.end()) {
                    sumMs3IntensityByMs2Mz.insert(make_pair(ms2MzKey, 0.0f));
                }
                if (matchingCoordsByMs2.find(ms2MzKey) == matchingCoordsByMs2.end()){
                    matchingCoordsByMs2.insert(make_pair(ms2MzKey, unordered_set<int>()));
                }

                sumMs3IntensityByMs2Mz[ms2MzKey] += it->second; //agglomerated measured intensity value
                matchingCoordsByMs2[ms2MzKey].insert(it->first.second); //coord in reference ms2 spectrum
            }

            map<int, int> ms3MatchesByMs2Mz{};
            for (auto it = matchingCoordsByMs2.begin(); it != matchingCoordsByMs2.end(); ++it){
                ms3MatchesByMs2Mz.insert(make_pair(it->first, it->second.size()));
            }

            map<int, float> intensityByMs3Mz{};

            float sumMs3MzIntensity = 0.0f;

            for (auto it = scanIntensitiesByMs3Mz.begin(); it != scanIntensitiesByMs3Mz.end(); ++it) {
                sort(it->second.begin(), it->second.end());

                float ms3MzIntensity = 0.0f;

                if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
                    ms3MzIntensity = accumulate(it->second.begin(), it->second.end(), 0.0f) / it->second.size();
                } else if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
                    ms3MzIntensity = median(it->second);
                }

                intensityByMs3Mz.insert(make_pair(it->first, ms3MzIntensity));

                sumMs3MzIntensity += ms3MzIntensity;

            }

            //Issue 296: skip IDs with too few params matching
            if (intensityByMs3Mz.size() < params->ms3MinNumMs3MzMatches) continue;

            Ms3SingleSampleMatch *ms3SingleSampleMatch = new Ms3SingleSampleMatch;
            ms3SingleSampleMatch->ms3Compound = ms3Compound;
            ms3SingleSampleMatch->sample = sample;
            ms3SingleSampleMatch->numMs3Matches = intensityByMs1Ms2Ms3Mzs.size();
            ms3SingleSampleMatch->numMs3MzMatches = intensityByMs3Mz.size();
            ms3SingleSampleMatch->observedMs1Intensity = observedMs1Intensity;
            ms3SingleSampleMatch->scanIntensitiesByMs1Ms2Ms3Mzs = scanIntensitiesByMs1Ms2Ms3Mzs;
            ms3SingleSampleMatch->scanIntensitiesByMs3Mz = scanIntensitiesByMs3Mz;
            ms3SingleSampleMatch->intensityByMs1Ms2Ms3Mzs = intensityByMs1Ms2Ms3Mzs;
            ms3SingleSampleMatch->intensityByMs3Mz = intensityByMs3Mz;
            ms3SingleSampleMatch->sumMs3MzIntensity = sumMs3MzIntensity;
            ms3SingleSampleMatch->matchingCoordsByMs2 = matchingCoordsByMs2;
            ms3SingleSampleMatch->ms3MatchesByMs2Mz = ms3MatchesByMs2Mz;
            ms3SingleSampleMatch->sumMs3IntensityByMs2Mz = sumMs3IntensityByMs2Mz;

            output.push_back(ms3SingleSampleMatch);
            if (debug) cout << ms3Compound->baseCompound->name << " " << ms3Compound->baseCompound->adductString << ": " << intensityByMs1Ms2Ms3Mzs.size() << " matches; observedMs1Intensity=" << ms3SingleSampleMatch->observedMs1Intensity << endl;
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

    //Issue 313
    if (directInfusionSearchSet->mapKeys.find(DirectInfusionSearchSet::getNoMs2ScansMapKey()) != directInfusionSearchSet->mapKeys.end()) {
        ms2ScansByBlockNumber.insert(make_pair(DirectInfusionSearchSet::getNoMs2ScansMapKey(), vector<Scan*>()));
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
                                                                          validMs1Scans,
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
                                       const vector<Scan*>& ms1Scans,
                                       const vector<Scan*>& ms2Scans,
                                       const Fragment *ms1Fragment,
                                       const vector<pair<Compound*, Adduct*>> library,
                                       const shared_ptr<DirectInfusionSearchParameters> params,
                                       const bool debug){

    if (debug) cout << "DirectInfusionProcessor::processBlock(): block #" << blockNum << endl;

    //need valid compounds to identify matches
    if (library.empty()) return nullptr;

    //build search spectrum
    Fragment *ms2Fragment = nullptr;
    Scan* representativeScan = nullptr;
    for (auto& scan : ms2Scans) {
        if (!ms2Fragment){
            ms2Fragment = new Fragment(scan,
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

            ms2Fragment->addFragment(brother);
        }
    }

    if(ms2Fragment) {
        ms2Fragment->buildConsensus(params->consensusPpmTolr,
                          params->consensusIntensityAgglomerationType,
                          params->consensusIsIntensityAvgByObserved,
                          params->consensusIsNormalizeTo10K,
                          params->consensusMinNumMs2Scans,
                          params->consensusMinFractionMs2Scans
                          );

        ms2Fragment->consensus->sortByMz();
    }

    vector<shared_ptr<DirectInfusionMatchData>> libraryMatches;

    //Compare data to library
    for (auto libraryEntry : library){

        Compound *compound = libraryEntry.first;
        Adduct *adduct = libraryEntry.second;

        int minNumMatches = params->ms2MinNumMatches;
        int minNumDiagnosticMatches = params->ms2MinNumDiagnosticMatches;
        int minNumSn1Matches = params->ms2sn1MinNumMatches;
        int minNumSn2Matches = params->ms2sn2MinNumMatches;

        string lipidClass;

        //Issue 316: ensure that compounds have all metadata by this point
        if (debug) {
            cout << "Metadata for "<< compound->name << ": " << endl;
            for (auto it = compound->metaDataMap.begin(); it != compound->metaDataMap.end(); ++it){
                cout << it->first << ": " << it->second << endl;
            }
        }

        //Issue 316: check for lipid class specific, or lipid class and adduct specific search criteria.
        if (compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != compound->metaDataMap.end()) {
            lipidClass = compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];

            string adductName = adduct->name;

            pair<string, string> lipidClassAndAdductKey = make_pair(lipidClass, adductName);
            pair<string, string> lipidClassKey = make_pair(lipidClass, "*");

            if (params->ms2MinNumMatchesByLipidClassAndAdduct.find(lipidClassAndAdductKey) != params->ms2MinNumMatchesByLipidClassAndAdduct.end()) {
                minNumMatches = params->ms2MinNumMatchesByLipidClassAndAdduct[lipidClassAndAdductKey];
            } else if (params->ms2MinNumMatchesByLipidClassAndAdduct.find(lipidClassKey) != params->ms2MinNumMatchesByLipidClassAndAdduct.end()) {
                minNumMatches = params->ms2MinNumMatchesByLipidClassAndAdduct[lipidClassKey];
            }

            if (debug) cout << "ms2MinNumMatches for lipidClass=" << lipidClass << ", adduct=" << adductName << ": " << minNumMatches << endl;

            if (params->ms2MinNumDiagnosticMatchesByLipidClassAndAdduct.find(lipidClassAndAdductKey) != params->ms2MinNumDiagnosticMatchesByLipidClassAndAdduct.end()) {
                minNumDiagnosticMatches = params->ms2MinNumDiagnosticMatchesByLipidClassAndAdduct[lipidClassAndAdductKey];
            } else if (params->ms2MinNumDiagnosticMatchesByLipidClassAndAdduct.find(lipidClassKey) != params->ms2MinNumDiagnosticMatchesByLipidClassAndAdduct.end()) {
                minNumDiagnosticMatches = params->ms2MinNumDiagnosticMatchesByLipidClassAndAdduct[lipidClassKey];
            }

            if (debug) cout << "ms2MinNumDiagnosticMatches for lipidClass=" << lipidClass << ", adduct=" << adductName << ": " << minNumDiagnosticMatches << endl;

            //Issue 359
            if (params->ms2sn1MinNumMatchesByLipidClassAndAdduct.find(lipidClassAndAdductKey) != params->ms2sn1MinNumMatchesByLipidClassAndAdduct.end()) {
                minNumSn1Matches = params->ms2sn1MinNumMatchesByLipidClassAndAdduct[lipidClassAndAdductKey];
            } else if (params->ms2sn1MinNumMatchesByLipidClassAndAdduct.find(lipidClassKey) != params->ms2sn1MinNumMatchesByLipidClassAndAdduct.end()) {
                minNumSn1Matches = params->ms2sn1MinNumMatchesByLipidClassAndAdduct[lipidClassKey];
            }

            if (debug) cout << "ms2MinNumSn1Matches for lipidClass=" << lipidClass << ", adduct=" << adductName << ": " << minNumSn1Matches << endl;

            if (params->ms2sn2MinNumMatchesByLipidClassAndAdduct.find(lipidClassAndAdductKey) != params->ms2sn2MinNumMatchesByLipidClassAndAdduct.end()) {
                minNumSn2Matches = params->ms2sn2MinNumMatchesByLipidClassAndAdduct[lipidClassAndAdductKey];
            } else if (params->ms2sn2MinNumMatchesByLipidClassAndAdduct.find(lipidClassKey) != params->ms2sn2MinNumMatchesByLipidClassAndAdduct.end()) {
                minNumSn2Matches = params->ms2sn2MinNumMatchesByLipidClassAndAdduct[lipidClassKey];
            }

            if (debug) cout << "ms2MinNumSn2Matches for lipidClass=" << lipidClass << ", adduct=" << adductName << ": " << minNumSn2Matches << endl;
        }

        unique_ptr<DirectInfusionMatchAssessment> matchAssessment = assessMatch(ms1Scans, ms1Fragment, ms2Fragment, libraryEntry, params, debug);

        FragmentationMatchScore s = matchAssessment->fragmentationMatchScore;
        float fragmentMaxObservedIntensity = matchAssessment->fragmentMaxObservedIntensity;
        float observedMs1Intensity = matchAssessment->observedMs1Intensity;
        int ms1IntensityCoord = matchAssessment->ms1IntensityCoord;
        float observedMs1ScanIntensity = matchAssessment->observedMs1ScanIntensity;

        //individual compound matches
        if (!matchAssessment->isDisqualifyThisMatch &&
                s.numMatches >= minNumMatches &&
                s.numDiagnosticMatches >= minNumDiagnosticMatches &&
                s.numSn1Matches >= minNumSn1Matches &&
                s.numSn2Matches >= minNumSn2Matches &&
                params->isDiagnosticFragmentMapAgreement(matchAssessment->diagnosticFragmentMatchMap)) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());

            directInfusionMatchData->compound = compound;
            directInfusionMatchData->adduct = adduct;
            directInfusionMatchData->fragmentationMatchScore = s;
            directInfusionMatchData->fragmentMaxObservedIntensity = fragmentMaxObservedIntensity;
            directInfusionMatchData->observedMs1Intensity = observedMs1Intensity;
            directInfusionMatchData->ms1IntensityCoord = ms1IntensityCoord;
            directInfusionMatchData->observedMs1ScanIntensity = observedMs1ScanIntensity;

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
        directInfusionAnnotation->fragmentationPattern = ms2Fragment;

        //Issue 113: Support ms1-only searches
        Fragment *ms2ConsensusSpectrum = nullptr;
        if (ms2Fragment && ms2Fragment->consensus) ms2ConsensusSpectrum = ms2Fragment->consensus;

        //determine fragment match maps, apply filters, and agglomerate compounds (if needed)
        unique_ptr<DirectInfusionMatchInformation> matchInfo = DirectInfusionProcessor::getMatchInformation(
                    libraryMatches,
                    ms2ConsensusSpectrum,
                    params,
                    debug);

        //Issue 288
        if (!ms2Scans.empty()) {
//            matchInfo->computeMs1PartitionFractions(ms2Scans, ms2Fragment, params, debug);

            //Issue 388: alternative computation approach
            matchInfo->computeMs1PartitionFractions2(ms2Fragment, params, debug);
        }

        directInfusionAnnotation->compounds = matchInfo->getCompounds();
        directInfusionAnnotation->matchInformation = move(matchInfo);

        return directInfusionAnnotation;
    }

    return nullptr;
}

unique_ptr<DirectInfusionMatchAssessment> DirectInfusionProcessor::assessMatch(
        const vector<Scan*>& ms1Scans,
        const Fragment *ms1Fragment,
        const Fragment *ms2Fragment,
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
    int ms1IntensityCoord = -1;

    float precMz = compound->precursorMz;

    if (!params->ms1IsRequireAdductPrecursorMatch) {

        //Compute this way instead of using compound->precursorMz to allow for possibility of matching compound to unexpected adduct
        MassCalculator massCalc;
        float compoundMz = adduct->computeAdductMass(static_cast<float>(massCalc.computeNeutralMass(compound->getFormula())));
        precMz = adduct->computeAdductMass(compoundMz);
    }

    if (ms1Fragment && ms1Fragment->consensus) {

        float minMz = precMz - precMz*params->ms1PpmTolr/1e6f;
        float maxMz = precMz + precMz*params->ms1PpmTolr/1e6f;

        auto lb = lower_bound(ms1Fragment->consensus->mzs.begin(), ms1Fragment->consensus->mzs.end(), minMz);

        auto pos = lb - ms1Fragment->consensus->mzs.begin();

        for (unsigned int i = pos; i < ms1Fragment->consensus->mzs.size(); i++) {
            if (ms1Fragment->consensus->mzs[i] <= maxMz) {
                if (ms1Fragment->consensus->intensity_array[i] > observedMs1Intensity) {
                    observedMs1Intensity = ms1Fragment->consensus->intensity_array[i];
                    ms1IntensityCoord = static_cast<int>(i);
                }
            } else {
                break;
            }
        }

        //Issue 303
        if (observedMs1Intensity > 0.0f && params->ms1IsRequireMonoisotopic) {

            float minMMinusOneMz = minMz - static_cast<float>(DirectInfusionUtils::C_13_MASS);
            float maxMMinusOneMz = maxMz + static_cast<float>(DirectInfusionUtils::C_13_MASS);

            auto lb = lower_bound(ms1Fragment->consensus->mzs.begin(), ms1Fragment->consensus->mzs.end(), minMMinusOneMz);

            auto pos = lb - ms1Fragment->consensus->mzs.begin();

            float maxMinusOneIntensity = -1.0f;
            for (unsigned int i = pos; i < ms1Fragment->consensus->mzs.size(); i++) {
                if (ms1Fragment->consensus->mzs[i] <= maxMMinusOneMz) {
                    if (ms1Fragment->consensus->intensity_array[i] > maxMinusOneIntensity) {
                        maxMinusOneIntensity = ms1Fragment->consensus->intensity_array[i];
                    }
                } else {
                    break;
                }
            }

            float mMinusOneToCandidateFraction = maxMinusOneIntensity / observedMs1Intensity;
            if (mMinusOneToCandidateFraction > params->ms1MMinusOnePeakMaxIntensityFraction) {
                observedMs1Intensity = 0.0f; //disqualify observed ms1 intensity based on presence of [M-1] peak
            }
        }
    }

    float observedMs1ScanIntensity = DirectInfusionUtils::getObservedMs1ScanIntensity(ms1Scans, precMz, params, debug);

    //an intensity of 0 always indicates absence, irrespective of the params->ms1MinIntensity value.
    bool isMs1InConsensusSpectrum = (observedMs1Intensity >= params->ms1MinIntensity && observedMs1Intensity > 0.0f);

    bool isFoundMs1PrecursorIonInScans = observedMs1ScanIntensity > 0.0f;

    if (params->ms1IsFindPrecursorIon && !(isMs1InConsensusSpectrum || isFoundMs1PrecursorIonInScans)){
        directInfusionMatchAssessment->isDisqualifyThisMatch = true;
        return directInfusionMatchAssessment; // will return with no matching fragments, 0 for every score.
     }

    directInfusionMatchAssessment->observedMs1ScanIntensity = observedMs1ScanIntensity;
    directInfusionMatchAssessment->observedMs1Intensity = observedMs1Intensity;
    directInfusionMatchAssessment->ms1IntensityCoord = ms1IntensityCoord;

    if (debug) {
        cout << compound->id << ": observedMs1Intensity=" << observedMs1Intensity << ", ms1IntensityCoord=" << ms1IntensityCoord << endl;
    }

    //=============================================== //
    //END COMPARE MS1
    //=============================================== //

    //=============================================== //
    //START COMPARE MS2
    //=============================================== //

    Fragment *observedSpectrum = nullptr;
    if (ms2Fragment && ms2Fragment->consensus) {
        observedSpectrum = ms2Fragment->consensus;
    }
    directInfusionMatchAssessment->computeMs2MatchAssessment(observedSpectrum, compound, params, debug);

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

       if (debug) cout << "DirectInfusionProcessor::getFragmentMatchMaps()" << endl;

       unique_ptr<DirectInfusionMatchInformation> matchInfo = unique_ptr<DirectInfusionMatchInformation>(new DirectInfusionMatchInformation());

       for (auto directInfusionMatchData : allCandidates) {

           Compound *compound = directInfusionMatchData->compound;
           FragmentationMatchScore fragmentationMatchScore = directInfusionMatchData->fragmentationMatchScore;

           vector<int> compoundFrags(static_cast<unsigned int>(fragmentationMatchScore.numMatches));

           unsigned int matchCounter = 0;
           for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

               int observedIndex = fragmentationMatchScore.ranks[i];

               //Issue 209: peaks may be unmatched based on intensity as well as ranks[] position
               if (!observedSpectrum || observedIndex == -1 || observedSpectrum->intensity_array[static_cast<unsigned long>(observedIndex)] < params->ms2MinIntensity) continue;

               if (debug) cout << "allCandidates[" << i << "]: " << compound->name << "|" << compound->adductString << " observedIndex=" << observedIndex << endl;

               long fragInt = mzToIntKey(compound->fragment_mzs[i]);

               compoundFrags[matchCounter] = static_cast<int>(fragInt);
               matchCounter++;

               pair<int, shared_ptr<DirectInfusionMatchData>> key = make_pair(fragInt, directInfusionMatchData);

               if (matchInfo->fragToMatchData.find(fragInt) == matchInfo->fragToMatchData.end()) {
                   matchInfo->fragToMatchData.insert(make_pair(fragInt, unordered_set<shared_ptr<DirectInfusionMatchData>>()));
               }

               matchInfo->fragToMatchData[fragInt].insert(directInfusionMatchData);

   //            if (debug) cerr << "allCandidates [end] i=" << i << ", ranks=" << fragmentationMatchScore.ranks.size() << endl;
           }

           if (compoundFrags.empty()) {
               matchInfo->compoundsNoFragMatches.push_back(directInfusionMatchData);
           } else {
               //Issue 270: need these to be sorted
               sort(compoundFrags.begin(), compoundFrags.end());
               matchInfo->matchDataToFrags.insert(make_pair(directInfusionMatchData, compoundFrags));
           }

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

        //Issue 314
        unordered_set<int> constituentMzs{};

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
        int ms1IntensityCoord = -1;
        float observedMs1ScanIntensity = 0.0f;

        for (auto matchData : compoundList) {

            compounds.push_back(matchData->compound);
            long mzIntKey = mzUtils::mzToIntKey(static_cast<double>(matchData->compound->precursorMz));
            constituentMzs.insert(static_cast<int>(mzIntKey));

            adduct = matchData->adduct;
            observedMs1Intensity += matchData->observedMs1Intensity;
            observedMs1ScanIntensity += matchData->observedMs1ScanIntensity;

            //Issue 288: intensity coordinate corresponds to the highest valid ms1 coord,
            //which will not always line up with the observed ms1 intensity.
            if (matchData->ms1IntensityCoord > ms1IntensityCoord){
                ms1IntensityCoord = matchData->ms1IntensityCoord;
            }

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
        observedMs1ScanIntensity /= compoundList.size();

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

            //Issue 305: include adduct information
            summarizedId += compounds.at(0)->adductString;

            summarizedCompound->adductString = compounds.at(0)->adductString;
            summarizedCompound->formula = compounds.at(0)->getFormula();
            summarizedCompound->precursorMz = compounds.at(0)->precursorMz;
            summarizedCompound->setExactMass(compounds.at(0)->getExactMass());
            summarizedCompound->charge = compounds.at(0)->charge;
            summarizedCompound->id = summarizedId;
            summarizedCompound->constituentMzsSet = constituentMzs;

            summarizedCompound->computeSummarizedData();

            summarizedMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());
            summarizedMatchData->compound = summarizedCompound;
            summarizedMatchData->adduct = adduct;

            if (observedSpectrum) {
                //Issue 313: Always use computeMs2MatchAssessment() method for computing ms2 match information
                unique_ptr<DirectInfusionMatchAssessment> directInfusionMatchAssessment = unique_ptr<DirectInfusionMatchAssessment>(new DirectInfusionMatchAssessment());
                directInfusionMatchAssessment->computeMs2MatchAssessment(observedSpectrum, summarizedCompound, params, debug);
                summarizedMatchData->fragmentationMatchScore = directInfusionMatchAssessment->fragmentationMatchScore;
            }

            summarizedMatchData->observedMs1Intensity = observedMs1Intensity;
            summarizedMatchData->ms1IntensityCoord = ms1IntensityCoord;
            summarizedMatchData->observedMs1ScanIntensity = observedMs1ScanIntensity;

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
            summarizedCompound->constituentMzsSet = constituentMzs;

            summarizedCompound->computeSummarizedData();

            summarizedMatchData = shared_ptr<DirectInfusionMatchData>(new DirectInfusionMatchData());
            summarizedMatchData->compound = summarizedCompound;
            summarizedMatchData->adduct = adduct;

            if (observedSpectrum) {
                //Issue 313: Always use computeMs2MatchAssessment() method for computing ms2 match information
                unique_ptr<DirectInfusionMatchAssessment> directInfusionMatchAssessment = unique_ptr<DirectInfusionMatchAssessment>(new DirectInfusionMatchAssessment());
                directInfusionMatchAssessment->computeMs2MatchAssessment(observedSpectrum, summarizedCompound, params, debug);
                summarizedMatchData->fragmentationMatchScore = directInfusionMatchAssessment->fragmentationMatchScore;
            }

            summarizedMatchData->observedMs1Intensity = observedMs1Intensity;
            summarizedMatchData->ms1IntensityCoord = ms1IntensityCoord;
            summarizedMatchData->observedMs1ScanIntensity = observedMs1ScanIntensity;
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

    //Issue 311: If number of unique fragments > 0, compounds with no frag matches must be discarded
    if (params->ms2MinNumUniqueMatches > 0) {
        matchInfo->compoundsNoFragMatches.clear();
    }

    return matchInfo;
}

unique_ptr<DirectInfusionMatchInformation> DirectInfusionProcessor::getMatchInformation(
        vector<shared_ptr<DirectInfusionMatchData>> allCandidates,
        Fragment *observedSpectrum,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    if (debug) cout << "DirectInfusionProcessor::getMatchInformation()" << endl;

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
        cout << "Fragments --> Compounds: (" << matchInfo->matchDataToFrags.size() << " passing compounds)" << endl;

        for (auto iterator = matchInfo->fragToMatchData.begin(); iterator != matchInfo->fragToMatchData.end(); ++iterator) {
            int frag = iterator->first;
            unordered_set<shared_ptr<DirectInfusionMatchData>> compounds = iterator->second;
            cout<< "frag= " << intKeyToMz(frag) << " m/z : ";
            for (auto matchData : compounds) {
                cout << matchData->compound->name << "|" << matchData->compound->adductString << " ";
            }
            cout << endl;
        }

        cout << "Compounds --> Fragments: (" << matchInfo->fragToMatchData.size() << " matched fragments)" << endl;

        for (auto iterator = matchInfo->matchDataToFrags.begin(); iterator != matchInfo->matchDataToFrags.end(); ++iterator) {

            shared_ptr<DirectInfusionMatchData> directInfusionMatchData = iterator->first;
            vector<int> frags = iterator->second;

            cout << "Compound= " << directInfusionMatchData->compound->name << "|" << directInfusionMatchData->compound->adductString << ": ";
            for (auto frag : frags){
                cout << intKeyToMz(frag) << " ";
            }
            cout << endl;
        }

        cout << "Compounds with no fragment matches:";
        if (matchInfo->compoundsNoFragMatches.empty()) {
            cout << " (none) ";
        }
        cout << endl;

        for (auto directInfusionMatchData : matchInfo->compoundsNoFragMatches) {
            cout << "Compound= " << directInfusionMatchData->compound->name << "|" << directInfusionMatchData->compound->adductString << endl;
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
        if (it->second && it->second->fragmentationPattern) delete(it->second->fragmentationPattern);
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
       groupMatchData->ms1IntensityCoord = matchData->ms1IntensityCoord;
       groupMatchData->observedMs1ScanIntensity = matchData->observedMs1ScanIntensity;

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

//returns -1 if no valid scans found.
float DirectInfusionUtils::findNormalizedIntensity(const vector<Scan*>& scans,
                                                   float queryMz,
                                                   float standardMz,
                                                   shared_ptr<DirectInfusionSearchParameters> params,
                                                   bool debug
                                                   ){
    vector<float> normalizedIntensities;

    map<int, vector<float>> normalizedIntensitiesByMassDiff{};
    vector<int> keys{};

    for (auto scan : scans) {

        float singleScanNormalizedIntensity = scan->findNormalizedIntensity(queryMz, standardMz, params->ms1PpmTolr, params->ms1MinScanIntensity);

        if (singleScanNormalizedIntensity < 0) continue;

        if (params->isPreferSmallestScanMassWindow) {
            int massDiff = static_cast<int>(round(scan->getMaxMz()-scan->getMinMz()));
            if (normalizedIntensitiesByMassDiff.find(massDiff) == normalizedIntensitiesByMassDiff.end()) {
                normalizedIntensitiesByMassDiff.insert(make_pair(massDiff, vector<float>()));
                keys.push_back(massDiff);
            }
            normalizedIntensitiesByMassDiff[massDiff].push_back(singleScanNormalizedIntensity);
        }

        if (debug) cout << "Scan #"
                        << scan->scannum << ", "
                        << scan->filterString
                        << ", minMz: " << scan->getMinMz()
                        << ", maxMz: " << scan->getMaxMz()
                        << ", mzWidth: " << (scan->getMaxMz()-scan->getMinMz())
                        << ": normalized intensity="
                        << to_string(singleScanNormalizedIntensity)
                        << endl;

        normalizedIntensities.push_back(singleScanNormalizedIntensity);
    }

    //Issue 481
    if (params->isPreferSmallestScanMassWindow) {
        sort(keys.begin(), keys.end());
        if (!keys.empty() && normalizedIntensitiesByMassDiff.find(keys[0]) != normalizedIntensitiesByMassDiff.end()) {

            if (debug) cout << "Using mzDiffKey "
                            << keys[0]
                            << endl;

            normalizedIntensities = normalizedIntensitiesByMassDiff[keys[0]];
        }
    }

    if (debug) cout << "Found " << normalizedIntensities.size() << " scans." << endl;

    if (normalizedIntensities.empty()) return -1.0f;

    if (params->consensusIntensityAgglomerationType == Fragment::Mean) {
        return accumulate(normalizedIntensities.begin(), normalizedIntensities.end(), 0.0f) / normalizedIntensities.size();
    } else if (params->consensusIntensityAgglomerationType == Fragment::Median) {
        return median(normalizedIntensities);
    } else { //unsupported type
        if (debug) cout << "Unsupported quant agglomeration type, DirectInfusionUtils::findNormalizedIntensity() returning -1.0f" << endl;
        return -1.0f;
    }
}

//returns -1.0f if not able to find scans
//expects MS1 scans to be sorted in increasing order by scan number
ScanQuantOutput DirectInfusionUtils::findNearestScanNormalizedIntensity(const vector<Scan*>& scans,
                                                              float queryMz,
                                                              float standardMz,
                                                              shared_ptr<DirectInfusionSearchParameters> params,
                                                              int scanWidthInDa,
                                                              bool debug){

    ScanQuantOutput scanQuantOutput;

    if (debug) cout << "DirectInfusionUtils::findNearestScanNormalizedIntensity()" << endl;

    vector<ScanIntensity> queryScans{};
    vector<ScanIntensity> standardScans{};

    for (auto scan : scans) {

        int scanWidth = static_cast<int>(round(scan->upperLimitMz - scan->lowerLimitMz));

        //only consider scans of a certain width (in Da), if argument provided.
        //otherwise, remove this scan from consideration.
        if (scanWidthInDa > 0 && scanWidth != scanWidthInDa) {
            if (debug) cout << "Scan #"
                            << scan->scannum << ", "
                            << scan->filterString
                            << " has scan width of " << scanWidth
                            << ", required scan width is " << scanWidthInDa
                            << ". skipping this scan and continuing."
                            << endl;
            continue;
        }

        float queryMzIntensityCandidate = scan->findClosestMzIntensity(queryMz, params->ms1PpmTolr);

        if (queryMzIntensityCandidate > 0.0f && queryMzIntensityCandidate >= params->ms1MinScanIntensity) {

            queryScans.push_back(ScanIntensity(scan, queryMzIntensityCandidate, ScanIntensityType::QUERY));

            if (debug) cout << "Scan #"
                            << scan->scannum << ", "
                            << scan->filterString
                            << ", minMz: " << scan->getMinMz()
                            << ", maxMz: " << scan->getMaxMz()
                            << ", mzWidth: " << (scan->getMaxMz()-scan->getMinMz())
                            << ": query m/z=" << queryMz
                            << ", query intensity="
                            << to_string(queryMzIntensityCandidate)
                            << endl;
        }

        float standardMzIntensityCandidate = scan->findClosestMzIntensity(standardMz, params->ms1PpmTolr);

        if (standardMzIntensityCandidate > 0.0f && standardMzIntensityCandidate >= params->ms1MinScanIntensity) {

            standardScans.push_back(ScanIntensity(scan, standardMzIntensityCandidate, ScanIntensityType::STANDARD));

            if (debug) cout << "Scan #"
                            << scan->scannum << ", "
                            << scan->filterString
                            << ", minMz: " << scan->getMinMz()
                            << ", maxMz: " << scan->getMaxMz()
                            << ", mzWidth: " << (scan->getMaxMz()-scan->getMinMz())
                            << ": standard m/z=" << standardMz
                            << ", standard intensity="
                            << to_string(standardMzIntensityCandidate)
                            << endl;
        }

    }

    vector<NearestScanIntensityPair> pairs = ScanIntensity::matchStandardScanIntensitiesToQueryScanIntensities(standardScans, queryScans, params, debug);

    vector<float> normalizedIntensities(pairs.size());

    for (unsigned int i = 0; i < pairs.size(); i++){
        normalizedIntensities[i] = pairs[i].getIntensity();
    }

    if (debug) cout << "Found " << normalizedIntensities.size() << " scans." << endl;

    if (normalizedIntensities.empty()) return scanQuantOutput;

    if (debug) {
        cout << " Intensities: ";
        for (unsigned int i = 0; i < normalizedIntensities.size(); i++) {
            cout << normalizedIntensities[i] << " ";
        }
        cout << endl;
    }

    float outputIntensity = -1.0f;

    float medianIntensity = median(normalizedIntensities);

    if (params->consensusIntensityAgglomerationType == Fragment::Mean) {
        outputIntensity = accumulate(normalizedIntensities.begin(), normalizedIntensities.end(), 0.0f) / normalizedIntensities.size();
    } else if (params->consensusIntensityAgglomerationType == Fragment::Median) {
        outputIntensity = medianIntensity;
    } else { //unsupported type
        if (debug) cout << "Unsupported quant agglomeration type, DirectInfusionUtils::findNormalizedIntensity() returning -1.0f" << endl;
        return scanQuantOutput;
    }

    scanQuantOutput.isValid = true;
    scanQuantOutput.intensity = outputIntensity;
    scanQuantOutput.numMeasurements = normalizedIntensities.size();

    vector<float> deviations(normalizedIntensities.size());
    for (unsigned int i = 0; i < normalizedIntensities.size(); i++) {
        deviations[i] = abs(normalizedIntensities[i]-medianIntensity);
    }

    scanQuantOutput.medianAbsoluteDeviation = median(deviations);
    scanQuantOutput.scanDiff = pairs[0].dist;
    scanQuantOutput.scanWidth = pairs[0].queryScan.scanWidth;

    return scanQuantOutput;
}

vector<NearestScanIntensityPair> ScanIntensity::matchStandardScanIntensitiesToQueryScanIntensities(vector<ScanIntensity> queryScans,
                                                                                                   vector<ScanIntensity> standardScans,
                                                                                                   shared_ptr<DirectInfusionSearchParameters> params,
                                                                                                   bool debug) {
    map<int, vector<ScanIntensity>> scansByWidth{};

    vector<int> allWidths{};

    for (auto queryScanIntensity : queryScans) {
        int width = queryScanIntensity.scanWidth;
        if (scansByWidth.find(width) == scansByWidth.end()) {
            scansByWidth.insert(make_pair(width, vector<ScanIntensity>{}));
            allWidths.push_back(width);
        }
        scansByWidth[width].push_back(queryScanIntensity);
    }

    for (auto standardScanIntensity : standardScans) {
        int width = standardScanIntensity.scanWidth;
        if (scansByWidth.find(width) == scansByWidth.end()) {
            scansByWidth.insert(make_pair(width, vector<ScanIntensity>{}));
            allWidths.push_back(width);
        }
        scansByWidth[width].push_back(standardScanIntensity);
    }

    sort(allWidths.begin(), allWidths.end());

    //pair<scan dist, scan width>
    map<pair<int, int>, vector<NearestScanIntensityPair>> validCandidatesMap{};
    vector<pair<int, int>> validCandidateKeys{};

    for (auto width : allWidths) {

        vector<NearestScanIntensityPair> pairs{};
        vector<ScanIntensity> allScans = scansByWidth[width];

        sort(allScans.begin(), allScans.end(), [](const ScanIntensity& lhs, const ScanIntensity& rhs){
            return lhs.scan->scannum < rhs.scan->scannum;
        });

        vector<ScanIntensityMatch> allMatches{};

        for (unsigned int i = 0; i < allScans.size(); i++) {

            ScanIntensity scanIntensity = allScans[i];

            if (scanIntensity.scanIntensityType == ScanIntensityType::QUERY) {

                ScanIntensityMatch scanIntensityMatch;

                scanIntensityMatch.queryScan = scanIntensity;

                int leftDiff = -1;
                int rightDiff = -1;

                //check left
                if (i >= 1 && allScans[i-1].scanIntensityType == ScanIntensityType::STANDARD) {
                    leftDiff = scanIntensity.scan->scannum - allScans[i-1].scan->scannum;
                }

                //check right
                if (i < allScans.size()-1 && allScans[i+1].scanIntensityType == ScanIntensityType::STANDARD) {
                    rightDiff = allScans[i+1].scan->scannum - scanIntensity.scan->scannum;
                }

                //case: only left diff valid
                if (leftDiff >= 0 && rightDiff < 0) {
                    scanIntensityMatch.standardScan = allScans[i-1];
                    scanIntensityMatch.dist = leftDiff;

                //case: only right diff valid
                } else if (leftDiff < 0 && rightDiff >= 0) {
                    scanIntensityMatch.standardScan = allScans[i+1];
                    scanIntensityMatch.dist = rightDiff;

                //case: both sides valid
                }  else if (leftDiff >= 0 && rightDiff >= 0) {
                    if (leftDiff <= rightDiff) {
                        scanIntensityMatch.standardScan = allScans[i-1];
                        scanIntensityMatch.dist = leftDiff;
                    } else {
                        scanIntensityMatch.standardScan = allScans[i+1];
                        scanIntensityMatch.dist = rightDiff;
                    }
                }

                if (scanIntensityMatch.standardScan.scan) {
                    allMatches.push_back(scanIntensityMatch);
                }
            }


        }

        if (debug) cout << "Identified " << allMatches.size() << " candidate matches for scans with width " << width << " Da." << endl;

        sort(allMatches.begin(), allMatches.end(), [](const ScanIntensityMatch& lhs, const ScanIntensityMatch& rhs){
            return lhs.dist < rhs.dist;
        });

        vector<int> queryScanNumsUsed{};
        vector<int> standardScanNumsUsed{};

        int previousDist = -1;
        int currentDist = -1;

        int numAtCurrentDist = 0;

        bool isWroteBestWidth = false;

        for (unsigned int i = 0; i < allMatches.size(); i++) {

            ScanIntensityMatch scanIntensityMatch = allMatches[i];
            currentDist = scanIntensityMatch.dist;

            int queryScanNum = scanIntensityMatch.queryScan.scan->scannum;
            int standardScanNum = scanIntensityMatch.standardScan.scan->scannum;

            bool isHasQueryScanNum = find(queryScanNumsUsed.begin(), queryScanNumsUsed.end(), queryScanNum) != queryScanNumsUsed.end();
            bool isHasStandardScanNum = find(standardScanNumsUsed.begin(), standardScanNumsUsed.end(), standardScanNum) != standardScanNumsUsed.end();

            if (!isHasQueryScanNum && !isHasStandardScanNum) {

                if ((previousDist == currentDist && previousDist != -1) || previousDist == -1) {
                    numAtCurrentDist++;
                }

                //a dist change - either exit here, or reset and continue
                if (previousDist != currentDist && previousDist != -1) {
                    if (numAtCurrentDist > 0 && numAtCurrentDist >= params->minNumScansNearestScanNormalizedIntensity) {

                        pair<int, int> key = make_pair(previousDist, width);
                        validCandidatesMap.insert(make_pair(key, pairs));
                        validCandidateKeys.push_back(key);
                        isWroteBestWidth = true;
                        break;

                    } else {

                        if (debug) {
                            cout << "numAtCurrentDist="
                                 << numAtCurrentDist
                                 << ", threshold="
                                 << params->minNumScansNearestScanNormalizedIntensity
                                 << endl;
                        }

                        //reset to prepare for new count
                        pairs.clear();
                        numAtCurrentDist = 1; // this entry starts a new count.
                    }

                }

                queryScanNumsUsed.push_back(queryScanNum);
                standardScanNumsUsed.push_back(standardScanNum);

                NearestScanIntensityPair nearestScanintensityPair(scanIntensityMatch.standardScan, scanIntensityMatch.queryScan);
                pairs.push_back(nearestScanintensityPair);

                if (debug) {
                    cout << "width=" << width
                         << ", dist=" << currentDist
                         << ", numAtCurrentDist= " << numAtCurrentDist
                         << ", threshold= " << params->minNumScansNearestScanNormalizedIntensity
                         << " "
                         << "Query Scan #" << nearestScanintensityPair.queryScan.scan->scannum
                         << " <--> "
                         <<  "Standard Scan #" << nearestScanintensityPair.standardScan.scan->scannum
                          << " intensity: " << nearestScanintensityPair.queryScan.intensity << "/" << nearestScanintensityPair.standardScan.intensity << " "
                          << nearestScanintensityPair.getIntensity() << endl;
                }

                previousDist = currentDist;
            }

        } // end for (allMatches)

        // handle case where all scans have the same diff (so not yet written to the map)
        if (!isWroteBestWidth && numAtCurrentDist > 0 && numAtCurrentDist >= params->minNumScansNearestScanNormalizedIntensity) {
            pair<int, int> key = make_pair(previousDist, width);
            validCandidatesMap.insert(make_pair(key, pairs));
            validCandidateKeys.push_back(key);
        }

    } // end for (allWidths)

    if (validCandidateKeys.empty()){
        if (debug) cout << "Found no valid candidate keys." << endl;
        return vector<NearestScanIntensityPair>{}; //invalid
    }

    sort(validCandidateKeys.begin(), validCandidateKeys.end(), [](const pair<int, int>& lhs, const pair<int, int>& rhs){
        if (lhs.first == rhs.first) {
            return lhs.second < rhs.second;
        } else {
            return lhs.first < rhs.first;
        }
    });

    if (debug) {
        for (auto key : validCandidateKeys) {
            cout << "valid key: " << key.first << " scan diff, " << key.second << " m/z window" << endl;
        }
    }

    return validCandidatesMap[validCandidateKeys[0]];
}


/**
 * @brief getCompounds
 *
 * method to compute current set of compounds, based on the current status of the
 * data within this struct
 *
 * This is intentionally left as a method instead of a field to avoid duplicated state
 * dependencies.
 *
 * Issue 311: Includes all compounds associated without any fragment matches, which are kept separate
 * from the rest of the matches.
 *
 * @return
 */
vector<shared_ptr<DirectInfusionMatchData>> DirectInfusionMatchInformation::getCompounds(){
    vector<shared_ptr<DirectInfusionMatchData>> compounds(matchDataToFrags.size() + compoundsNoFragMatches.size());

    unsigned int i = 0;
    for (auto it = matchDataToFrags.begin(); it != matchDataToFrags.end(); ++it){
        compounds[i] = it->first;
        i++;
    }

    for (auto compound : compoundsNoFragMatches) {
        compounds[i] = compound;
        i++;
    }

    return compounds;
}

//WARNING: not rounded to nearest precision, just output that many digits
string DirectInfusionMatchInformation::getFragmentGroupId(shared_ptr<DirectInfusionMatchData> compound, int precision){

    //guard to avoid nonsensical output
    if (precision < 0) {
        precision = 0;
    } else if (precision > 6) {
        precision = 6;
    }

    if (matchDataToFrags.find(compound) != matchDataToFrags.end()) {
        vector<int> frags = matchDataToFrags[compound];

        stringstream s;
        s << std::fixed << setprecision(precision);
        s << "(";

        for (unsigned int i = 0; i < frags.size(); i++) {

            if (i > 0) {
                s << ", ";
            }

            s << mzUtils::intKeyToMz(frags[i]);
        }

        s << ")";

        return s.str();
    } else {
        return "";
    }
}

void DirectInfusionMatchInformation::computeMs1PartitionFractions2(
        const Fragment *ms2Fragment,
        const shared_ptr<DirectInfusionSearchParameters> params,
        const bool debug){

    if (debug) cout << "DirectInfusionMatchInformation::computeMs1PartitionFractions2()" << endl;

    if (!ms2Fragment || !ms2Fragment->consensus) return;

    map<shared_ptr<DirectInfusionMatchData>, float> partitionFractions{};

    // <compoundName, fragId> = division of this fragment between multiple compounds in same window
    map<pair<shared_ptr<DirectInfusionMatchData>, int>, float> fragToSAFMultiplier;
    set<shared_ptr<DirectInfusionMatchData>> compoundsWithAdjustedSAFs{};

    map<long, vector<shared_ptr<DirectInfusionMatchData>>> ms1MzToCompoundNames{};

    for (auto compound : getCompounds()) {
        long coord = mzUtils::mzToIntKey(static_cast<double>(compound->observedMs1ScanIntensity), 1L);

        //If the observed intensity is 0, there is nothing to partition
        if (coord == 0L) continue;

        if (ms1MzToCompoundNames.find(coord) == ms1MzToCompoundNames.end()) {
            ms1MzToCompoundNames.insert(make_pair(coord, vector<shared_ptr<DirectInfusionMatchData>>{}));
        }
        ms1MzToCompoundNames[coord].push_back(compound);
    }

    //STEP 1: determine SAF adjustments

    for (auto it = fragToMatchData.begin(); it != fragToMatchData.end(); ++it) {

        int fragId = it->first;
        unordered_set<shared_ptr<DirectInfusionMatchData>> compounds = it->second;

        set<long> ms1Ids{};

        if (compounds.size() > 1) {

            for (auto compound : compounds) {
                for (auto it2 = ms1MzToCompoundNames.begin(); it2 != ms1MzToCompoundNames.end(); ++it2) {
                    long ms1Id = it2->first;
                    vector<shared_ptr<DirectInfusionMatchData>> ms1Compounds = it2->second;

                    if (ms1Id == 0L) continue;

                    if (std::find(ms1Compounds.begin(), ms1Compounds.end(), compound) != ms1Compounds.end()){
                        ms1Ids.insert(ms1Id);
                    }
                }
            }

            //SAF is only a concern when a fragment is involved in multiple compounds and multiple ms1 m/zs
            if (ms1Ids.size() > 1) {

                float totalIntensity = 0.0f;
                for (auto ms1Id : ms1Ids) {
                    totalIntensity += static_cast<float>(mzUtils::intKeyToMz(ms1Id, 1L));
                }

                for (auto ms1Id : ms1Ids) {

                    float ms1IdIntensity = static_cast<float>(mzUtils::intKeyToMz(ms1Id, 1L));
                    vector<shared_ptr<DirectInfusionMatchData>> compoundsWithMatchingId = ms1MzToCompoundNames.at(ms1Id);

                    for (auto matchData : compoundsWithMatchingId) {

                        //even if a compound doesn't directly have an ambigous fragment,
                        //if it shares an MS1 intensity peak with a compound that
                        //does have an ambiguous fragment, it is also affected by the SAF
                        compoundsWithAdjustedSAFs.insert(matchData);

                        if (debug) cout << matchData->compound->name << " "
                                        << matchData->compound->adductString << " "
                                        << "has ambiguous fragment m/z="
                                        << mzUtils::intKeyToMz(fragId)
                                        << endl;

                        pair<shared_ptr<DirectInfusionMatchData>, long> key = make_pair(matchData, fragId);

                        if (fragToSAFMultiplier.find(key) == fragToSAFMultiplier.end()) {
                            float SAFpartition = ms1IdIntensity/totalIntensity;
                            fragToSAFMultiplier.insert(make_pair(key, SAFpartition));
                        }

                    }
                }

            //multiple compounds, but single ms1 m/z case: no need for SAF
            } else {
                for (auto compound : compounds) {
                    pair<shared_ptr<DirectInfusionMatchData>, long> key = make_pair(compound, fragId);
                    if (fragToSAFMultiplier.find(key) == fragToSAFMultiplier.end()) {
                        fragToSAFMultiplier.insert(make_pair(key, 1.0f));
                    }
                }
            }

        } else {
            pair<shared_ptr<DirectInfusionMatchData>, long> key = make_pair(*(compounds.begin()), fragId);
            if (fragToSAFMultiplier.find(key) == fragToSAFMultiplier.end()) {
                fragToSAFMultiplier.insert(make_pair(key, 1.0f));
            }
        }
    }

    //STEP 2: divide ms1 intensities between all compounds, based on fragments

    for (auto it = ms1MzToCompoundNames.begin(); it != ms1MzToCompoundNames.end(); ++it) {

        vector<shared_ptr<DirectInfusionMatchData>> compoundNames = it->second;

        float totalEfectiveFragIntensityAllCompounds = 0.0f;

        map<shared_ptr<DirectInfusionMatchData>, float> compoundToTotalEffectiveFragIntensity{};

        for (auto matchData : compoundNames) {

            if (debug) {
                cout << "ms1 intensity: " << it->first
                     << ", compound: " << matchData->compound->name
                     << ", adduct: " << matchData->compound->adductString
                     << endl;
            }

            //not all IDs have MS2 fragment matches
            vector<int> fragMzs{};
            map<int, float> fragMzToIntensity{};

            if (matchDataToFrags.find(matchData) != matchDataToFrags.end()) {

                fragMzs = matchDataToFrags.at(matchData);

                vector<int> ranks = matchData->fragmentationMatchScore.ranks;

                for (unsigned int i = 0; i < ranks.size(); i++) {

                    int y = ranks[i];
                    if (y == -1) continue;

                    float fragObservedIntensity = ms2Fragment->consensus->intensity_array[static_cast<unsigned int>(y)];

                    if (fragObservedIntensity >= params->ms2MinIntensity) {

                        string fragmentLabel = matchData->compound->fragment_labels[i];

                        vector<string> fragmentLabelTags = DirectInfusionMatchAssessment::getFragmentLabelTags(fragmentLabel, params, false);

                        if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn1FragmentLabelTag") != fragmentLabelTags.end() ||
                            find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn2FragmentLabelTag") != fragmentLabelTags.end()) {

                            int mzKey = static_cast<int>(mzUtils::mzToIntKey(static_cast<double>(matchData->compound->fragment_mzs[i])));

                            if (debug) cout << "Acyl chain fragment: " << fragmentLabel
                                            <<   ", m/z=" << matchData->compound->fragment_mzs[i]
                                            << ", intensity=" << fragObservedIntensity
                                            << endl;

                            fragMzToIntensity.insert(make_pair(mzKey, fragObservedIntensity));

                        }
                    }
                }
            }

            float effectiveFragIntensityCompound = 0.0f;

            for (auto fragId : fragMzs) {

                if (fragMzToIntensity.find(fragId) == fragMzToIntensity.end()) continue;

                pair<shared_ptr<DirectInfusionMatchData>, int> safKey = make_pair(matchData, fragId);

                if (fragToSAFMultiplier.find(safKey) == fragToSAFMultiplier.end()) continue;

                float effectiveFragIntensity = (fragMzToIntensity.at(fragId) * fragToSAFMultiplier.at(safKey));

                effectiveFragIntensityCompound += effectiveFragIntensity;
                totalEfectiveFragIntensityAllCompounds += effectiveFragIntensity;

            }

            compoundToTotalEffectiveFragIntensity.insert(make_pair(matchData, effectiveFragIntensityCompound));

        }

        for (auto it2 = compoundToTotalEffectiveFragIntensity.begin(); it2 != compoundToTotalEffectiveFragIntensity.end(); ++it2) {

            auto matchData = it2->first;
            float effectiveFragIntensityCompound = it2->second;

            float partitionFraction = 1.0f/compoundNames.size();
            if (totalEfectiveFragIntensityAllCompounds > 0.0f) {
                partitionFraction = effectiveFragIntensityCompound/totalEfectiveFragIntensityAllCompounds;
            }

            partitionFractions.insert(make_pair(matchData, partitionFraction));
        }

    }

    //STEP 3: disqualify SAFs, if appropriate

    for (auto it = partitionFractions.begin(); it != partitionFractions.end(); ++it) {

        auto matchData = it->first;
        matchData->ms1PartitionFractionSplitAmbiguousFragments = it->second;

        if (std::find(compoundsWithAdjustedSAFs.begin(), compoundsWithAdjustedSAFs.end(), matchData) != compoundsWithAdjustedSAFs.end()) {
            matchData->ms1PartitionFraction = -1.0f;
        } else {
            matchData->ms1PartitionFraction = it->second;
        }

        //TODO: handle scans, or delete these fields
        matchData->ms1PartitionFractionByScan = -1.0f;
        matchData->ms1PartitionFractionByScanSplitAmbiguousFragments = -1.0f;

        if (debug) cout << matchData->compound->name << " "
                        << matchData->adduct->name << ":"
                        << " partitionFraction= " << matchData->ms1PartitionFraction
                        << ", partitionFractionSAF= " << matchData->ms1PartitionFractionSplitAmbiguousFragments
                        << endl;
    }
}

//Issue 288, 292
void DirectInfusionMatchInformation::computeMs1PartitionFractions(const vector<Scan*>& ms2Scans,
                                                                  const Fragment *ms2Fragment,
                                                                  const shared_ptr<DirectInfusionSearchParameters> params,
                                                                  const bool debug) {

    if (debug) cout << "DirectInfusionMatchInformation::computeMs1PartitionFractions()" << endl;

    if (!ms2Fragment || !ms2Fragment->consensus) return;
    if (ms2Scans.empty()) return;

    //Issue 318
    map<int, float> fragMzToSumObservedMs1ScanIntensity= getFragToSumObservedMs1ScanIntensity(debug);

    map<int, vector<shared_ptr<DirectInfusionMatchData>>> partitionMap = getPrecMzPartitionMap(debug);

    map<int, set<int>> partitionMapMzs{};

    for (auto it = partitionMap.begin(); it != partitionMap.end(); ++it){

        float allFragIntensity = 0.0f;
        float allFragIntensitySAF = 0.0f;

        map<shared_ptr<DirectInfusionMatchData>, float> totalFragIntensityByCompound{};
        map<shared_ptr<DirectInfusionMatchData>, float> totalFragIntensitySAFByCompound{};

        //Issue 292
        map<Scan*, float> totalFragIntensityByScan{};
        map<Scan*, map<shared_ptr<DirectInfusionMatchData>, float>> compoundFragIntensityByScan{};

        //Issue 318
        map<Scan*, float> totalFragIntensitySAFByScan{};
        map<Scan*, map<shared_ptr<DirectInfusionMatchData>, float>> compoundFragIntensitySAFByScan{};

        for (auto matchData : it->second) {

            if (debug) {
                cout << "compound: " << matchData->compound->name << ", adduct: " << matchData->compound->adductString << endl;
            }

            float compoundFragIntensity = 0.0f;
            float compoundFragIntensitySAF = 0.0f;

            vector<int> ranks = matchData->fragmentationMatchScore.ranks;

            for (unsigned int i = 0; i < ranks.size(); i++) {
                int y = ranks[i];

                if (y == -1) continue;

                float fragObservedIntensity = ms2Fragment->consensus->intensity_array[y];

                if (fragObservedIntensity >= params->ms2MinIntensity) {

                    string fragmentLabel = matchData->compound->fragment_labels[i];

                    //Issue 388: Disable debugging for easier readout
                    vector<string> fragmentLabelTags = DirectInfusionMatchAssessment::getFragmentLabelTags(fragmentLabel, params, false);

                    bool isAcylFragment = false;

                    //Issue 313: use marked sn1/sn2 fragment label tags instead of explicit list of fragment labels
                    if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn1FragmentLabelTag") != fragmentLabelTags.end() ||
                        find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn2FragmentLabelTag") != fragmentLabelTags.end()) {

                        isAcylFragment = true;
                        compoundFragIntensity += fragObservedIntensity;

                        //Issue 318: split ambiguous fragments to avoid multiple-counting errors
                        int fragKey = static_cast<int>(mzUtils::mzToIntKey(static_cast<double>(matchData->compound->fragment_mzs[i])));
                        float sumObservedMs1ScanIntensity = fragMzToSumObservedMs1ScanIntensity[fragKey];

                        compoundFragIntensitySAF += (fragObservedIntensity * (matchData->observedMs1ScanIntensity/sumObservedMs1ScanIntensity));

                        if (partitionMapMzs.find(it->first) == partitionMapMzs.end()) {
                            partitionMapMzs.insert(make_pair(it->first, set<int>()));
                        }

                        long key = mzUtils::mzToIntKey(static_cast<double>(matchData->compound->fragment_mzs[i]));
                        partitionMapMzs[it->first].insert(static_cast<int>(key));
                    }

                    //Issue 388: print info about each fragment
                    if (debug) {
                        cout << (isAcylFragment ? "ACYL FRAGMENT " : " NON-ACYL FRAGMENT ")
                             << "fragment label: " << fragmentLabel
                             << ", fragment intensity=" << fragObservedIntensity
                             << ", compoundFragIntensity=" << compoundFragIntensity
                             << ", compoundFragIntensitySAF=" << compoundFragIntensitySAF
                             << ", allFragIntensity=" << allFragIntensity
                             << ", allFragIntensitySAF=" << allFragIntensitySAF
                             << endl;
                    }

                }
            }

            totalFragIntensityByCompound.insert(make_pair(matchData, compoundFragIntensity));
            allFragIntensity += compoundFragIntensity;

            totalFragIntensitySAFByCompound.insert(make_pair(matchData, compoundFragIntensitySAF));
            allFragIntensitySAF += compoundFragIntensitySAF;

            //SCAN-BASED QUANT APPROACH

            map<int, vector<float>> partitionFractionsByNumFragmentsPresent{};

            vector<float> partitionFragmentMzs;
            for (unsigned int i = 0; i < matchData->compound->fragment_labels.size(); i++) {

                string fragmentLabel = matchData->compound->fragment_labels[i];
                float fragmentMz = matchData->compound->fragment_mzs[i];

                //Issue 313: use marked sn1/sn2 fragment label tags instead of explicit list of fragment labels
                vector<string> fragmentLabelTags = DirectInfusionMatchAssessment::getFragmentLabelTags(fragmentLabel, params, false);

                if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn1FragmentLabelTag") != fragmentLabelTags.end() ||

                    find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn2FragmentLabelTag") != fragmentLabelTags.end()) {

                    partitionFragmentMzs.push_back(fragmentMz);

                    if (debug) {
                        cout << "fragment label: " << fragmentLabel
                             << ", fragmentMz: " << fragmentMz
                             << endl;
                    }
                }
            }

            if (debug) {
                cout << matchData->compound->name << " partitionMzs: ";
                for (auto partitionMz : partitionFragmentMzs) {
                    cout << partitionMz << " ";
                }
                cout << endl;
            }

            for (auto scan : ms2Scans) {

                if (debug) cout << "scan #" << scan->scannum << ": ";

                float scanSumIntensity = 0.0f;
                float scanSumSAFIntensity = 0.0f;

                for (auto queryMz : partitionFragmentMzs) {
                    if (debug) cout << "mz=" << queryMz << ": ";
                    float queryIntensity = scan->findClosestMzIntensity(queryMz, params->ms2PpmTolr);
                    if (queryIntensity > 0) {

                        long fragKey = mzUtils::mzToIntKey(static_cast<double>(queryMz));
                        float sumObservedMs1ScanIntensity = fragMzToSumObservedMs1ScanIntensity[static_cast<int>(fragKey)];

                        float queryIntensitySAF = queryIntensity * (matchData->observedMs1ScanIntensity/sumObservedMs1ScanIntensity);

                        scanSumIntensity += queryIntensity;
                        scanSumSAFIntensity += queryIntensitySAF;

                        if (debug) cout << "intensity=" << queryIntensity << ", SAF intensity=" << queryIntensitySAF << " ";
                    } else {
                        if (debug) cout << "NA ";
                    }
                }

                //Issue 292
                if (totalFragIntensityByScan.find(scan) == totalFragIntensityByScan.end()) {
                    totalFragIntensityByScan.insert(make_pair(scan, 0.0f));
                }

                if (scanSumIntensity > 0.0f) {
                    totalFragIntensityByScan[scan] += scanSumIntensity;

                    if (compoundFragIntensityByScan.find(scan) == compoundFragIntensityByScan.end()) {
                        compoundFragIntensityByScan.insert(make_pair(scan, map<shared_ptr<DirectInfusionMatchData>, float>()));
                    }
                    compoundFragIntensityByScan[scan].insert(make_pair(matchData, scanSumIntensity));

                    if (debug) cout << ", compound total: " << scanSumIntensity << ", all total: " << totalFragIntensityByScan[scan] << endl;
                }

                //Issue 318
                if (totalFragIntensitySAFByScan.find(scan) == totalFragIntensitySAFByScan.end()) {
                    totalFragIntensitySAFByScan.insert(make_pair(scan, 0.0f));
                }

                if (scanSumSAFIntensity > 0.0f) {
                    totalFragIntensitySAFByScan[scan] += scanSumSAFIntensity;

                    if (compoundFragIntensitySAFByScan.find(scan) == compoundFragIntensitySAFByScan.end()) {
                        compoundFragIntensitySAFByScan.insert(make_pair(scan, map<shared_ptr<DirectInfusionMatchData>, float>()));
                    }
                    compoundFragIntensitySAFByScan[scan].insert(make_pair(matchData, scanSumSAFIntensity));

                    if (debug) cout << ", SAF compound total: " << scanSumSAFIntensity << ", all SAF total: " << totalFragIntensitySAFByScan[scan] << endl;
                }

            }

        }

        computeMs1PartitionFractionFromMap(totalFragIntensityByCompound, allFragIntensity, false, debug);
        computeMs1PartitionFractionFromMap(totalFragIntensitySAFByCompound, allFragIntensitySAF, true, debug);

//        //Issue 292
//        map<shared_ptr<DirectInfusionMatchData>, vector<float>> scanPartitionFractions{};

//        for (auto it2 = compoundFragIntensityByScan.begin(); it2 != compoundFragIntensityByScan.end(); ++it2){

//            Scan* scan = it2->first;
//            map<shared_ptr<DirectInfusionMatchData>, float> fragIntensityByCompound = it2->second;
//            float totalScanIntensity = totalFragIntensityByScan[scan];

//            if (totalScanIntensity > 0.0f) {
//                for (auto it2 = fragIntensityByCompound.begin(); it2 != fragIntensityByCompound.end(); ++it2) {

//                    shared_ptr<DirectInfusionMatchData> matchData = it2->first;
//                    float compoundTotalIntensity = it2->second;

//                    float scanPartitionFraction = compoundTotalIntensity / totalScanIntensity;

//                    if (scanPartitionFractions.find(matchData) == scanPartitionFractions.end()) {
//                        scanPartitionFractions.insert(make_pair(matchData, vector<float>()));
//                    }
//                    scanPartitionFractions[matchData].push_back(scanPartitionFraction);

//                    //Issue 292
//                    if (debug) {
//                        cout << "compound: " << matchData->compound->name
//                             << ", scan #" << scan->scannum
//                             << ", compoundFragIntensity/allFragIntensity = " << compoundTotalIntensity << "/" << totalScanIntensity
//                             << ", fraction = " << scanPartitionFraction
//                             << endl;
//                    }
//                }
//            }

//        } // END for (auto it2 = compoundFragIntensityByScan.begin(); it2 != compoundFragIntensityByScan.end(); ++it2)

//        for (auto it2 = scanPartitionFractions.begin(); it2 != scanPartitionFractions.end(); ++it2) {

//            vector<float> intensities = it2->second;

//            if (debug) {
//                cout << "compound: " << it2->first->compound->name
//                     << ", scan partition value ";
//            }

//            if (!intensities.empty()) {
//                if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
//                    it2->first->ms1PartitionFractionByScan = accumulate(intensities.begin(), intensities.end(), 0.0f) / intensities.size();

//                    if (debug) cout << "[mean] = ";

//                } else if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
//                    it2->first->ms1PartitionFractionByScan = median(intensities);

//                    if (debug) cout << "[median] = ";
//                }
//            }

//            if (debug) cout << it2->first->ms1PartitionFractionByScan << endl;

//        } // END for (auto it2 = scanPartitionFractions.begin(); it2 != scanPartitionFractions.end(); ++it2)

        computeScanMs1PartitionFractionFromMap(totalFragIntensityByScan, compoundFragIntensityByScan, false, params, debug);
        computeScanMs1PartitionFractionFromMap(totalFragIntensitySAFByScan, compoundFragIntensitySAFByScan, true, params, debug);

        //If only one compound exists in the IDed compounds group, it gets all of the partitioned intensity.
        if (it->second.size() == 1) {
            it->second[0]->ms1PartitionFraction = 1;
            it->second[0]->ms1PartitionFractionByScan = 1;
            it->second[0]->ms1PartitionFractionSplitAmbiguousFragments = 1;
            it->second[0]->ms1PartitionFractionByScanSplitAmbiguousFragments = 1;
        }

    } // END for (auto it = partitionMap.begin(); it != partitionMap.end(); ++it)

    //Issue 314: If any m/zs are shared between sets,
    //all compounds associated with the shared m/z are invalidated -
    //that is, the partition fraction cannot be accurately determined.
    //in that case, the partition fraction is returned as -1 (indicating an uncomputable value).
    //
    //Additionally, if any set is associated with multiple distinct m/zs, the partition fraction
    //for that set cannot be determined, and is returned as -1 (indicating an uncomputable value).
    //
    //Any set that shares a fragment with such a multiple distinct m/z set is also invalidated.

    //reflects the true m/z <==> compounds mapping, prior to mutation by summarization.
    map<int, vector<shared_ptr<DirectInfusionMatchData>>> constituentMzExpandedPartitionMap{};

    //disqualify multiple sets, add singletons to map
    for (auto it = partitionMap.begin(); it != partitionMap.end(); ++it){

        vector<shared_ptr<DirectInfusionMatchData>> matchDataVector = it->second;

        for (auto matchData : matchDataVector) {
            vector<int> constituentMzs = matchData->compound->getConstituentMzs();

            if (constituentMzs.size() > 1) {
                //invalidate this entry
                matchData->ms1PartitionFraction = -1;
                matchData->ms1PartitionFractionByScan = -1;
            }

            //re-organize compounds based on all constituent m/z values.
            for (auto consituentMz : constituentMzs) {
                if (constituentMzExpandedPartitionMap.find(consituentMz) == constituentMzExpandedPartitionMap.end()) {
                    constituentMzExpandedPartitionMap.insert(make_pair(consituentMz, vector<shared_ptr<DirectInfusionMatchData>>{}));
                }
                constituentMzExpandedPartitionMap[consituentMz].push_back(matchData);
            }
        }
    }

    if (partitionMapMzs.size() > 1) {

        set<int> invalidatedSets{};

        for (auto it = partitionMapMzs.begin(); it != partitionMapMzs.end(); ++it) {

            int keyI = it->first;
            set<int> ithSet = it->second;

            for (auto it2 = partitionMapMzs.begin(); it2 != partitionMapMzs.end(); ++it2) {

                int keyJ = it2->first;
                set<int> jthSet = it2->second;

                if (keyI == keyJ) continue;

                vector<int> intersectionVector{};
                set_intersection(ithSet.begin(), ithSet.end(), jthSet.begin(), jthSet.end(), back_inserter(intersectionVector));
                if (!intersectionVector.empty()) {
                    invalidatedSets.insert(static_cast<int>(keyI));
                    invalidatedSets.insert(static_cast<int>(keyJ));
                }

            }
        }

        //invalidate, when appropriate
        for (auto it = constituentMzExpandedPartitionMap.begin(); it != constituentMzExpandedPartitionMap.end(); ++it){

            if (invalidatedSets.find(it->first) != invalidatedSets.end()) {

                //Issue 314: Groups should only be invalidated because of overlapping m/zs
                //when the group itself contains multiple candidate compounds.
                //Singleton groups don't require any partitioning, so ambiguity in fragment m/zs has no effect
                if (it->second.size() > 1) {

                    for (auto matchData : it->second) {
                        matchData->ms1PartitionFraction = -1;
                        matchData->ms1PartitionFractionByScan = -1;

                        if (debug) {
                            cout << matchData->compound->id
                                 << " ms1PartitionFraction and ms1PartitionFractionByScan were set to -1, indicating that a partition fraction could not be determined."
                                 << endl;
                        }

                    }

                }

            }
        }

    }
}

/**
 * @brief DirectInfusionMatchInformation::getFragToSumObservedMs1ScanIntensity
 *
 * Given a fragment m/z, return the total of all observedMs1ScanIntensity values
 * associated with all compound matches that contain this fragment m/z.
 *
 * However, to avoid multiple-counting observedMs1ScanIntensity values, only consider
 * a compound with a given precursor m/z value one time.  At high tolerance, this should
 * correspond to different observedMs1ScanIntensity values.
 *
 * The refMs1Mz value is used instead of the observedMs1ScanIntensity value directly
 * to handle general summarized compounds, that may contain more than one constituent m/z.
 *
 * Note also that if two distinct refMs1Mz values map the same observedMs1ScanIntensity value,
 * the downstream consequences will be a 50/50 split of the fragment intensity (which is
 * a fair estimate given no other information about the compounds or fragments).
 *
 * These values should only be computing using real precursor m/z values, not synthetic averages
 * from summarization.
 *
 * @param debug
 * @return
 */
map<int, float> DirectInfusionMatchInformation::getFragToSumObservedMs1ScanIntensity(
        const bool debug){

    map<pair<int, int>, float> precFragToSumObservedMs1ScanIntensity{};

    map<int, float> fragMzToSumObservedMs1ScanIntensity{};
    map<int, unordered_set<int>> fragMzToPrecMzs{};

    for (auto matchData : getCompounds()){

            //split based on real, observed prec m/zs, not on synthetic averaged m/zs
            vector<int> realPrecMzs = matchData->compound->getConstituentMzs();

            for (auto fragMzVal : matchData->compound->fragment_mzs) {

                int fragMz = static_cast<int>(mzUtils::mzToIntKey(static_cast<double>(fragMzVal)));

                for (int precMz : realPrecMzs) {

                    bool isAddIntensity = true;

                    if (fragMzToPrecMzs.find(fragMz) != fragMzToPrecMzs.end()) {
                         unordered_set<int> precMzs = fragMzToPrecMzs[fragMz];
                         if (precMzs.find(precMz) != precMzs.end()) {
                             isAddIntensity = false;
                         } else {
                             fragMzToPrecMzs[precMz].insert(precMz);
                         }
                    } else {
                        fragMzToPrecMzs.insert(make_pair(fragMz, unordered_set<int>{precMz}));
                    }

                    if (isAddIntensity) {
                        if (fragMzToSumObservedMs1ScanIntensity.find(fragMz) == fragMzToSumObservedMs1ScanIntensity.end()) {
                            fragMzToSumObservedMs1ScanIntensity.insert(make_pair(fragMz, 0.0f));
                        }
                        fragMzToSumObservedMs1ScanIntensity[fragMz] += matchData->observedMs1ScanIntensity;
                    }
            }

        }
    }

    if (debug) {
        cout << "DirectInfusionMatchInformation::getFragToSumObservedMs1ScanIntensity(): " << endl;
        for (auto it = fragMzToSumObservedMs1ScanIntensity.begin(); it != fragMzToSumObservedMs1ScanIntensity.end(); ++it) {
            double fragMz = mzUtils::intKeyToMz(it->first);
            cout << "m/z=" << fragMz << ", intensity=" << it->second << endl;
        }
    }

    return fragMzToSumObservedMs1ScanIntensity;
}

void DirectInfusionMatchInformation::computeMs1PartitionFractionFromMap(
        const map<shared_ptr<DirectInfusionMatchData>, float>& totalFragIntensityByCompound,
        const float allFragIntensity,
        bool isSAF,
        const bool debug){

    for (auto it2 = totalFragIntensityByCompound.begin(); it2 != totalFragIntensityByCompound.end(); ++it2) {

        float compoundFragIntensity = it2->second;
        if (allFragIntensity > 0.0f) {

            float partitionFraction = compoundFragIntensity / allFragIntensity;

            if (isSAF) {
                it2->first->ms1PartitionFractionSplitAmbiguousFragments = partitionFraction;
            } else {
                it2->first->ms1PartitionFraction = partitionFraction;
            }

            if (debug) {
                cout << "compound: " << it2->first->compound->name
                     << ", adduct: " << it2->first->compound->adductString
                     << ", compoundFragIntensity/allFragIntensity = " << compoundFragIntensity << "/" << allFragIntensity
                     << " = " << partitionFraction
                     << endl;
            }
        }
    }
}

void DirectInfusionMatchInformation::computeScanMs1PartitionFractionFromMap(
        const map<Scan*, float>& totalFragIntensityByScan,
        const map<Scan*, map<shared_ptr<DirectInfusionMatchData>, float>>& compoundFragIntensityByScan,
        bool isSAF,
        const shared_ptr<DirectInfusionSearchParameters> params,
        const bool debug){

    //Issue 292
    map<shared_ptr<DirectInfusionMatchData>, vector<float>> scanPartitionFractions{};

    for (auto it2 = compoundFragIntensityByScan.begin(); it2 != compoundFragIntensityByScan.end(); ++it2){

        Scan* scan = it2->first;
        map<shared_ptr<DirectInfusionMatchData>, float> fragIntensityByCompound = it2->second;
        float totalScanIntensity = totalFragIntensityByScan.at(scan);

        if (totalScanIntensity > 0.0f) {
            for (auto it2 = fragIntensityByCompound.begin(); it2 != fragIntensityByCompound.end(); ++it2) {

                shared_ptr<DirectInfusionMatchData> matchData = it2->first;
                float compoundTotalIntensity = it2->second;

                float scanPartitionFraction = compoundTotalIntensity / totalScanIntensity;

                if (scanPartitionFractions.find(matchData) == scanPartitionFractions.end()) {
                    scanPartitionFractions.insert(make_pair(matchData, vector<float>()));
                }
                scanPartitionFractions[matchData].push_back(scanPartitionFraction);

                //Issue 292
                if (debug) {
                    cout << "compound: " << matchData->compound->name
                         << ", scan #" << scan->scannum
                         << ", compoundFragIntensity/allFragIntensity = " << compoundTotalIntensity << "/" << totalScanIntensity
                         << ", fraction = " << scanPartitionFraction
                         << endl;
                }
            }
        }

    } // END for (auto it2 = compoundFragIntensityByScan.begin(); it2 != compoundFragIntensityByScan.end(); ++it2)

    for (auto it2 = scanPartitionFractions.begin(); it2 != scanPartitionFractions.end(); ++it2) {

        vector<float> intensities = it2->second;

        if (debug) {
            cout << "compound: " << it2->first->compound->name
                 << ", scan partition value ";
        }

        if (!intensities.empty()) {
            if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {

                if (isSAF) {
                    it2->first->ms1PartitionFractionByScanSplitAmbiguousFragments = accumulate(intensities.begin(), intensities.end(), 0.0f) / intensities.size();
                } else {
                    it2->first->ms1PartitionFractionByScan = accumulate(intensities.begin(), intensities.end(), 0.0f) / intensities.size();
                }

                if (debug) cout << "[mean] = ";

            } else if (params->consensusIntensityAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
                if (isSAF) {
                    it2->first->ms1PartitionFractionByScanSplitAmbiguousFragments = median(intensities);
                } else {
                    it2->first->ms1PartitionFractionByScan = median(intensities);
                }

                if (debug) cout << "[median] = ";
            }
        }

        if (debug) cout << it2->first->ms1PartitionFractionByScan << endl;

    } // END for (auto it2 = scanPartitionFractions.begin(); it2 != scanPartitionFractions.end(); ++it2)

}

map<int, vector<shared_ptr<DirectInfusionMatchData>>> DirectInfusionMatchInformation::getPrecMzPartitionMap(const bool debug){

    map<int, vector<shared_ptr<DirectInfusionMatchData>>> partitionMap{};

    for (auto it = matchDataToFrags.begin(); it != matchDataToFrags.end(); ++it){

        //Issue 314: Prefer theoretical m/z as key over observed intensity from consensus ms1
        //This allows fractions to be computed when ms1 intensity can be identified in ms1 scans,
        //but not in the consensus ms1 spectrum.
        int key = -1;
        if (it->first->compound->precursorMz > 0){
            key = static_cast<int>(mzUtils::mzToIntKey(static_cast<double>(it->first->compound->precursorMz)));
        } else {
            key = it->first->ms1IntensityCoord;
        }

        //theoretical precursor m/z not provided in compound, and no intensity detected in consensus ms1 spectrum
        if (key == -1) continue;

        if (partitionMap.find(key) == partitionMap.end()) {
            partitionMap.insert(make_pair(key, vector<shared_ptr<DirectInfusionMatchData>>()));
        }
        partitionMap[key].push_back(it->first);
    }

    //Issue 311: any matches without fragments that have the same theoretical m/z
    //as something else with matches
    //will always return a partition fraction of 0 (all intensity is taken by other compounds).
    //
    //if there is only one theoretical m/z, partition fraction is 1.
    for (auto compound : compoundsNoFragMatches) {

        //Issue 314: Prefer theoretical m/z as key over observed intensity from consensus ms1
        //This allows fractions to be computed when ms1 intensity can be identified in ms1 scans,
        //but not in the consensus ms1 spectrum.
        int key = -1;
        if (compound->compound->precursorMz > 0){
            key = static_cast<int>(mzUtils::mzToIntKey(static_cast<double>(compound->compound->precursorMz)));
        } else {
            key = compound->ms1IntensityCoord;
        }

        //theoretical precursor m/z not provided in compound, and no intensity detected in consensus ms1 spectrum
        if (key == -1) continue;

        if (partitionMap.find(key) == partitionMap.end()) {
            partitionMap.insert(make_pair(key, vector<shared_ptr<DirectInfusionMatchData>>()));
        }
        partitionMap[key].push_back(compound);
    }

    return partitionMap;
}

void DirectInfusionMatchAssessment::computeMs2MatchAssessment(
                               Fragment *observedSpectrum,
                               const Compound *compound,
                               const shared_ptr<DirectInfusionSearchParameters> params,
                               const bool debug){

    Fragment t;
    t.precursorMz = compound->precursorMz;
    t.mzs = compound->fragment_mzs;
    t.intensity_array = compound->fragment_intensity;
    t.fragment_labels = compound->fragment_labels;

    float maxDeltaMz = (params->ms2PpmTolr * static_cast<float>(t.precursorMz))/ 1000000;

    if (observedSpectrum) {
        fragmentationMatchScore.ranks = Fragment::findFragPairsGreedyMz(&t, observedSpectrum, maxDeltaMz);
    } else {
        //Issue 303: downstream analysis expects the ranks vector to exist and be the same size as the compound fragment vectors
        fragmentationMatchScore.ranks = vector<int>(compound->fragment_mzs.size(),-1);
    }

    bool isHasLabels = compound->fragment_labels.size() == fragmentationMatchScore.ranks.size();

    for (unsigned long i=0; i < fragmentationMatchScore.ranks.size(); i++) {

        int y = fragmentationMatchScore.ranks[i];

        if (y != -1 && observedSpectrum->intensity_array[y] >= params->ms2MinIntensity) {

            float fragmentObservedIntensity = observedSpectrum->intensity_array[y];

            if (fragmentObservedIntensity > fragmentMaxObservedIntensity) {
                fragmentMaxObservedIntensity = fragmentObservedIntensity;
            }

            fragmentationMatchScore.numMatches++;

            if (!isHasLabels) continue;

            vector<string> fragmentLabelTags = DirectInfusionMatchAssessment::getFragmentLabelTags(compound->fragment_labels[i], params, debug);

            if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2DiagnosticFragmentLabelTag") != fragmentLabelTags.end()) {
                fragmentationMatchScore.numDiagnosticMatches++;
            }

            if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn1FragmentLabelTag") != fragmentLabelTags.end()) {
                fragmentationMatchScore.numSn1Matches++;
            }

            if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn2FragmentLabelTag") != fragmentLabelTags.end()) {
                fragmentationMatchScore.numSn2Matches++;
            }

            if (abs(compound->fragment_mzs[i] - compound->precursorMz) < maxDeltaMz) {
                fragmentationMatchScore.isHasPrecursorMatch = true;
            }
        }
    }

}

//Issue 313
vector<string> DirectInfusionMatchAssessment::getFragmentLabelTags(string fragmentLabel,
                                                                   const shared_ptr<DirectInfusionSearchParameters> params,
                                                                   const bool debug){

    if (debug) {
        cout << "fragmentLabel: " << fragmentLabel << endl;
    }

    vector<string> singleFrags{};

    unsigned long posPrevious = 0;
    unsigned long posCurrent = 0;

    if (fragmentLabel.at(fragmentLabel.size()-1) != '/') {
        fragmentLabel.append("/");
    }

    while ((posCurrent = fragmentLabel.find("/", posPrevious)) != string::npos) {

        string singleFragmentLabel = fragmentLabel.substr(posPrevious, posCurrent-posPrevious);
        posPrevious = posCurrent + 1;

        singleFrags.push_back(singleFragmentLabel);
    }

    if (debug) {
        cout << "singleFrags:" << endl;
        for (auto lbl : singleFrags) {
            cout << "\t" << lbl << endl;
        }
    }

    bool isMs2DiagnosticFragmentLabel = false;
    bool isMs2sn1FragmentLabel = false;
    bool isMs2sn2FragmentLabel = false;

    for (auto singleFrag : singleFrags) {

        size_t posDiagnosticFragmentLabelTag = singleFrag.find(params->ms2DiagnosticFragmentLabelTag);
        size_t posSn1FragmentLabelTag = singleFrag.find(params->ms2sn1FragmentLabelTag);
        size_t posSn2FragmentLabelTag = singleFrag.find(params->ms2sn2FragmentLabelTag);

        if (debug) {

            string diagTag = to_string(posDiagnosticFragmentLabelTag);
            if (posDiagnosticFragmentLabelTag == string::npos) {
                diagTag = "NONE";
            }

            string sn1Tag = to_string(posSn1FragmentLabelTag);
            if (posSn1FragmentLabelTag == string::npos) {
                sn1Tag = "NONE";
            }

            string sn2Tag = to_string(posSn2FragmentLabelTag);
            if (posSn2FragmentLabelTag == string::npos) {
                sn2Tag = "NONE";
            }

            cout << "SINGLE FRAG:\""<< singleFrag << "\"; TAG POSITIONS: diagnostic=" << diagTag
                 << ", sn1=" << sn1Tag
                 << ", sn2=" << sn2Tag
                 << endl;
        }
        if (posDiagnosticFragmentLabelTag == 0) {
            isMs2DiagnosticFragmentLabel = true;

            if (posSn1FragmentLabelTag == params->ms2DiagnosticFragmentLabelTag.size()){
                isMs2sn1FragmentLabel = true;

                if (posSn2FragmentLabelTag == (params->ms2DiagnosticFragmentLabelTag.size()+params->ms2sn1FragmentLabelTag.size())){
                    isMs2sn2FragmentLabel = true;
                    break;
                }
            }

            if (posSn2FragmentLabelTag == params->ms2DiagnosticFragmentLabelTag.size()) {
                isMs2sn2FragmentLabel = true;

                if (posSn1FragmentLabelTag == (params->ms2DiagnosticFragmentLabelTag.size()+params->ms2sn2FragmentLabelTag.size())){
                    isMs2sn1FragmentLabel = true;
                    break;
                }
            }
        }

        if (posSn1FragmentLabelTag == 0) {
            isMs2sn1FragmentLabel = true;

            if (posDiagnosticFragmentLabelTag == params->ms2sn1FragmentLabelTag.size()) {
                isMs2DiagnosticFragmentLabel = true;

                if (posSn2FragmentLabelTag == (params->ms2sn1FragmentLabelTag.size() + params->ms2DiagnosticFragmentLabelTag.size())) {
                    isMs2sn2FragmentLabel = true;
                    break;
                }
            }

            if (posSn2FragmentLabelTag == params->ms2sn1FragmentLabelTag.size()) {
                isMs2sn2FragmentLabel = true;

                if (posDiagnosticFragmentLabelTag == (params->ms2sn1FragmentLabelTag.size() + params->ms2sn2FragmentLabelTag.size())) {
                    isMs2DiagnosticFragmentLabel = true;
                    break;
                }
            }
        }

        if (posSn2FragmentLabelTag == 0) {
            isMs2sn2FragmentLabel = true;

            if (posDiagnosticFragmentLabelTag == params->ms2sn2FragmentLabelTag.size()) {
                isMs2DiagnosticFragmentLabel = true;

                if (posSn1FragmentLabelTag == (params->ms2sn2FragmentLabelTag.size() + params->ms2DiagnosticFragmentLabelTag.size())) {
                    isMs2sn1FragmentLabel = true;
                    break;
                }
            }

            if (posSn1FragmentLabelTag == params->ms2sn2FragmentLabelTag.size()) {
                isMs2sn1FragmentLabel = true;

                if (posDiagnosticFragmentLabelTag == (params->ms2sn2FragmentLabelTag.size() + params->ms2sn1FragmentLabelTag.size())) {
                    isMs2DiagnosticFragmentLabel = true;
                    break;
                }
            }
        }

    }

    vector<string> fragmentLabelTags{};

    if (isMs2DiagnosticFragmentLabel) {
        fragmentLabelTags.push_back("ms2DiagnosticFragmentLabelTag");
    }
    if (isMs2sn1FragmentLabel) {
        fragmentLabelTags.push_back("ms2sn1FragmentLabelTag");
    }
    if (isMs2sn2FragmentLabel) {
        fragmentLabelTags.push_back("ms2sn2FragmentLabelTag");
    }

    return fragmentLabelTags;
}

string DirectInfusionMatchAssessment::getFragmentLabelWithoutTags(string fragmentLabel,
                                                                  const shared_ptr<DirectInfusionSearchParameters> params,
                                                                  const bool debug){

    if (debug) {
        cout << "fragmentLabel: " << fragmentLabel << endl;
    }

    vector<string> singleFrags{};

    unsigned long posPrevious = 0;
    unsigned long posCurrent = 0;

    if (fragmentLabel.at(fragmentLabel.size()-1) != '/') {
        fragmentLabel.append("/");
    }

    while ((posCurrent = fragmentLabel.find("/", posPrevious)) != string::npos) {

        string singleFragmentLabel = fragmentLabel.substr(posPrevious, posCurrent-posPrevious);
        posPrevious = posCurrent + 1;

        singleFrags.push_back(singleFragmentLabel);
    }

    if (debug) {
        cout << "singleFrags:" << endl;
        for (auto lbl : singleFrags) {
            cout << "\t" << lbl << endl;
        }
    }

    string cleanedFragmentLabel;
    for (unsigned int i = 0; i < singleFrags.size(); i++) {

        if (i > 0) {
            cleanedFragmentLabel.append("/");
        }

        string singleFrag = singleFrags[i];

        size_t posDiagnosticFragmentLabelTag = singleFrag.find(params->ms2DiagnosticFragmentLabelTag);
        size_t posSn1FragmentLabelTag = singleFrag.find(params->ms2sn1FragmentLabelTag);
        size_t posSn2FragmentLabelTag = singleFrag.find(params->ms2sn2FragmentLabelTag);

        vector<pair<size_t, string>> positions{};
        if (posDiagnosticFragmentLabelTag != string::npos) {
            positions.push_back(make_pair(posDiagnosticFragmentLabelTag, params->ms2DiagnosticFragmentLabelTag));
        }
        if (posSn1FragmentLabelTag != string::npos) {
            positions.push_back(make_pair(posSn1FragmentLabelTag, params->ms2sn1FragmentLabelTag));
        }
        if (posSn2FragmentLabelTag != string::npos) {
            positions.push_back(make_pair(posSn2FragmentLabelTag, params->ms2sn2FragmentLabelTag));
        }

        sort(positions.begin(), positions.end());

        if (debug) {
            cout << "positions: ";
            for (auto position : positions) {
                cout << position.first << " (" << position.second << ") ";
            }
            cout << endl;
        }

        unsigned long substringPosition = 0;
        for (auto position : positions) {
            if (position.first != substringPosition) {
                break;
            } else {
                substringPosition = substringPosition + position.second.size();
            }
        }

        string singleFragCleaned = singleFrag.substr(substringPosition);

        cleanedFragmentLabel.append(singleFragCleaned);
    }

    return cleanedFragmentLabel;
}

float DirectInfusionUtils::getObservedMs1ScanIntensity(
        const vector<Scan*>& validMs1Scans,
        float queryMz,
        shared_ptr<DirectInfusionSearchParameters> params,
        bool debug){

    float observedMs1ScanIntensity = 0.0f;
    vector<float> ms1ScanIntensities{};

    for (auto scan : validMs1Scans) {
        float queryMzIntensityCandidate = scan->findClosestMzIntensity(queryMz, params->ms1PpmTolr);
        if (queryMzIntensityCandidate >= params->ms1MinIntensity && queryMzIntensityCandidate > 0.0f) {
            ms1ScanIntensities.push_back(queryMzIntensityCandidate);
        }
    }

    if (!ms1ScanIntensities.empty()) {
        if (params->consensusIntensityAgglomerationType == Fragment::Mean) {
            observedMs1ScanIntensity = accumulate(ms1ScanIntensities.begin(), ms1ScanIntensities.end(), 0.0f) / ms1ScanIntensities.size();
        } else if (params->consensusIntensityAgglomerationType == Fragment::Median) {
            observedMs1ScanIntensity = median(ms1ScanIntensities);
        }
    }

    return observedMs1ScanIntensity;
}
