#include "secprocessor.h"
#include <numeric>

string SECSearchParameters::encodeParams() {

    string encodedParams;

    //program level
    encodedParams = encodedParams + "searchVersion" + "=" + searchVersion + ";";

    //SEC trace filling
    encodedParams = encodedParams + "traceMissingIntensityFill" + "=" + to_string(traceMissingIntensityFill) + ";";
    encodedParams = encodedParams + "traceMinFractionNumber" + "=" + to_string(traceMinFractionNumber) + ";";
    encodedParams = encodedParams + "traceMaxFractionNumber" + "=" + to_string(traceMaxFractionNumber) + ";";
    encodedParams = encodedParams + "traceNormalizeToSumIntensity" + "=" + to_string(traceNormalizeToSumIntensity) + ";";

    //Peak Picking

    encodedParams = encodedParams + "traceSmoothingType" + "=";
    if (traceSmoothingType == EIC::SmootherType::GAUSSIAN) {
        encodedParams = encodedParams + "GAUSSIAN";
    } else if (traceSmoothingType == EIC::SmootherType::AVG) {
        encodedParams = encodedParams + "AVG";
    } else if (traceSmoothingType == EIC::SmootherType::SAVGOL) {
        encodedParams = encodedParams + "SAVGOL";
    } else {
        encodedParams = encodedParams + "OTHER";
    }
    encodedParams = encodedParams + ";";

    encodedParams = encodedParams + "traceWindowSize" + "=" + to_string(traceWindowSize) + ";";
    encodedParams = encodedParams + "traceMinPeakIntensity" + "=" + to_string(traceMinPeakIntensity) + ";";
    encodedParams = encodedParams + "traceMinSmoothedIntensity" + "=" + to_string(traceMinSmoothedIntensity) + ";";
    encodedParams = encodedParams + "traceMinFracTopPeakIntensity" + "=" + to_string(traceMinFracTopPeakIntensity) + ";";
    encodedParams = encodedParams + "traceMinFracTopSmoothedIntensity" + "=" + to_string(traceMinFracTopSmoothedIntensity) + ";";
    encodedParams = encodedParams + "traceMinPeakSN" + "=" + to_string(traceMinPeakSN) + ";";
    encodedParams = encodedParams + "traceMinPeakWidth" + "=" + to_string(traceMinPeakWidth) + ";";
    encodedParams = encodedParams + "traceBaselineDropTopX" + "=" + to_string(traceBaselineDropTopX) + ";";
    encodedParams = encodedParams + "tracePeakBoundsMaxIntensityFraction" + "=" + to_string(tracePeakBoundsMaxIntensityFraction) + ";";
    encodedParams = encodedParams + "traceRtBoundsSlopeThreshold" + "=" + to_string(traceRtBoundsSlopeThreshold) + ";";

    // Fragment
    encodedParams = encodedParams + "fragmentIsSmoothedIntensity" + "=" + to_string(fragmentIsSmoothedIntensity) +";";

    // Trace Similarity Scoring
    encodedParams = encodedParams + "similarityMinNumPeaks" + "=" + to_string(similarityMinNumPeaks) + ";";
    encodedParams = encodedParams + "similarityFractionDiffTol" + "=" + to_string(similarityFractionDiffTol) + ";";

    // Peak Similarity Scoring
    encodedParams = encodedParams + "peakSimMaxCenterDiff" + "=" + to_string(peakSimMaxCenterDiff) + ";";

    return encodedParams;
}

shared_ptr<SECSearchParameters> SECSearchParameters::decode(string encodedParams) {

    shared_ptr<SECSearchParameters> secSearchParameters = shared_ptr<SECSearchParameters>(new SECSearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams);

    //program level
    if (decodedMap.find("searchVersion") != decodedMap.end()) {
        secSearchParameters->searchVersion = decodedMap["searchVersion"];
    }

    //SEC trace filling
    if (decodedMap.find("traceMissingIntensityFill") != decodedMap.end()) {
        secSearchParameters->traceMissingIntensityFill = stof(decodedMap["traceMissingIntensityFill"]);
    }
    if (decodedMap.find("traceMinFractionNumber") != decodedMap.end()) {
        secSearchParameters->traceMinFractionNumber = stoi(decodedMap["traceMinFractionNumber"]);
    }
    if (decodedMap.find("traceMaxFractionNumber") != decodedMap.end()) {
        secSearchParameters->traceMaxFractionNumber = stoi(decodedMap["traceMaxFractionNumber"]);
    }
    if (decodedMap.find("traceNormalizeToSumIntensity") != decodedMap.end()) {
        secSearchParameters->traceNormalizeToSumIntensity = decodedMap["traceNormalizeToSumIntensity"] == "1";
    }

    //Peak Picking
    if (decodedMap.find("traceSmoothingType") != decodedMap.end()) {
        string traceSmoothingType = decodedMap["traceSmoothingType"];
        if (traceSmoothingType == "GAUSSIAN") {
            secSearchParameters->traceSmoothingType = EIC::SmootherType::GAUSSIAN;
        } else if (traceSmoothingType == "AVG") {
            secSearchParameters->traceSmoothingType = EIC::SmootherType::AVG;
        } else if (traceSmoothingType == "SAVGOL") {
            secSearchParameters->traceSmoothingType = EIC::SmootherType::SAVGOL;
        }
    }
    if (decodedMap.find("traceWindowSize") != decodedMap.end()) {
        secSearchParameters->traceWindowSize = stoi(decodedMap["traceWindowSize"]);
    }
    if (decodedMap.find("traceMinPeakIntensity") != decodedMap.end()) {
        secSearchParameters->traceMinPeakIntensity = stof(decodedMap["traceMinPeakIntensity"]);
    }
    if (decodedMap.find("traceMinSmoothedIntensity") != decodedMap.end()) {
        secSearchParameters->traceMinSmoothedIntensity = stof(decodedMap["traceMinSmoothedIntensity"]);
    }
    if (decodedMap.find("traceMinFracTopPeakIntensity") != decodedMap.end()) {
        secSearchParameters->traceMinFracTopPeakIntensity = stof(decodedMap["traceMinFracTopPeakIntensity"]);
    }
    if (decodedMap.find("traceMinFracTopSmoothedIntensity") != decodedMap.end()) {
        secSearchParameters->traceMinFracTopSmoothedIntensity = stof(decodedMap["traceMinFracTopSmoothedIntensity"]);
    }
    if (decodedMap.find("traceMinPeakSN") != decodedMap.end()) {
        secSearchParameters->traceMinPeakSN = stof(decodedMap["traceMinPeakSN"]);
    }
    if (decodedMap.find("traceMinPeakWidth") != decodedMap.end()) {
        secSearchParameters->traceMinPeakWidth = stoi(decodedMap["traceMinPeakWidth"]);
    }
    if (decodedMap.find("traceBaselineDropTopX") != decodedMap.end()) {
        secSearchParameters->traceBaselineDropTopX = stoi(decodedMap["traceBaselineDropTopX"]);
    }
    if (decodedMap.find("tracePeakBoundsMaxIntensityFraction") != decodedMap.end()) {
        secSearchParameters->tracePeakBoundsMaxIntensityFraction = stof(decodedMap["tracePeakBoundsMaxIntensityFraction"]);
    }
    if (decodedMap.find("traceRtBoundsSlopeThreshold") != decodedMap.end()) {
        secSearchParameters->traceRtBoundsSlopeThreshold = stof(decodedMap["traceRtBoundsSlopeThreshold"]);
    }

    // Fragment
    if (decodedMap.find("fragmentIsSmoothedIntensity") != decodedMap.end()) {
        secSearchParameters->fragmentIsSmoothedIntensity = decodedMap["fragmentIsSmoothedIntensity"] == "1";
    }

    // Similarity Scoring
    if (decodedMap.find("similarityMinNumPeaks") != decodedMap.end()) {
        secSearchParameters->similarityMinNumPeaks = stoi(decodedMap["similarityMinNumPeaks"]);
    }
    if (decodedMap.find("similarityFractionDiffTol") != decodedMap.end()) {
        secSearchParameters->similarityFractionDiffTol = stoi(decodedMap["similarityFractionDiffTol"]);
    }

    // Peak Similarity Scoring
    if (decodedMap.find("peakSimMaxCenterDiff") != decodedMap.end()) {
        secSearchParameters->peakSimMaxCenterDiff = stoi(decodedMap["peakSimMaxCenterDiff"]);
    }

    return secSearchParameters;
}

SECTrace::SECTrace(string id,
                   SECTraceType type,
                   vector<int> fractionNums,
                   vector<float> rawIntensities,
                   shared_ptr<SECSearchParameters> params,
                   bool debug) {

    if (fractionNums.size() != rawIntensities.size()) {
        cerr << "SECTrace() requires same length for fractionNums, rawIntensities. Exiting." << endl;
        abort();
    }

    this->id = id;
    this->params = params;
    this->type = type;

    int N = params->traceMaxFractionNumber - params->traceMinFractionNumber + 1;

    this->fractionNums = vector<int>(static_cast<unsigned long>(N));
    this->rawIntensities = vector<float>(static_cast<unsigned long>(N));

    map<int, float> rawDataMap{};
    for (unsigned int i = 0; i < fractionNums.size(); i++) {
        int fracNum = fractionNums[i];
        float rawIntensity = rawIntensities[i];
        rawDataMap.insert(make_pair(fracNum, rawIntensity));
    }

    unsigned int counter = 0;

    int fracNum = params->traceMinFractionNumber;

    vector<float> pseudoRt(static_cast<unsigned long>(N));

    while (fracNum <= params->traceMaxFractionNumber) {
        this->fractionNums[counter] = fracNum;
        pseudoRt[counter] = static_cast<float>(fracNum);

        float intensityVal = params->traceMissingIntensityFill;

        if (rawDataMap.find(fracNum) != rawDataMap.end()) {
            intensityVal = rawDataMap.at(fracNum);
        }

        this->rawIntensities[counter] = intensityVal;

        counter++;
        fracNum++;
    }

    if (params->traceNormalizeToSumIntensity) {
        float intensitySum = accumulate(this->rawIntensities.begin(), this->rawIntensities.end(), 0.0f);

        //guard against divide by 0
        if (intensitySum == 0.0f) {
            intensitySum = 1.0f;
        }

        for (unsigned int i = 0; i < this->rawIntensities.size(); i++) {
            this->rawIntensities[i] = this->rawIntensities[i]/intensitySum;
        }
    }

    EIC *eic = new EIC();

    eic->intensity = this->rawIntensities;
    eic->rt = pseudoRt; // x axis: fractions
    eic->mz = pseudoRt; // x axis: fractions
    eic->scannum = this->fractionNums;

    eic->setSmootherType(params->traceSmoothingType);
    eic->setBaselineSmoothingWindow(params->traceWindowSize);
    eic->setBaselineDropTopX(params->traceBaselineDropTopX);

    eic->getPeakPositionsC(
                params->traceWindowSize,
                debug,
                true,
                params->tracePeakBoundsMaxIntensityFraction,
                params->traceRtBoundsSlopeThreshold);

    this->smoothedIntensities = eic->spline;

    float maxRawIntensity = *max_element(this->rawIntensities.begin(), this->rawIntensities.end());
    float maxSmoothedIntensity = *max_element(this->smoothedIntensities.begin(), this->smoothedIntensities.end());

    float rawIntensityThreshold = max(params->traceMinFracTopPeakIntensity * maxRawIntensity, params->traceMinPeakIntensity);
    float smoothedIntensityThreshold = max(params->traceMinFracTopSmoothedIntensity * maxSmoothedIntensity, params->traceMinSmoothedIntensity);

    //float minPeakIntensity = maxRawIntensity * params->trace
    for (auto p : eic->peaks) {

        //Issue 598: peak width has a different meaning for LC data, involving noise estimates
        //here, it simply means the total # of fractions the peak spans
        //fractions are numbered as integers, so position difference can be used.
        p.width = p.maxpos - p.minpos;

        float peakRawIntensity = this->rawIntensities[p.pos];
        float peakSmoothedIntensity = this->smoothedIntensities[p.pos];

        if (peakRawIntensity >= rawIntensityThreshold && peakSmoothedIntensity >= smoothedIntensityThreshold && p.signalBaselineRatio >= params->traceMinPeakSN) {
            this->peaks.push_back(p);
        }
    }

    delete(eic);
}

vector<string> SECTrace::getPeakSummaryString(
        string empty,
        string leftPrefix,
        string maxPrefix,
        string rightPrefix){

    vector<string> summaryString(this->fractionNums.size(), empty);

    for (unsigned int i = 0; i < this->peaks.size(); i++){

        Peak p = this->peaks[i];

        string left = leftPrefix + to_string(i);
        string max = maxPrefix + to_string(i);
        string right = rightPrefix + to_string(i);

        //Issue 598 testing: use saved indices of smoothed intensity vector instead of from fraction values
        unsigned int left_coord = p.minpos;
        unsigned int max_coord = p.pos;
        unsigned int right_coord = p.maxpos;

        if (summaryString.at(left_coord) != empty) {
            summaryString[left_coord] = summaryString[left_coord] + ", " + left;
        } else {
            summaryString[left_coord] = left;
        }

        if (summaryString.at(max_coord) != empty) {
            summaryString[max_coord] = summaryString[max_coord] + ", " + max;
        } else {
            summaryString[max_coord] = max;
        }

        if (summaryString.at(right_coord) != empty) {
            summaryString[right_coord] = summaryString[right_coord] + ", " + right;
        } else {
            summaryString[right_coord] = right;
        }
    }

    return summaryString;
}

Fragment* SECTrace::getFragment(shared_ptr<SECSearchParameters> params, bool debug) {

    if (fragment) return fragment;

    if (debug) {
        cout << "[SECTrace::getFragment()]: Computing fragment." << endl;
    }

    fragment = new Fragment();

    for (auto p : peaks) {
        fragment->mzs.push_back(fractionNums[p.pos]);

        float intensityVal = params->fragmentIsSmoothedIntensity ? smoothedIntensities[p.pos] : rawIntensities[p.pos];
        fragment->intensity_array.push_back(intensityVal);
    }

    return fragment;
}

float SECTraceSimilarity::getSimilarity(bool debug){
    if (debug) cout << "SECTraceSimilarity::getSimilarity()" << endl;
    return similarity;
}

SECTraceSimilarityCosine::SECTraceSimilarityCosine(SECTrace* first, SECTrace* second, shared_ptr<SECSearchParameters> params) {
    this->first = first;
    this->second = second;
    this->params = params;
    this->compareId = first->id + "_" + second->id;
}

float SECTraceSimilarityCosine::getSimilarity(bool debug) {
    if (!first || !second || similarity > -1.0f) return similarity;

    unsigned int minPeaks = static_cast<unsigned int>(params->similarityMinNumPeaks);
    if (first->peaks.size() < minPeaks || second->peaks.size() < minPeaks) return similarity;

    Fragment *f1 = first->getFragment(params, debug);
    Fragment *f2 = second->getFragment(params, debug);

    float productPpmTolr = params->similarityFractionDiffTol + 0.001f; // avoid rounding errors

    this->ranks = Fragment::compareRanks(f1, f2, productPpmTolr);

    //all scores
    for(int rank: ranks){
        if(rank != -1) numPeakMatches++;
    }
    fractionPeaksMatched = 2.0f * numPeakMatches / (f1->mzs.size() + f2->mzs.size());
    cosineScore = static_cast<float>(Fragment::normCosineScore(f1, f2, ranks));
    matchedPeakCosineScore = static_cast<float>(Fragment::matchedPeakCosineScore(f1, f2, ranks));
    pearsonCorrelationRaw = mzUtils::correlation(first->rawIntensities, second->rawIntensities);
    pearsonCorrelationSmoothed = mzUtils::correlation(first->smoothedIntensities, second->smoothedIntensities);

    //main score
    similarity = matchedPeakCosineScore;

    return similarity;
}

vector<SECTraceSimilarityCosine> SECTraceCosineSimilarityScorer::scoreTraces(
        vector<SECTrace*> traces,
        shared_ptr<SECSearchParameters> params,
        bool debug){

    vector<SECTraceSimilarityCosine> similarityScores{};

    sort(traces.begin(), traces.end(), [](SECTrace* lhs, SECTrace* rhs){
       return lhs->id < rhs->id;
    });

    for (unsigned int i = 0; i < traces.size(); i++) {

        auto ithTrace = traces[i];

        if (static_cast<int>(ithTrace->peaks.size()) < params->similarityMinNumPeaks) continue;

        for (unsigned int j = i+1; j < traces.size(); j++) {

            auto jthTrace = traces[j];

            if (static_cast<int>(jthTrace->peaks.size()) < params->similarityMinNumPeaks) continue;

            auto similarityScore = SECTraceSimilarityCosine(ithTrace, jthTrace, params);
            similarityScore.getSimilarity(debug);

            similarityScores.push_back(similarityScore);

        }
    }

    sort(similarityScores.begin(), similarityScores.end(), [](SECTraceSimilarityCosine& lhs, SECTraceSimilarityCosine& rhs){
       if (lhs.matchedPeakCosineScore != rhs.matchedPeakCosineScore) {
           return lhs.matchedPeakCosineScore < rhs.matchedPeakCosineScore;
       }  else if (lhs.cosineScore != rhs.cosineScore) {
           return lhs.cosineScore < rhs.cosineScore;
       } else {
           return lhs.compareId < rhs.compareId;
       }
    });

    return similarityScores;
}

bool SECTracePeak::isValid() {
    if (!trace || peakNum < 0) return false;
    if (static_cast<unsigned int>(peakNum) >= trace->peaks.size()) return false;

    return true;
}

int SECTracePeak::getPeakFractionNum(){
    if (!isValid()) return -1;
    return static_cast<int>(trace->peaks[static_cast<unsigned int>(peakNum)].rt);
}

int SECTracePeak::getMinFractionNum(){
    if (!isValid()) return -1;
    return static_cast<int>(trace->peaks[static_cast<unsigned int>(peakNum)].rtmin);
}

int SECTracePeak::getMaxFractionNum(){
    if (!isValid()) return -1;
    return static_cast<int>(trace->peaks[static_cast<unsigned int>(peakNum)].rtmax);
}

int SECTracePeak::getPeakIndex() {
    if (!isValid()) return -1;
    Peak p = trace->peaks[static_cast<unsigned int>(peakNum)];
    return static_cast<int>(p.pos - p.minpos);
}

vector<float> SECTracePeak::getSmoothedIntensities(pair<int, int> comparableRange){
    if (!isValid()) return vector<float>{};

    Peak p = trace->peaks[static_cast<unsigned int>(peakNum)];
    vector<float> peakSmoothedIntensities(p.width+1);

    unsigned int index = 0;
    for (unsigned int i = p.minpos; i <= p.maxpos; i++) {
        peakSmoothedIntensities[index] = trace->smoothedIntensities[i];
        index++;
    }

    if (comparableRange.first >= 0 && comparableRange.second >= 0) {
        peakSmoothedIntensities = getComparableRangeSubset(peakSmoothedIntensities, comparableRange);
    }

    return peakSmoothedIntensities;
}

vector<float> SECTracePeak::getRawIntensities(pair<int, int> comparableRange){
    if (!isValid()) return vector<float>{};

    Peak p = trace->peaks[static_cast<unsigned int>(peakNum)];
    vector<float> peakRawIntensities(p.width+1);

    unsigned int index = 0;
    for (unsigned int i = p.minpos; i <= p.maxpos; i++) {
        peakRawIntensities[index] = trace->rawIntensities[i];
        index++;
    }

    if (comparableRange.first >= 0 && comparableRange.second >= 0) {
        peakRawIntensities = getComparableRangeSubset(peakRawIntensities, comparableRange);
    }

    return peakRawIntensities;
}

vector<int> SECTracePeak::getFractionNums(pair<int, int> comparableRange) {
    if (!isValid()) return vector<int>{};
    Peak p = trace->peaks[static_cast<unsigned int>(peakNum)];

    int left_coord = static_cast<int>(p.rtmin);
    int right_coord = static_cast<int>(p.rtmax);

    vector<int> peakFractions(p.width+1);

    unsigned int index = 0;
    for (int i = left_coord; i <= right_coord; i++) {
        peakFractions[index] = i;
        index++;
    }

    if (comparableRange.first >= 0 && comparableRange.second >= 0) {
        peakFractions = getComparableRangeSubset(peakFractions, comparableRange);
    }

    return peakFractions;
}

string SECTracePeak::getPeakId() {
    if (!isValid()) return "no_peak";
    return trace->id + "_" + to_string(peakNum);
}

string SECTracePeakComparison::getPeakComparisonId() {
    return first.getPeakId() + "_" + second.getPeakId();
}

SECTracePeak::SECTracePeak(SECTrace *trace, int peakNum){
    this->trace = trace;
    this->peakNum = peakNum;
}

vector<float> SECTracePeak::getComparableRangeSubset(const vector<float>&x, pair<int, int> comparableRange){
    unsigned int N = comparableRange.first + comparableRange.second + 1;

    vector<float> subset(N);

    unsigned int index = 0;

    unsigned int start = static_cast<unsigned int>(getPeakIndex() - comparableRange.first);
    unsigned int stop = static_cast<unsigned int>(getPeakIndex() + comparableRange.second);

    for (unsigned int i = start; i <= stop; i++) {
        subset[index] = x[i];
        index++;
    }

    return subset;
}

vector<int> SECTracePeak::getComparableRangeSubset(const vector<int>& x, pair<int, int> comparableRange) {
    unsigned int N = comparableRange.first + comparableRange.second + 1;

    vector<int> subset(N);

    unsigned int index = 0;

    unsigned int start = static_cast<unsigned int>(getPeakIndex() - comparableRange.first);
    unsigned int stop = static_cast<unsigned int>(getPeakIndex() + comparableRange.second);

    for (unsigned int i = start; i <= stop; i++) {
        subset[index] = x[i];
        index++;
    }

    return subset;
}

void SECTracePeakComparison::computeComparableRange(){
    int firstPeakIndex = first.getPeakIndex();
    int secondPeakIndex = second.getPeakIndex();

    int firstPointsLeftOfMax = firstPeakIndex;
    int firstPointsRightOfMax = first.getRawIntensities().size()-firstPeakIndex-1; // exclude max

    int secondPointsLeftOfMax = secondPeakIndex;
    int secondPointsRightOfMax = second.getRawIntensities().size()-secondPeakIndex-1;

    int finalPointsLeftOfMax = min(firstPointsLeftOfMax, secondPointsLeftOfMax);
    int finalPointsRightOfMax = min(firstPointsRightOfMax, secondPointsRightOfMax);

    comparableRange = make_pair(finalPointsLeftOfMax, finalPointsRightOfMax);
}

void SECTracePeakComparison::computeSecFractionJaccard(){

    int intersectMinFrac = max(first.getMinFractionNum(), second.getMinFractionNum());
    int intersectMaxFrac = min(first.getMaxFractionNum(), second.getMaxFractionNum());

    int minFrac = min(first.getMinFractionNum(), second.getMinFractionNum());
    int maxFrac = max(first.getMaxFractionNum(), second.getMaxFractionNum());

    int numIntersect = intersectMaxFrac - intersectMinFrac + 1;
    int numUnion = maxFrac - minFrac + 1;

    secFractionJaccard = static_cast<float>(numIntersect)/static_cast<float>(numUnion);
}

SECTracePeakComparison::SECTracePeakComparison(SECTrace *firstTrace, int firstPeakNum, SECTrace *secondTrace, int secondPeakNum){
    this->first = SECTracePeak(firstTrace, firstPeakNum);
    this->second = SECTracePeak(secondTrace, secondPeakNum);

    computeComparableRange();

    pearsonCorrelationSmoothed = mzUtils::correlation(first.getSmoothedIntensities(comparableRange), second.getSmoothedIntensities(comparableRange));
    pearsonCorrelationRaw = mzUtils::correlation(first.getRawIntensities(comparableRange), second.getSmoothedIntensities(comparableRange));

    secFractionOverlap = mzUtils::checkOverlap(
                first.getMinFractionNum(), first.getMaxFractionNum(),
                second.getMinFractionNum(), second.getMaxFractionNum());

    computeSecFractionJaccard();

    peakCenterDistance = abs(first.getPeakFractionNum() - second.getPeakFractionNum());

}

void SECTracePeakComparison::printSummary() {
    cout << getPeakComparisonId() << ":\n"
         << "\tpearsonCorrelationSmoothed: " << pearsonCorrelationSmoothed << "\n"
         << "\tpearsonCorrelationRaw: " << pearsonCorrelationRaw << "\n"
         << "\tsecFractionOverlap: " << secFractionOverlap << "\n"
         << "\tsecFractionJaccard: " << secFractionJaccard << "\n"
         << "\tpeakCenterDistance: " << peakCenterDistance
         << endl;
}

vector<SECTracePeakComparison> SECTracePeakScorer::scorePeaks(
        vector<SECTrace*> traces,
        shared_ptr<SECSearchParameters> params,
        bool debug){

    vector<SECTracePeakComparison> peakComparisons{};

    sort(traces.begin(), traces.end(), [](SECTrace* lhs, SECTrace* rhs){
       return lhs->id < rhs->id;
    });

    for (unsigned int i = 0; i < traces.size(); i++) {

        auto ithTrace = traces[i];

        if (ithTrace->peaks.empty()) continue;

        for (unsigned int j = i+1; j < traces.size(); j++) {

            auto jthTrace = traces[j];

            if (jthTrace->peaks.empty()) continue;

            for (unsigned int k = 0; i < ithTrace->peaks.size(); k++) {
                Peak peakI = ithTrace->peaks.at(k);
                for (unsigned int l = 0; l < jthTrace->peaks.size(); l++) {
                    Peak peakJ = jthTrace->peaks.at(l);
                    int fracDiff = static_cast<int>(abs(peakI.rt - peakJ.rt));
                    if (fracDiff < params->peakSimMaxCenterDiff) {
                        SECTracePeakComparison comparison = SECTracePeakComparison(
                                    ithTrace, static_cast<int>(k),
                                    jthTrace, static_cast<int>(l));
                        peakComparisons.push_back(comparison);
                    }
                }
            }
        }
    }

    sort(peakComparisons.begin(), peakComparisons.end(), [](SECTracePeakComparison& lhs, SECTracePeakComparison& rhs){
        if (lhs.pearsonCorrelationSmoothed == rhs.pearsonCorrelationSmoothed) {
            if (lhs.pearsonCorrelationRaw == rhs.pearsonCorrelationRaw) {
                return lhs.getPeakComparisonId() <  rhs.getPeakComparisonId();
            } else {
                return lhs.pearsonCorrelationRaw < rhs.pearsonCorrelationRaw;
            }
        } else {
            return lhs.pearsonCorrelationSmoothed < rhs.pearsonCorrelationSmoothed;
        }
    });

    return peakComparisons;
}
