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
    encodedParams = encodedParams + "traceIsPickEdgePeaks" + "=" + to_string(traceIsPickEdgePeaks) + ";";

    // Grouping
    encodedParams = encodedParams + "groupMaxFracDiff" + "=" + to_string(groupMaxFracDiff) + ";";
    encodedParams = encodedParams + "groupMergeOverlap" + "=" + to_string(groupMergeOverlap) + ";";
    encodedParams = encodedParams + "groupIsMergeOverlappingPeakGroups" + "=" + to_string(groupIsMergeOverlappingPeakGroups) + ";";
    encodedParams = encodedParams + "groupMinNumPeaks" + "=" + to_string(groupMinNumPeaks) + ";";

    // Fragment
    encodedParams = encodedParams + "fragmentIsSmoothedIntensity" + "=" + to_string(fragmentIsSmoothedIntensity) +";";

    // Trace Similarity Scoring
    encodedParams = encodedParams + "similarityMinNumPeaks" + "=" + to_string(similarityMinNumPeaks) + ";";
    encodedParams = encodedParams + "similarityFractionDiffTol" + "=" + to_string(similarityFractionDiffTol) + ";";

    // Peak Similarity Scoring
    encodedParams = encodedParams + "peakSimMaxCenterDiff" + "=" + to_string(peakSimMaxCenterDiff) + ";";
    encodedParams = encodedParams + "peakSimMinSecFractionOverlap" + "=" + to_string(peakSimMinSecFractionOverlap) + ";";
    encodedParams = encodedParams + "peakSimMinSecFractionJaccard" + "=" + to_string(peakSimMinSecFractionJaccard) + ";";
    encodedParams = encodedParams + "peakSimMinSmoothedCorrelation" + "=" + to_string(peakSimMinSmoothedCorrelation) + ";";
    encodedParams = encodedParams + "peakSimMinRawCorrelation" + "=" + to_string(peakSimMinRawCorrelation) + ";";
    encodedParams = encodedParams + "peakSimMinFractionNumber" + "=" + to_string(peakSimMinFractionNumber) + ";";
    encodedParams = encodedParams + "peakSimMaxFractionNumber" + "=" + to_string(peakSimMaxFractionNumber) + ";";

    return encodedParams;
}

shared_ptr<PeakPickingAndGroupingParameters> SECSearchParameters::toPeakPickingAndGroupingParams() {
    shared_ptr<PeakPickingAndGroupingParameters> peakPickingAndGroupingParams = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());

    // peak picking

    peakPickingAndGroupingParams->peakSmoothingWindow = traceWindowSize;
    peakPickingAndGroupingParams->peakRtBoundsMaxIntensityFraction = tracePeakBoundsMaxIntensityFraction;
    peakPickingAndGroupingParams->peakRtBoundsSlopeThreshold = traceRtBoundsSlopeThreshold;
    peakPickingAndGroupingParams->peakBaselineSmoothingWindow = traceWindowSize;
    peakPickingAndGroupingParams->peakBaselineDropTopX = traceBaselineDropTopX;
    peakPickingAndGroupingParams->peakRtBoundsMaxIntensityFraction = tracePeakBoundsMaxIntensityFraction;
    peakPickingAndGroupingParams->peakIsPickEdgePeaks = traceIsPickEdgePeaks;

    peakPickingAndGroupingParams->peakIsComputeBounds = true;
    peakPickingAndGroupingParams->peakIsReassignPosToUnsmoothedMax = false;

    // peak grouping
    peakPickingAndGroupingParams->mergedSmoothingWindow = traceWindowSize;
    peakPickingAndGroupingParams->mergedPeakRtBoundsMaxIntensityFraction = tracePeakBoundsMaxIntensityFraction;
    peakPickingAndGroupingParams->mergedPeakRtBoundsSlopeThreshold = traceRtBoundsSlopeThreshold;
    peakPickingAndGroupingParams->mergedBaselineSmoothingWindow = traceWindowSize;
    peakPickingAndGroupingParams->mergedBaselineDropTopX = traceBaselineDropTopX;

    peakPickingAndGroupingParams->groupMaxRtDiff = groupMaxFracDiff;
    peakPickingAndGroupingParams->groupMergeOverlap = groupMergeOverlap;
    peakPickingAndGroupingParams->groupIsMergeOverlappingPeakGroups = groupIsMergeOverlappingPeakGroups;

    return peakPickingAndGroupingParams;
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
    if (decodedMap.find("traceIsPickEdgePeaks") != decodedMap.end()) {
        secSearchParameters->traceIsPickEdgePeaks = decodedMap["traceIsPickEdgePeaks"] == "1";
    }

    //Grouping
    if (decodedMap.find("groupMaxFracDiff") != decodedMap.end()) {
        secSearchParameters->groupMaxFracDiff = stof(decodedMap["groupMaxFracDiff"]);
    }
    if (decodedMap.find("groupMergeOverlap") != decodedMap.end()) {
        secSearchParameters->groupMergeOverlap = stof(decodedMap["groupMergeOverlap"]);
    }
    if (decodedMap.find("groupIsMergeOverlappingPeakGroups") != decodedMap.end()) {
        secSearchParameters->groupIsMergeOverlappingPeakGroups = decodedMap["groupIsMergeOverlappingPeakGroups"] == "1";
    }
    if (decodedMap.find("groupMinNumPeaks") != decodedMap.end()) {
        secSearchParameters->groupMinNumPeaks = stoi(decodedMap["groupMinNumPeaks"]);
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
    if (decodedMap.find("peakSimMinSecFractionOverlap") != decodedMap.end()) {
        secSearchParameters->peakSimMinSecFractionOverlap = stof(decodedMap["peakSimMinSecFractionOverlap"]);
    }
    if (decodedMap.find("peakSimMinSecFractionJaccard") != decodedMap.end()) {
        secSearchParameters->peakSimMinSecFractionJaccard = stof(decodedMap["peakSimMinSecFractionJaccard"]);
    }
    if (decodedMap.find("peakSimMinSmoothedCorrelation") != decodedMap.end()) {
        secSearchParameters->peakSimMinSmoothedCorrelation = stof(decodedMap["peakSimMinSmoothedCorrelation"]);
    }
    if (decodedMap.find("peakSimMinRawCorrelation") != decodedMap.end()) {
        secSearchParameters->peakSimMinRawCorrelation = stof(decodedMap["peakSimMinRawCorrelation"]);
    }
    if (decodedMap.find("peakSimMinFractionNumber") != decodedMap.end()) {
        secSearchParameters->peakSimMinFractionNumber = stoi(decodedMap["peakSimMinFractionNumber"]);
    }
    if (decodedMap.find("peakSimMaxFractionNumber") != decodedMap.end()) {
        secSearchParameters->peakSimMaxFractionNumber = stoi(decodedMap["peakSimMaxFractionNumber"]);
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

    this->computeTraceData(fractionNums, rawIntensities, params, debug);
    this->pickPeaks(debug);
}

SECTraceDiff::SECTraceDiff(SECTrace *compare, SECTrace *reference, bool debug) {

    if (!compare || !reference) {
        cerr << "SECTraceDiff must intake two non-null SECTrace*. Exiting." << endl;
        abort();
    }

    if (compare->fractionNums.empty() || reference->fractionNums.empty()) {
        cerr << "input SECTrace objects contain no fractionNums. Exiting." << endl;
        abort();
    }

    if (compare->fractionNums.size() != reference->fractionNums.size()) {
        cerr << "compare and reference must have a fractionNums vector of the same size. Exiting." << endl;
        abort();
    }

    for (unsigned int i = 0; i < reference->fractionNums.size(); i++) {
        if (reference->fractionNums[i] != compare->fractionNums[i]) {
            cerr << "reference->fractionNums[" << i << "]="
                 << reference->fractionNums[i]
                 << " != "
                 << "compare->fractionNums[" << i << "]="
                 << compare->fractionNums[i]
                 << "."
                 << endl;
            cerr << "compare and reference must have identical fractionNums vectors. Exiting." << endl;
            abort();
        }
    }

    if (compare->smoothedIntensities.size() != compare->rawIntensities.size()) {
        cerr << "compare has not undergone peak-picking. Exiting." << endl;
        abort();
    }
    if (reference->smoothedIntensities.size() != reference->rawIntensities.size()) {
        cerr << "reference has not undergone peak-picking. Exiting." << endl;
        abort();
    }

    string id = reference->id;
    if (reference->id != compare->id) {
        id = "ref=" + reference->id + "_vs_comp=" + compare->id;
    }

    this->id = id;
    this->type = reference->type;
    this->params = reference->params;
    this->fractionNums = reference->fractionNums;

    unsigned long N = reference->fractionNums.size();

    this->diffRawIntensities = vector<float>(N);
    this->diffSmoothedIntensities = vector<float>(N);
    this->rawIntensities = vector<float>(N); //absolute difference between traces

    for (unsigned int i = 0; i < N; i++) {

        float rawDiff = compare->rawIntensities[i] - reference->rawIntensities[i];

        this->diffRawIntensities[i] = rawDiff;
        this->rawIntensities[i] = abs(rawDiff);

        this->diffSmoothedIntensities[i] = compare->smoothedIntensities[i] - reference->smoothedIntensities[i];
    }

    this->pickPeaks(debug);

    this->similarityScore = new SECTraceSimilarityCosine(compare, reference, params);
    this->similarityScore->getSimilarity(debug); //side effects: fill out fields
}

void SECTrace::computeTraceData(
        vector<int> fractionNums,
        vector<float> rawIntensities,
        shared_ptr<SECSearchParameters> params,
        bool debug){

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

    while (fracNum <= params->traceMaxFractionNumber) {
        this->fractionNums[counter] = fracNum;

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
}

void SECTrace::pickPeaks(bool debug) {
    eic = new EIC();

    vector<float> pseudoRt(fractionNums.size());
    for (unsigned int i = 0; i < fractionNums.size(); i++) {
        pseudoRt[i] = static_cast<float>(fractionNums[i]);
    }

    eic->intensity = this->rawIntensities;
    eic->rt = pseudoRt; // x axis: fractions
    eic->mz = pseudoRt; // x axis: fractions
    eic->scannum = this->fractionNums;

    eic->rtmin = pseudoRt[0];
    eic->rtmax = pseudoRt[pseudoRt.size()-1];

    eic->setSmootherType(params->traceSmoothingType);
    eic->setBaselineSmoothingWindow(params->traceWindowSize);
    eic->setBaselineDropTopX(params->traceBaselineDropTopX);

    //Issue 740: Update from EIC::getPeakPositionsD() to EIC::getPeakPositions()
    eic->getPeakPositionsD(params->toPeakPickingAndGroupingParams(), debug);

    this->smoothedIntensities = eic->spline;

    float maxRawIntensity = *max_element(this->rawIntensities.begin(), this->rawIntensities.end());
    float maxSmoothedIntensity = *max_element(this->smoothedIntensities.begin(), this->smoothedIntensities.end());

    float rawIntensityThreshold = max(params->traceMinFracTopPeakIntensity * maxRawIntensity, params->traceMinPeakIntensity);
    float smoothedIntensityThreshold = max(params->traceMinFracTopSmoothedIntensity * maxSmoothedIntensity, params->traceMinSmoothedIntensity);

    if (debug) {
        cout << "SECTrace::pickPeaks():\n"
             << "\tmaxRawIntensity: " << maxRawIntensity << "\n"
             << "\tmaxSmoothedIntensity: " << maxSmoothedIntensity << "\n"
             << "\trawIntensityThreshold: " << rawIntensityThreshold << "\n"
             << "\tsmoothedIntensityThreshold: " << smoothedIntensityThreshold << "\n"
             << endl;
    }

    //float minPeakIntensity = maxRawIntensity * params->trace

    //EIC peaks are subject to additional constraints
    vector<Peak> validPeaks{};

    for (auto & p : eic->peaks) {

        //Issue 598: peak width has a different meaning for LC data, involving noise estimates
        //here, it simply means the total # of fractions the peak spans
        //fractions are numbered as integers, so position difference can be used.
        p.width = p.maxpos - p.minpos;

        float peakRawIntensity = this->rawIntensities[p.pos];
        float peakSmoothedIntensity = this->smoothedIntensities[p.pos];

        if (debug) {
            cout << "Peak @ pos=" << p.pos
                 << ": width=" << p.width
                 << ", peakRawIntensity=" << peakRawIntensity
                 << ", peakSmoothedIntensity=" << peakSmoothedIntensity;
        }

        if (peakRawIntensity >= rawIntensityThreshold
                && peakSmoothedIntensity >= smoothedIntensityThreshold
                && p.signalBaselineRatio >= params->traceMinPeakSN
                && p.width >= static_cast<unsigned int>(params->traceMinPeakWidth)) {
            validPeaks.push_back(p);
            if (debug) cout << " --> passing";
        }

        if (debug) cout << endl;
    }

    this->peaks = validPeaks;
    eic->peaks = validPeaks;

    //Issue 759: Retain for grouping
    //delete(eic);
}

//Issue 759
void SECTraceGroups::computePeakGroups(bool debug) {

    //reset to prepare for new computation
    delete_all(samples);
    groups.clear();

    // Prepare EICs
    vector<EIC*> eics{};

    map<mzSample*, SECTrace*> sampleToTrace{};

    unsigned long traceCounter = 0;
    for (SECTrace *trace : secTraces){
        if (trace && trace->eic) {
            string sampleName = to_string(traceCounter);
            if (!trace->id.empty()) {
                sampleName = trace->id;
            }

            mzSample *traceSample = new mzSample();
            traceSample->setSampleId(traceCounter);
            traceSample->setSampleName(sampleName);
            traceSample->isBlank = false;

            //Necessary for EIC::calculateBlankBackground() calculation
            trace->eic->sample = traceSample;

            //Needed for re-assignment of corrected peaks
            sampleToTrace.insert(make_pair(traceSample, trace));

            samples.push_back(traceSample);
            if (debug) cout << "Sample: '" << sampleName << "': ";

            for (Peak& p : trace->eic->peaks) {
                p.sample = traceSample;
                if (debug) cout << "(" << p.rt << ", " << p.peakIntensity << ") ";
            }
            if (debug) cout << endl;

            eics.push_back(trace->eic);

            //peaks will be reassigned after grouping, with possible merges.
            trace->peaks.clear();

            traceCounter++;
        }
    }

    vector<PeakGroup> unfilteredGroups = EIC::groupPeaksE(eics, params->toPeakPickingAndGroupingParams(), debug);

    //SECTrace peaks are updated after grouping, as grouping may have merged peaks
    //A group is only retained if has sufficiently many peaks

    //Issue 781: groups should be re-numbered to only include those groups that pass filters.
    //This step is necessary to keep group IDs aligned
    unsigned int groupNum = 0;

    for (PeakGroup& group : unfilteredGroups) {

        if (group.peakCount() >= params->groupMinNumPeaks) {
            group.groupId = groupNum;

            for (Peak& p : group.peaks) {
                p.groupNum = groupNum;
                if (sampleToTrace.find(p.getSample()) != sampleToTrace.end()) {
                    SECTrace *trace = sampleToTrace[p.getSample()];
                    trace->peaks.push_back(p);

                    if (debug) {
                        cout << "(" << trace->analyteId
                             << ", " << trace->biologicalId
                             << "): group #"
                             << groupNum
                             << endl;
                    }
                }
            }

            groups.push_back(group);
            groupNum++;

        }

    }

    //avoid memory leaks
    delete_all(samples);
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

vector<int> SECTraceGroups::getGroupIdsVector(SECTrace* trace, unsigned long groupIdOffset) {

    if (!trace) return vector<int>{};

    vector<int> groupIds(trace->fractionNums.size(), -1);

    for (unsigned int i = 0; i < trace->peaks.size(); i++) {
        unsigned int max_coord = trace->peaks[i].pos;
        groupIds[max_coord] = trace->peaks[i].groupNum + groupIdOffset;
    }

    return groupIds;
}

vector<float> SECTraceGroups::getGroupsQuantVector(SECTrace *trace, string quantType) {
    if (!trace) return vector<float>{};

    vector<float> groupsQuant(trace->fractionNums.size(), -1);

    for (unsigned int i = 0; i < trace->peaks.size(); i++) {
        unsigned int coord = trace->peaks[i].pos;

        if (quantType == "abundance") {
            groupsQuant[coord] = trace->peaks[i].peakIntensity;
        } else if (quantType == "smoothed") {
            groupsQuant[coord] = trace->peaks[i].smoothedIntensity;
        } else if (quantType == "area") {
            groupsQuant[coord] = trace->peaks[i].peakArea;
        } else if (quantType == "smoothed_area") {
            groupsQuant[coord] = trace->peaks[i].smoothedPeakArea;
        } else {
            groupsQuant[coord] = trace->peaks[i].getQuantByName(quantType);
        }
    }

    return groupsQuant;
}

//Issue 740: return peak positions based on fraction nums.
string SECTrace::getPeakPositionsString(){
    string peakStr = "";

    for (unsigned int i = 0; i < peaks.size(); i++) {
        if (i>0) {
            peakStr = peakStr + ", ";
        }
        peakStr = peakStr + to_string(fractionNums[peaks[i].pos]);
    }

    return peakStr;
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

    float maxMzDiff = params->similarityFractionDiffTol + 0.001f; // avoid rounding errors

    this->ranks = Fragment::findFragPairsGreedyMz(f1, f2, maxMzDiff);

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

int SECTracePeakComparison::getMinFractionNum() {
    return min(first.getMinFractionNum(), second.getMinFractionNum());
}

int SECTracePeakComparison::getMaxFractionNum() {
    return max(first.getMaxFractionNum(), second.getMaxFractionNum());
}

SECTracePeakComparison::SECTracePeakComparison(SECTrace *firstTrace, int firstPeakNum, SECTrace *secondTrace, int secondPeakNum, shared_ptr<SECSearchParameters> params){

    //explicitly set this to false initially so that early returns indicate failure to pass parameters.
    isPassesParameterFilters = false;

    this->first = SECTracePeak(firstTrace, firstPeakNum);
    this->second = SECTracePeak(secondTrace, secondPeakNum);

    //this comparison is made prior to the creation of SECTracePeakComparison object
    peakCenterDistance = abs(first.getPeakFractionNum() - second.getPeakFractionNum());

    secFractionOverlap = mzUtils::checkOverlap(
                first.getMinFractionNum(), first.getMaxFractionNum(),
                second.getMinFractionNum(), second.getMaxFractionNum());

    if (secFractionOverlap < params->peakSimMinSecFractionOverlap) return;

    computeSecFractionJaccard();

    if (secFractionJaccard < params->peakSimMinSecFractionJaccard) return;

    computeComparableRange();

    pearsonCorrelationSmoothed = mzUtils::correlation(first.getSmoothedIntensities(comparableRange), second.getSmoothedIntensities(comparableRange));

    if (pearsonCorrelationSmoothed < params->peakSimMinSmoothedCorrelation) return;

    pearsonCorrelationRaw = mzUtils::correlation(first.getRawIntensities(comparableRange), second.getSmoothedIntensities(comparableRange));

    if (pearsonCorrelationRaw < params->peakSimMinRawCorrelation) return;

    //by default, this is set to false, early returns from this constructor fail to set this to true.
    isPassesParameterFilters = true;
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
        if (debug)  cout << "i=" << i << ": " << ithTrace->id << ", " << ithTrace->peaks.size() << " peaks." << endl;

        if (ithTrace->peaks.empty()) continue;

        for (unsigned int j = i+1; j < traces.size(); j++) {

            auto jthTrace = traces[j];
            if (debug)  cout << "j=" << j << ": " << jthTrace->id << ", " << jthTrace->peaks.size() << " peaks." << endl;

            if (jthTrace->peaks.empty()) continue;

            for (unsigned int k = 0; k < ithTrace->peaks.size(); k++) {
                Peak peakI = ithTrace->peaks.at(k);

                int ithFractionNum = static_cast<int>(peakI.rt);
                if (params->peakSimMinFractionNumber > 0 && ithFractionNum < params->peakSimMinFractionNumber) continue;
                if (params->peakSimMaxFractionNumber > 0 && ithFractionNum > params->peakSimMaxFractionNumber) continue;

                for (unsigned int l = 0; l < jthTrace->peaks.size(); l++) {
                    Peak peakJ = jthTrace->peaks.at(l);

                    int jthFractionNum = static_cast<int>(peakJ.rt);
                    if (params->peakSimMinFractionNumber > 0 && jthFractionNum < params->peakSimMinFractionNumber) continue;
                    if (params->peakSimMaxFractionNumber > 0 && jthFractionNum > params->peakSimMaxFractionNumber) continue;

                    int fracDiff = static_cast<int>(abs(peakI.rt - peakJ.rt));
                    if (debug) cout << "(i[k], j[l]): " << "(" << i << "[" << k << "], " << j << "[" << l << "]) fracDiff = " << fracDiff;

                    if (fracDiff <= params->peakSimMaxCenterDiff) {

                        if (debug) cout << " (comparison)" << endl;

                        SECTracePeakComparison comparison = SECTracePeakComparison(
                                    ithTrace, static_cast<int>(k),
                                    jthTrace, static_cast<int>(l),
                                    params);

                        if (debug) comparison.printSummary();

                        if (comparison.isPassesParameterFilters) {
                            peakComparisons.push_back(comparison);
                            if (debug) cout << "PASSES peakSim FILTERS" << endl;
                        } else if (debug) {
                            cout << "FAILS peakSim FILTERS" << endl;
                        }

                    } else if (debug) {
                        cout << " (no comparison)" << endl;
                    }
                }
            }

            if (debug) cout << "Completed j=" << j << ": " << jthTrace->id << endl;
        }

        if (debug) cout << "Completed i=" << i << ": " << ithTrace->id << endl;
    }

    if (debug) cout << "Retained " << peakComparisons.size() << " peak comparisons." << endl;

    sort(peakComparisons.begin(), peakComparisons.end(), [](SECTracePeakComparison& lhs, SECTracePeakComparison& rhs){
        if (lhs.pearsonCorrelationSmoothed == rhs.pearsonCorrelationSmoothed) {
            if (lhs.pearsonCorrelationRaw == rhs.pearsonCorrelationRaw) {
                return lhs.getPeakComparisonId() < rhs.getPeakComparisonId();
            } else {
                return lhs.pearsonCorrelationRaw < rhs.pearsonCorrelationRaw;
            }
        } else {
            return lhs.pearsonCorrelationSmoothed < rhs.pearsonCorrelationSmoothed;
        }
    });

    return peakComparisons;
}

/**
 * @brief SECTraceDiffGenerator::generateSECTraceDiffs
 * The SECTrace.id field is used to map reference traces to compare traces.
 *
 * @param referenceTraces
 * @param compareTraces
 * @param debug
 * @return
 */
vector<SECTraceDiff*> SECTraceDiffGenerator::generateSECTraceDiffs(
        vector<SECTrace*> referenceTraces,
        vector<SECTrace*> compareTraces,
        bool debug) {

    vector<SECTraceDiff*> secTraceDiffs{};

    set<string> referenceTraceIds{};
    map<string, SECTrace*> compareTracesMap{};

    //Ensure that each SECTrace ID is seen only once in the vector of compare traces
    for (auto compareTrace : compareTraces) {
        string id = compareTrace->id;
        if (compareTracesMap.find(id) != compareTracesMap.end()) {
            cerr << "Duplicate ID: '" << id << "' in compareTraces. This is illegal, exiting." << endl;
            //Issue 643: avoid aborting for debugging
            //abort();
        }
        compareTracesMap.insert(make_pair(id, compareTrace));
    }

    //Ensure that each SECTrace ID is seen only once in the vector of reference traces
    for (auto refTrace : referenceTraces) {
        string id = refTrace->id;
        if (referenceTraceIds.find(id) != referenceTraceIds.end()) {
            cerr << "Duplicate ID: '" << id << "' in referenceTraces. This is illegal, exiting." << endl;
            //Issue 643: avoid aborting for debugging
            //abort();
        }
        referenceTraceIds.insert(id);

        if (compareTracesMap.find(id) != compareTracesMap.end()) {
            SECTrace *compareTrace = compareTracesMap[id];
            SECTraceDiff *diff = new SECTraceDiff(refTrace, compareTrace, debug);
            secTraceDiffs.push_back(diff);
        }
    }

    return secTraceDiffs;
}
