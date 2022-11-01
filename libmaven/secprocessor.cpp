#include "secprocessor.h"

string SECSearchParameters::encodeParams() {

    string encodedParams;

    //program level
    encodedParams = encodedParams + "searchVersion" + "=" + searchVersion + ";";

    //SEC trace filling
    encodedParams = encodedParams + "traceMissingIntensityFill" + "=" + to_string(traceMissingIntensityFill) + ";";
    encodedParams = encodedParams + "traceMinFractionNumber" + "=" + to_string(traceMinFractionNumber) + ";";
    encodedParams = encodedParams + "traceMaxFractionNumber" + "=" + to_string(traceMaxFractionNumber) + ";";

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
    encodedParams = encodedParams + "traceMinPeakSN" + "=" + to_string(traceMinPeakSN) + ";";
    encodedParams = encodedParams + "traceBaselineDropTopX" + "=" + to_string(traceBaselineDropTopX) + ";";
    encodedParams = encodedParams + "tracePeakBoundsMaxIntensityFraction" + "=" + to_string(tracePeakBoundsMaxIntensityFraction) + ";";

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
    if (decodedMap.find("traceMinPeakSN") != decodedMap.end()) {
        secSearchParameters->traceMinPeakSN = stof(decodedMap["traceMinPeakSN"]);
    }
    if (decodedMap.find("traceBaselineDropTopX") != decodedMap.end()) {
        secSearchParameters->traceBaselineDropTopX = stoi(decodedMap["traceBaselineDropTopX"]);
    }
    if (decodedMap.find("tracePeakBoundsMaxIntensityFraction") != decodedMap.end()) {
        secSearchParameters->tracePeakBoundsMaxIntensityFraction = stof(decodedMap["tracePeakBoundsMaxIntensityFraction"]);
    }

    return secSearchParameters;
}

SECTrace::SECTrace(SECTraceType type,
                   vector<int> fractionNums,
                   vector<float> rawIntensities,
                   shared_ptr<SECSearchParameters> params,
                   bool debug) {

    if (fractionNums.size() != rawIntensities.size()) {
        cerr << "SECTrace() requires same length for fractionNums, rawIntensities. Exiting." << endl;
        abort();
    }

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

    EIC *eic = new EIC();

    eic->intensity = this->rawIntensities;
    eic->rt = pseudoRt;
    eic->mz = vector<float>(static_cast<unsigned long>(N));
    eic->scannum = this->fractionNums;

    eic->setSmootherType(params->traceSmoothingType);
    eic->setBaselineSmoothingWindow(params->traceWindowSize);
    eic->setBaselineDropTopX(params->traceBaselineDropTopX);

    eic->getPeakPositionsC(params->traceWindowSize, debug, true, params->tracePeakBoundsMaxIntensityFraction);

    this->smoothedIntensities = eic->spline;

    for (auto p : eic->peaks) {
        if (this->rawIntensities[p.pos] >= params->traceMinPeakIntensity && p.signalBaselineRatio >= params->traceMinPeakSN) {
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

        unsigned int left_coord = static_cast<unsigned int>(p.rtmin - params->traceMinFractionNumber + 0.00001f);
        unsigned int max_coord = static_cast<unsigned int>(p.rt - params->traceMinFractionNumber + 0.00001f);
        unsigned int right_coord = static_cast<unsigned int>(p.rtmax - params->traceMinFractionNumber + 0.00001f);

        if (summaryString.at(left_coord) != "") {
            summaryString[left_coord] = summaryString[left_coord] + ", " + left;
        } else {
            summaryString[left_coord] = left;
        }

        if (summaryString.at(max_coord) != "") {
            summaryString[max_coord] = summaryString[max_coord] + ", " + max;
        } else {
            summaryString[max_coord] = max;
        }

        if (summaryString.at(right_coord) != "") {
            summaryString[right_coord] = summaryString[right_coord] + ", " + right;
        } else {
            summaryString[right_coord] = right;
        }
    }

    return summaryString;
}
