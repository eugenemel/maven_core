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
    encodedParams = encodedParams + "traceMinPeakSN" + "=" + to_string(traceMinPeakSN);
    encodedParams = encodedParams + "traceBaselineDropTopX" + "=" + to_string(traceBaselineDropTopX);
    encodedParams = encodedParams + "tracePeakBoundsMaxIntensityFraction" + "=" + to_string(tracePeakBoundsMaxIntensityFraction);

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
                   shared_ptr<SECSearchParameters> params) {
    //TODO
}
