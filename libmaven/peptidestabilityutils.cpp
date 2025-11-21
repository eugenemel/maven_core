#include "mzSample.h"

//default constructor
PeptideStabilitySearchParameters::PeptideStabilitySearchParameters() {

    //Enforce synchrony between two shared_ptr<PeakPickingAndGroupingParameters>
    this->isotopeParameters.peakPickingAndGroupingParameters=this->peakPickingAndGroupingParameters;
}

string PeptideStabilitySearchParameters::encodeParams() {

    string encodedParams = baseParams();

    encodedParams = encodedParams + peakPickingAndGroupingParameters->getEncodedPeakParameters();

    encodedParams = encodedParams + "isPullIsotopes=" + to_string(isPullIsotopes) + ";";
    encodedParams = encodedParams + "minNumIsotopes=" + to_string(minNumIsotopes) + ";";

    //exclude peak picking and grouping parameters from isotope parameters encoding
    //only using the isotope-specific information
    encodedParams = encodedParams + isotopeParameters.encodeParams(false);

    return encodedParams;
}

shared_ptr<PeptideStabilitySearchParameters> PeptideStabilitySearchParameters::decode(string encodedParams){
    shared_ptr<PeptideStabilitySearchParameters> params = shared_ptr<PeptideStabilitySearchParameters>(new PeptideStabilitySearchParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams); //use semicolon (default)

    params->fillInBaseParams(decodedMap);

    params->peakPickingAndGroupingParameters = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());
    params->peakPickingAndGroupingParameters->fillInPeakParameters(decodedMap);

    IsotopeParameters isotopeParameters = IsotopeParameters::decode(encodedParams);
    params->isotopeParameters = isotopeParameters;

    //enforce synchrony between two copies of PeakPickingAndGroupingParameters
    params->isotopeParameters.peakPickingAndGroupingParameters->fillInPeakParameters(decodedMap);

    // START PeptideStabilitySearchParameters

    if (decodedMap.find("isPullIsotopes") != decodedMap.end()) {
        params->isPullIsotopes = decodedMap["isPullIsotopes"] == "1";
    }

    if (decodedMap.find("minNumIsotopes") != decodedMap.end()) {
        params->minNumIsotopes = stoi(decodedMap["minNumIsotopes"]);
    }

    // END PeptideStabilitySearchParameters

    return params;
}

vector<PeakGroup> PeptideStabilityProcessor::filterPeptideStabilitySet(
    vector<PeakGroup>& peptideStabilityGroupSet,
    shared_ptr<PeptideStabilitySearchParameters> params,
    bool debug
    ) {

    vector<PeakGroup> groups{};

    // TODO

    return groups;
}

void labelPeptideStabilitySet(
    vector<PeakGroup>& peptideStabilityGroupSet,
    shared_ptr<PeptideStabilitySearchParameters> params,
    bool debug
    ) {
    // TODO
}
