#include "mzSample.h"

//default constructor
PeptideStabilitySearchParameters::PeptideStabilitySearchParameters() {

    //Enforce synchrony between two shared_ptr<PeakPickingAndGroupingParameters>
    this->isotopeParameters.peakPickingAndGroupingParameters=this->peakPickingAndGroupingParameters;
}

string PeptideStabilitySearchParameters::encodeParams() {

    string encodedParams = baseParams();

    encodedParams = encodedParams + peakPickingAndGroupingParameters->getEncodedPeakParameters();

    encodedParams = encodedParams + "peptideName=" + peptideName + ";";
    encodedParams = encodedParams + "peptideExactMass=" + to_string(peptideExactMass) + ";";
    encodedParams = encodedParams + "peptideAdducts=" + peptideAdducts + ";";
    encodedParams = encodedParams + "peptideRtTol=" + to_string(peptideRtTol) + ";";

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
    if (decodedMap.find("peptideName") != decodedMap.end()) {
        params->peptideName = decodedMap["peptideName"];
    }
    if (decodedMap.find("peptideExactMass") != decodedMap.end()) {
        params->peptideExactMass = stod(decodedMap["peptideExactMass"]);
    }
    if (decodedMap.find("peptideAdducts") != decodedMap.end()) {
        params->peptideAdducts = decodedMap["peptideAdducts"];
    }
    if (decodedMap.find("peptideRtTol") != decodedMap.end()) {
        params->peptideRtTol = stof(decodedMap["peptideRtTol"]);
    }

    if (decodedMap.find("isPullIsotopes") != decodedMap.end()) {
        params->isPullIsotopes = decodedMap["isPullIsotopes"] == "1";
    }

    if (decodedMap.find("minNumIsotopes") != decodedMap.end()) {
        params->minNumIsotopes = stoi(decodedMap["minNumIsotopes"]);
    }

    // END PeptideStabilitySearchParameters

    return params;
}

vector<PeakGroup> PeptideStabilityProcessor::processPeptideStabilitySet(
    vector<PeakGroup>& peptideStabilityGroupSet,
    const vector<mzSample*>& samples,
    shared_ptr<PeptideStabilitySearchParameters> params,
    bool debug
    ) {

    vector<PeakGroup> groups{};

    for (PeakGroup g : peptideStabilityGroupSet) {
        if (params->isPullIsotopes) {
            g.pullIsotopes(
                params->isotopeParameters,
                samples,
                debug);
        }

        if (g.getChildren().size() < params->minNumIsotopes) continue;

        // any peak groups that survive all filters are saved.
        // peak groups may be transformed along the way (for example, isotopes may be added).
        groups.push_back(g);
    }

    return groups;
}

void PeptideStabilityProcessor::labelPeptideStabilitySet(
    vector<PeakGroup>& peptideStabilityGroupSet,
    shared_ptr<PeptideStabilitySearchParameters> params,
    bool debug
    ) {

    //sort peak groups in descending order by intensity
    sort(peptideStabilityGroupSet.begin(), peptideStabilityGroupSet.end(), [](PeakGroup & lhs, PeakGroup & rhs) {
        return lhs.maxIntensity > rhs.maxIntensity;
    });

    float peptideRt = -1.0f;
    set<string> labeledAdducts{};

    //Associate peptides with chosen candidate peptide, based on proximity to max intensity peptide.
    for (PeakGroup & g : peptideStabilityGroupSet) {

        // block is only triggered on the first peak group (most abundant peak group)
        if (peptideRt < 0) {
            peptideRt = g.maxPeakRtVal;
            labeledAdducts.insert(g.compound->adductString);
            g.addLabel('e');
        } else {
            // compare all subsequent groups to data defined from first group
            // candidate peak group must be in RT range
            float rtDiff = abs(g.maxPeakRtVal - peptideRt);
            string adductName = g.compound->adductString;
            if (rtDiff < params->peptideRtTol && labeledAdducts.find(adductName) == labeledAdducts.end()) {
                labeledAdducts.insert(adductName);
                g.addLabel('e');
            }
        }
    }
}
