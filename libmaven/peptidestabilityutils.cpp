#include "mzSample.h"

string PeptideStabilitySearchParameters::encodeParams() {
    return "TODO";
}

shared_ptr<PeptideStabilitySearchParameters> PeptideStabilitySearchParameters::decode(string encodedParams){
    shared_ptr<PeptideStabilitySearchParameters> params = shared_ptr<PeptideStabilitySearchParameters>(new PeptideStabilitySearchParameters());

    // TODO

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
