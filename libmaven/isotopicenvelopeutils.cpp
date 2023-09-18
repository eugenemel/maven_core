#include "isotopicenvelopeutils.h"

//private constructor: this is a singleton
IsotopeProcessorOptions::IsotopeProcessorOptions(){}

IsotopeProcessorOptions& IsotopeProcessorOptions::instance() {
    static IsotopeProcessorOptions options;
    return options;
}

void IsotopeProcessorOptions::setOptions(string config_file) {
    this->config_file = config_file;

    //TODO: import config file
}

shared_ptr<IsotopicExtractionParameters> IsotopeProcessorOptions::getExtractionParameters() {
    shared_ptr<IsotopicExtractionParameters> params = shared_ptr<IsotopicExtractionParameters>(new IsotopicExtractionParameters());

    //TODO: fill in params

    return params;
}

double IsotopicEnvelope::getTotalIntensity() {
    if (totalIntensity < 0) {
        totalIntensity = std::accumulate(intensities.begin(), intensities.end(), 0.0);
    }
    return totalIntensity;
}

void IsotopicEnvelopeGroup::print() {
    if (!group) {
        cout << "Peak group is null." << endl;
        return;
    }

    if (!group->compound) {
        cout << "Compound is null." << endl;
        return;
    }

    if (!group->adduct) {
        cout << "Adduct is null." << endl;
        return;
    }

    cout << group->compound->name << " "
         << group->adduct->name << " (m/z, rt) = ("
         << group->compound->precursorMz << ", "
         << group->medianRt()
         << ")"
         << endl;

    cout << "Isotopes: ";
    cout << "{";
    for (unsigned int i = 0; i < isotopes.size(); i++) {
        Isotope isotope = isotopes.at(i);
        if (i > 0) cout << ", ";
        cout << isotope.name;
    }
    cout << "}\n";

    for (auto it = envelopeBySample.begin(); it != envelopeBySample.end(); ++it) {
        mzSample *sample = it->first;
        IsotopicEnvelope envelope = it->second;
        cout << sample->sampleName << ": ";
        envelope.print();
        cout << "\n";
    }

    cout << endl;
}

IsotopicEnvelope& IsotopicEnvelope::none() {
    static IsotopicEnvelope none;
    none.source = "NONE";
    return none;
}

void IsotopicEnvelope::print() {

    stringstream ss;
    ss << std::fixed << setprecision(4);

    ss << "{";
    for (unsigned int i = 0; i < intensities.size(); i++) {
        if (i > 0) ss << ", ";
        ss << intensities.at(i);
    }
    ss << "}; Total=";
    ss << getTotalIntensity();

    cout << ss.str();
}

string IsotopicExtractionParameters::getAlgorithmName(IsotopicExtractionAlgorithm algorithm) {

    if (algorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS) {
        return "peak-full-rt-bounds";
    } else if (algorithm == IsotopicExtractionAlgorithm::PEAK_SHRINKING_RT_BOUNDS) {
        return "peak-shrinking-rt-bounds";
    }

    return "unknown";
}

IsotopicEnvelope IsotopicEnvelopeExtractor::extractEnvelope(mzSample *sample, Peak *peak, vector<Isotope> &isotopes, shared_ptr<IsotopicExtractionParameters> params) {
    IsotopicEnvelope envelope;

    if (params->algorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS) {
        envelope = extractEnvelopePeakFullRtBounds(sample, peak, isotopes, params);
    } else if (params->algorithm == IsotopicExtractionAlgorithm::PEAK_SHRINKING_RT_BOUNDS) {
        envelope = extractEnvelopePeakShrinkingRtBounds(sample, peak, isotopes, params);
    }

    envelope.source = IsotopicExtractionParameters::getAlgorithmName(params->algorithm);
    return envelope;
}

IsotopicEnvelope IsotopicEnvelopeExtractor::extractEnvelopePeakFullRtBounds(mzSample* sample, Peak *peak, vector<Isotope>& isotopes, shared_ptr<IsotopicExtractionParameters> params){

    IsotopicEnvelope envelope;
    vector<double> intensities(isotopes.size());

    if (peak) {
        for (unsigned int i = 0; i < isotopes.size(); i++) {

            Isotope isotope = isotopes.at(i);

            EIC *eic = sample->getEIC(
                        static_cast<float>(isotope.mz - params->mzTol),
                        static_cast<float>(isotope.mz + params->mzTol),
                        peak->rtmin,
                        peak->rtmax,
                        1);

            intensities.at(i) = std::accumulate(eic->intensity.begin(), eic->intensity.end(), 0.0);
        }

       envelope.intensities = intensities;
       envelope.getTotalIntensity();

    }

    return envelope;
}

IsotopicEnvelope IsotopicEnvelopeExtractor::extractEnvelopePeakShrinkingRtBounds(mzSample* sample, Peak *peak, vector<Isotope>& isotopes, shared_ptr<IsotopicExtractionParameters> params){
    IsotopicEnvelope envelope;

    //TODO

    return envelope;
}

string IsotopicExtractionParameters::encodeParams() {
    // TODO
    return "";
}

shared_ptr<IsotopicExtractionParameters> IsotopicExtractionParameters::decode(string encodedIsotopicExtractionParameters){
    shared_ptr<IsotopicExtractionParameters> params = shared_ptr<IsotopicExtractionParameters>(new IsotopicExtractionParameters());
    //TODO
    return params;
}
