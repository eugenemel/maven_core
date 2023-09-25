#include "isotopicenvelopeutils.h"

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

    if (algorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA) {
        return "PEAK_FULL_RT_BOUNDS_AREA";
    } else if (algorithm == IsotopicExtractionAlgorithm::PEAK_SHRINKING_RT_BOUNDS_AREA) {
        return "PEAK_SHRINKING_RT_BOUNDS_AREA";
    }

    return "unknown";
}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopes(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<Isotope>& isotopes,
    vector<mzSample*> samples,
    shared_ptr<IsotopicExtractionParameters> params,
    bool debug){

    IsotopicEnvelopeGroup envelopeGroup;

    if (params->algorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA) {
        envelopeGroup = extractEnvelopesPeakFullRtBounds(compound, adduct, group, isotopes, params, debug);
    } else if (params->algorithm == IsotopicExtractionAlgorithm::PEAK_SHRINKING_RT_BOUNDS_AREA) {
        //TODO
    }

    envelopeGroup.extractionAlgorithmName = IsotopicExtractionParameters::getAlgorithmName(params->algorithm);
    return envelopeGroup;
}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopesPeakFullRtBounds(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<Isotope>& isotopes,
    shared_ptr<IsotopicExtractionParameters> params,
    bool debug) {

    IsotopicEnvelopeGroup envelopeGroup;

    envelopeGroup.compound = compound;
    envelopeGroup.adduct = adduct;
    envelopeGroup.group = group;
    envelopeGroup.isotopes = isotopes;

    //initialize peak groups
    envelopeGroup.isotopePeakGroups = vector<PeakGroup>(isotopes.size());
    for (unsigned int i = 0; i < envelopeGroup.isotopePeakGroups.size(); i++) {
        Isotope isotope = isotopes[i];

        PeakGroup g;
        g.meanMz = isotope.mz;
        g.tagString = isotope.name;
        g.expectedAbundance = isotope.abundance;
        g.isotopeC13count= isotope.C13;
        g.isotopicIndex = i;
        g.compound = compound;
        g.adduct = adduct;

        envelopeGroup.isotopePeakGroups[i] = g;
    }

    for (auto & peak : group->peaks) {

        IsotopicEnvelope envelope;
        envelope.intensities = vector<double>(isotopes.size());

        auto sample = peak.sample;

        for (unsigned int i = 0; i < isotopes.size(); i++) {

            Isotope isotope = isotopes.at(i);

            EIC *eic = sample->getEIC(
                static_cast<float>(isotope.mz - params->isotopicTheoreticalMzTolerance),
                static_cast<float>(isotope.mz + params->isotopicTheoreticalMzTolerance),
                peak.rtmin,
                peak.rtmax,
                1);

            float mzmin = eic->mzmin;
            float mzmax = eic->mzmax;

            double intensity = std::accumulate(eic->intensity.begin(), eic->intensity.end(), 0.0);

            delete(eic);
            eic = nullptr;

            envelope.intensities.at(i) = intensity;

            Peak p;
            p.scan = peak.scan;
            p.minscan = peak.minscan;
            p.maxscan = peak.maxscan;
            p.pos = peak.scan;
            p.minpos = peak.minpos;
            p.maxpos = peak.maxpos;
            p.width = (peak.maxpos-peak.minpos);

            p.rtmin = peak.rtmin;
            p.rtmax = peak.rtmax;
            p.mzmin = mzmin;
            p.mzmax = mzmax;

            p.peakArea = intensity;

            p.sample = sample;

            envelopeGroup.isotopePeakGroups[i].addPeak(p);
        }

        envelope.getTotalIntensity();

        envelopeGroup.envelopeBySample.insert(make_pair(sample, envelope));
    }

    return envelopeGroup;
}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopesVersion1(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        shared_ptr<IsotopicExtractionParameters> params,
        bool debug) {

    IsotopicEnvelopeGroup envelopeGroup;

    //TODO

    return envelopeGroup;
}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopesPeakShrinkingRtBounds(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<Isotope>& isotopes,
    shared_ptr<IsotopicExtractionParameters> params,
    bool debug){

    IsotopicEnvelopeGroup envelopeGroup;

    //TODO

    return envelopeGroup;
}

string IsotopicExtractionParameters::encodeParams() {
    string encodedParams;

    //extraction algorithm
    string algorithmStr = "UNKNOWN";
    if (algorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA) {
       algorithmStr = "PEAK_FULL_RT_BOUNDS_AREA";
    } else if (algorithm == IsotopicExtractionAlgorithm::PEAK_SHRINKING_RT_BOUNDS_AREA) {
       algorithmStr = "PEAK_SHRINKING_RT_BOUNDS_AREA";
    }
    encodedParams = encodedParams + "algorithm" + "=" + algorithmStr + ";";

    //theoretical isotope handling
    encodedParams = encodedParams + "isotopicTheoreticalMzTolerance" + "="+ to_string(isotopicTheoreticalMzTolerance) + ";";
    string isotopicTheoreticalMzToleranceTypeStr = "UNKNOWN";
    if (isotopicTheoreticalMzToleranceType == IsotopicTheoreticalMzToleranceType::Da) {
       isotopicTheoreticalMzToleranceTypeStr = "Da";
    } else if (isotopicTheoreticalMzToleranceType == IsotopicTheoreticalMzToleranceType::ppm) {
       isotopicTheoreticalMzToleranceTypeStr = "ppm";
    }
    encodedParams = encodedParams + "isotopicTheoreticalMzToleranceType" + "=" + isotopicTheoreticalMzToleranceTypeStr + ";";

    return encodedParams;
}

shared_ptr<IsotopicExtractionParameters> IsotopicExtractionParameters::decode(string encodedParams){
    shared_ptr<IsotopicExtractionParameters> params = shared_ptr<IsotopicExtractionParameters>(new IsotopicExtractionParameters());

    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams);

    //extraction algorithm
    if (decodedMap.find("algorithm") != decodedMap.end()) {
       string algorithmStr = decodedMap["algorithm"];
       if (algorithmStr == "PEAK_FULL_RT_BOUNDS_AREA") {
            params->algorithm = IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA;
       } else if (algorithmStr == "PEAK_SHRINKING_RT_BOUNDS_AREA") {
            params->algorithm = IsotopicExtractionAlgorithm::PEAK_SHRINKING_RT_BOUNDS_AREA;
       }
    }

    //theoretical isotope handling
    if (decodedMap.find("isotopicTheoreticalMzTolerance") != decodedMap.end()) {
       params->isotopicTheoreticalMzTolerance = stod(decodedMap["isotopicTheoreticalMzTolerance"]);
    }
    if (decodedMap.find("isotopicTheoreticalMzToleranceType") != decodedMap.end()) {
       string isotopicTheoreticalMzToleranceTypeStr = decodedMap["isotopicTheoreticalMzToleranceType"];
       if (isotopicTheoreticalMzToleranceTypeStr == "Da") {
            params->isotopicTheoreticalMzToleranceType = IsotopicTheoreticalMzToleranceType::Da;
       } else if (isotopicTheoreticalMzToleranceTypeStr == "ppm") {
            params->isotopicTheoreticalMzToleranceType = IsotopicTheoreticalMzToleranceType::ppm;
       }
    }

    return params;
}

vector<Isotope> IsotopicEnvelopeAdjuster::condenseTheoreticalIsotopes(
    vector<Isotope> defaultIsotopes,
    shared_ptr<IsotopicExtractionParameters> params,
    bool debug){

    //TODO: implement
    return defaultIsotopes;
}
