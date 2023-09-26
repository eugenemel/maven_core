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

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopes(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<mzSample*> samples,
    IsotopeParameters params,
    bool debug){

    vector<Isotope> theoreticalIsotopes = MassCalculator::computeIsotopes(
        compound->formula,
        adduct,
        params.maxIsotopesToExtract,
        params.isC13Labeled,
        params.isN15Labeled,
        params.isS34Labeled,
        params.isD2Labeled
        );

    vector<Isotope> isotopes = theoreticalIsotopes;
    if (params.isCondenseTheoreticalIsotopes) {
        isotopes = IsotopicEnvelopeAdjuster::condenseTheoreticalIsotopes(
            theoreticalIsotopes,
            params,
            false
            );
    }

    IsotopicEnvelopeGroup envelopeGroup;

    if (params.isotopicExtractionAlgorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA) {
        envelopeGroup = extractEnvelopesPeakFullRtBounds(compound, adduct, group, isotopes, params, debug);
    } else if (params.isotopicExtractionAlgorithm == IsotopicExtractionAlgorithm::PEAK_SHRINKING_RT_BOUNDS_AREA) {
        //TODO
    }

    envelopeGroup.extractionAlgorithmName = IsotopeParameters::getAlgorithmName(params.isotopicExtractionAlgorithm);
    return envelopeGroup;

}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopesPeakFullRtBounds(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<Isotope>& isotopes,
    IsotopeParameters params,
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
        g.setType(PeakGroup::GroupType::IsotopeType);

        envelopeGroup.isotopePeakGroups[i] = g;
    }

    for (auto & peak : group->peaks) {

        IsotopicEnvelope envelope;
        envelope.intensities = vector<double>(isotopes.size());

        auto sample = peak.sample;

        for (unsigned int i = 0; i < isotopes.size(); i++) {

            Isotope isotope = isotopes.at(i);


            EIC *eic = sample->getEIC(
                static_cast<float>(isotope.mz - isotope.mz/1e6f*params.ppm),
                static_cast<float>(isotope.mz + isotope.mz/1e6f*params.ppm),
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

            //RT related
            p.scan = peak.scan;
            p.minscan = peak.minscan;
            p.maxscan = peak.maxscan;
            p.pos = peak.scan;
            p.minpos = peak.minpos;
            p.maxpos = peak.maxpos;
            p.width = (peak.maxpos-peak.minpos);
            p.rtmin = peak.rtmin;
            p.rtmax = peak.rtmax;
            p.rt = peak.rt;

            //m/z related
            p.mzmin = mzmin;
            p.mzmax = mzmax;
            p.peakMz = isotope.mz;
            p.baseMz = isotope.mz;
            p.medianMz = isotope.mz;

            //intensity related
            p.peakArea = intensity;
            //TODO: use peakAreaCorrected as natural-abundance corrected abundance

            //avoid writing junk into mzrollDB
            p.gaussFitR2 = 0;
            p.noNoiseObs = 0;
            p.signalBaselineRatio = 0;

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
        IsotopeParameters params,
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
    IsotopeParameters params,
    bool debug){

    IsotopicEnvelopeGroup envelopeGroup;

    //TODO

    return envelopeGroup;
}

vector<Isotope> IsotopicEnvelopeAdjuster::condenseTheoreticalIsotopes(
    vector<Isotope> defaultIsotopes,
    IsotopeParameters params,
    bool debug){

    //TODO: implement
    return defaultIsotopes;
}
