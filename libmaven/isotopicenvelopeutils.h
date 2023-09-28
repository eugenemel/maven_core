#pragma once

#ifndef ISOTOPICENVELOPEUTILS_H
#define ISOTOPICENVELOPEUTILS_H

#include "mzSample.h"
#include <numeric>

/**
 * @brief The IsotopicEnvelope class
 * Measurements for a set of isotopes in a single sample.
 * the corresponding isotope for measurements in the intensities field
 * should be stored in a separate IsotopicEnvelopeGroup class.
 */
class IsotopicEnvelope {
public:

    static IsotopicEnvelope& none();

    string source = "unknown";
    double totalIntensity = -1.0;
    vector<double> intensities;

    double getTotalIntensity();
    void print();
};

/**
 * @brief The IsotopicEnvelopeGroup class
 * All measurements must follow the isotopic pattern defined in the
 * isotopes field.
 * Note that this means that all sample measurements must match this pattern -
 * cannot have some samples with one set of isotopes and other samples with a different set
 * of isotopes.
 *
 * This supports that isotopes can be a compound-specific thing - but not a sample specific thing.
 * In general, all samples should ve generated with the same resolving power,
 * so you shouldn't expect different sets of isotopes for different samples.
 */
class IsotopicEnvelopeGroup {
public:

    /** Input fields **/

    Compound *compound = nullptr;
    Adduct *adduct = nullptr;
    PeakGroup *group = nullptr;
    vector<Isotope> isotopes{};

    /** Computed Values **/

    //View cross-isotope measurements by sample
    //Individual (Sample, Isotope) measurements are stored by elements
    //in the IsotopicEnvelope.intensities vector (same indices as IsotopicEnvelopeGroup.isotopes).
    map<mzSample*, IsotopicEnvelope> envelopeBySample{};

    //View cross-sample measurements by isotope (index matches IsotopicEnvelopeGroup.isotopes)
    //Individual (Isotope, Sample) measurements are stored as peaks in the PeakGroup
    vector<PeakGroup> isotopePeakGroups{};

    //Usually provided by IsotopicExtractionParameters
    string extractionAlgorithmName = "";

    //convenience method that sets children of the IsotopicEnvelope.group field
    //to the IsotopicEnvelopeGroup.children (same indices as IsotopicEnvelopeGroup.isotopes).
    void setIsotopesToChildrenPeakGroups(Classifier *classifier = nullptr);

    void print();
};

class IsotopicEnvelopeExtractor {
public:

    /**
     * @brief extractEnvelope
     * Most common entry point. The specific extraction algorithm will be selected based on the input params.
     *
     * Note that sample information is passed in - this supports cases where a sample of interest has no Peak in
     * inthe input PeakGroup.
     *
     * Isotopes are computed based on the input parameters and compound information.
     *
     * @param compound
     * @param adduct
     * @param group
     * @param samples
     * @param params
     * @param debug
     * @return
     */
    static IsotopicEnvelopeGroup extractEnvelopes(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<mzSample*>& samples,
        IsotopeParameters params,
        bool debug=false);

    //usually called from IsotopicEnvelopeExtractor::extractEnvelope()
    static IsotopicEnvelopeGroup extractEnvelopesPeakFullRtBounds(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        IsotopeParameters params,
        bool debug=false);

    //TODO: not yet implemented
    static IsotopicEnvelopeGroup extractEnvelopesPeakShrinkingRtBounds(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        IsotopeParameters params,
        bool debug=false);

    //The original approach implemented in MAVEN up through version 2.0 2023-09-25
    static IsotopicEnvelopeGroup extractEnvelopesVersion1(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        IsotopeParameters params,
        bool debug=false);
};

class IsotopicEnvelopeAdjuster {
public:

    /**
     * @brief condenseTheoreticalIsotopes
     * Based on the resolving power of the instrument and other parameters, condense a vector
     * of isotopes produced by mzMassCalculator::computeIsotopes() into a new vector of isotopes.
     *
     * For example, for low-res data, condense separate 13C and 15N isotopes into a single "[M+1]" isotopes.
     *
     * @param uncondensedIsotopes
     * @param params
     * @param debug
     * @return corrected vector<Isotope> with new isotopic names, theretical m/z, and theoretical abundances.
     */
    static vector<Isotope> condenseTheoreticalIsotopes(
        vector<Isotope> uncondensedIsotopes,
        IsotopeParameters params,
        bool debug
        );
};


#endif // ISOTOPICENVELOPEUTILS_H
