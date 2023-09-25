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

    //View cross-sample measurements by isotope
    //Individual (Isotope, Sample) measurements are stored as peaks in the PeakGroup
    map<Isotope, PeakGroup> peakGroupByIsotope{};

    //Usually provided by IsotopicExtractionParameters
    string extractionAlgorithmName;

    //convenience method that sets children of the IsotopicEnvelope.group field
    //to the IsotopicEnvelopeGroup.children (same indices as IsotopicEnvelopeGroup.isotopes).
    void setIsotopesToChildrenPeakGroups();

    void print();
};


enum IsotopicTheoreticalMzToleranceType{Da, ppm};
enum IsotopicExtractionAlgorithm{PEAK_FULL_RT_BOUNDS_AREA, PEAK_SHRINKING_RT_BOUNDS_AREA};

/**
 * @brief The IsotopicExtractionParameters class
 * Container class to hold all parameters associated with isotopic extraction.
 * Distinct from IsotopeProcessorOptions, which is specific to isotopeprocessor.
 * This set of parameters is also applicable to the MAVEN GUI.
 */
class IsotopicExtractionParameters {
public:

    IsotopicExtractionAlgorithm algorithm = IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA;

    double isotopicTheoreticalMzTolerance = 0.01;
    IsotopicTheoreticalMzToleranceType isotopicTheoreticalMzToleranceType = IsotopicTheoreticalMzToleranceType::Da;

    string encodeParams();
    static shared_ptr<IsotopicExtractionParameters> decode(string encodedIsotopicExtractionParameters);
    static string getAlgorithmName(IsotopicExtractionAlgorithm algorithm);
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
     * @param compound
     * @param adduct
     * @param group
     * @param isotopes
     * @param samples
     * @param params
     * @param debug
     * @return
     */
    static IsotopicEnvelopeGroup extractEnvelopes(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        vector<mzSample*> samples,
        shared_ptr<IsotopicExtractionParameters> params,
        bool debug=false);

    //usually called from IsotopicEnvelopeExtractor::extractEnvelope()
    static IsotopicEnvelopeGroup extractEnvelopesPeakFullRtBounds(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        shared_ptr<IsotopicExtractionParameters> params,
        bool debug=false);

    //TODO: not yet implemented
    static IsotopicEnvelopeGroup extractEnvelopesPeakShrinkingRtBounds(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        shared_ptr<IsotopicExtractionParameters> params,
        bool debug=false);

    //The original approach implemented in MAVEN up through version 2.0 2023-09-25
    static IsotopicEnvelopeGroup extractEnvelopesVersion1(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        shared_ptr<IsotopicExtractionParameters> params,
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
        shared_ptr<IsotopicExtractionParameters> params,
        bool debug
        );
};


#endif // ISOTOPICENVELOPEUTILS_H
