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

    Compound *compound = nullptr;
    Adduct *adduct = nullptr;
    PeakGroup *group = nullptr;

    vector<Isotope> isotopes{};

    map<mzSample*, IsotopicEnvelope> envelopeBySample{};

    //TODO: implement these as needed

//    double getQuant(mzSample* sample, Isotope isotope);
//    double getQuant(mzSample* sample, int isotopeIndex);

//    IsotopicEnvelope getEnvelope(mzSample* sample);
//    map<mzSample*, double> getIsotopeQuant(Isotope isotope);

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

    //most common entrypoint
    static IsotopicEnvelope extractEnvelope(mzSample* sample, Peak *peak, vector<Isotope>& isotopes, shared_ptr<IsotopicExtractionParameters> params);

    //usually called from IsotopicEnvelopeExtractor::extractEnvelope()
    static IsotopicEnvelope extractEnvelopePeakFullRtBounds(mzSample* sample, Peak *peak, vector<Isotope>& isotopes, shared_ptr<IsotopicExtractionParameters> params);
    static IsotopicEnvelope extractEnvelopePeakShrinkingRtBounds(mzSample* sample, Peak *peak, vector<Isotope>& isotopes, shared_ptr<IsotopicExtractionParameters> params);

    //The original approach implemented in MAVEN up through version 2.0 2023-09-25
    static IsotopicEnvelope extractEnvelopeVersion1(mzSample* sample, Peak *peak, vector<Isotope>& isotopes);
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
