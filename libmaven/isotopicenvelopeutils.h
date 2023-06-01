#pragma once

#ifndef ISOTOPICENVELOPEUTILS_H
#define ISOTOPICENVELOPEUTILS_H

#include "mzSample.h"
#include <numeric>

class IsotopicExtractionParameters;

class IsotopeProcessorOptions {

public:
    static IsotopeProcessorOptions& instance();

    void setOptions(string config_file);

    void printOptions();

    shared_ptr<IsotopicExtractionParameters> getExtractionParameters();

    //fields below this point
    string config_file = "";

private:
    //Singleton: constructor is private
    IsotopeProcessorOptions();

    //Singleton: Disable copy constructor
    IsotopeProcessorOptions(const IsotopeProcessorOptions&) = delete;

    //Singleton: Disable assignment operator
    IsotopeProcessorOptions& operator= (const IsotopeProcessorOptions&) = delete;
};

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


enum IsotopicExtractionAlgorithm{PEAK_FULL_RT_BOUNDS};

/**
 * @brief The IsotopicExtractionParameters class
 * Container class to hold all parameters associated with isotopic extraction.
 * Distinct from IsotopeProcessorOptions, which is specific to isotopeprocessor.
 * This set of parameters is also applicable to the MAVEN GUI.
 */
class IsotopicExtractionParameters {
public:

    IsotopicExtractionAlgorithm algorithm = IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS;

    double mzTol = 0.01;

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

};


#endif // ISOTOPICENVELOPEUTILS_H
