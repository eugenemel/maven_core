#pragma once

#ifndef ISOTOPICENVELOPEUTILS_H
#define ISOTOPICENVELOPEUTILS_H

#include "mzSample.h"
#include <numeric>

class IsotopeProcessorOptions {

public:
    static IsotopeProcessorOptions& instance();

    void setOptions(string config_file);

    void printOptions();

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

#endif // ISOTOPICENVELOPEUTILS_H
