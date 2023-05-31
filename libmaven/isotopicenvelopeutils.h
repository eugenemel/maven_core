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

class IsotopicEnvelope {
public:

    string source = "unknown";
    vector<Isotope> isotopes{};
    PeakGroup *group = nullptr;

    double totalIntensity = -1.0;
    vector<double> isotopeIntensity;

    double getIntensity();
    void print();
};

#endif // ISOTOPICENVELOPEUTILS_H
