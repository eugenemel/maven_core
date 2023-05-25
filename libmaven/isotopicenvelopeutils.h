#pragma once

#ifndef ISOTOPICENVELOPEUTILS_H
#define ISOTOPICENVELOPEUTILS_H

#include "mzSample.h"

class IsotopeProcessorOptions {

public:
    static IsotopeProcessorOptions& instance();

    void setOption(string key, string value);

    void printOptions();

private:
    //Singleton: constructor is private
    IsotopeProcessorOptions();

    //Singleton: Disable copy constructor
    IsotopeProcessorOptions(const IsotopeProcessorOptions&) = delete;

    //Singleton: Disable assignment operator
    IsotopeProcessorOptions& operator= (const IsotopeProcessorOptions&) = delete;
};

#endif // ISOTOPICENVELOPEUTILS_H
