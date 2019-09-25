#pragma once

#include "mzSample.h"
#include "mzUtils.h"

class mzSample;

class DirectInfusionProcessor {

    public:

    explicit DirectInfusionProcessor(const vector<mzSample*>& samples, const vector<Compound*>& compounds){
        this->samples = samples;
        this->compounds = compounds;
    }

    private:
        vector<mzSample*> samples;
        vector<Compound*> compounds;

        void process();
};
