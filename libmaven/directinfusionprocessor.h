#pragma once

#include "mzSample.h"
#include "mzUtils.h"

class mzSample;

class DirectInfusionProcessor {

    public:
        static void processSingleSample(mzSample* sample, const vector<Compound*>& compounds);
};
