#pragma once

#ifndef MALDESI_H
#define MALDESI_H

#include "mzSample.h"

class MaldesiIonList {
public:
    vector<double> searchMz{};
    vector<double> minMz{};
    vector<double> maxMz{};
    vector<string> ionName{};
    vector<int> chg{};
    bool isSinglePrecursorMz = false;
};

class MaldesiIonListGenerator {
public:
    static MaldesiIonList getMaldesiIonList(
        double compoundMz,
        vector<Adduct>& adducts,
        bool ms1UseDaTol,
        double ms1PpmTol,
        double ms1DaTol,
        bool debug=false
        );

    static MaldesiIonList getLargePeptideProteinBindingAssayIonList(
        string peptideSequence,
        vector<Adduct>& adducts,
        double boundLigandExactMass,
        int maxNumBoundLigand,
        double peptidePredictedIsotopeRatioThreshold,
        bool ms1UseDaTol,
        double ms1PpmTol,
        double ms1DaTol,
        bool debug=false
        );
};


#endif // MALDESI_H
