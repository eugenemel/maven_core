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

    //only used by large peptide protein binding assay ion list
    vector<string> adductName{};
    vector<string> isotopeCode{};
    vector<double> isotopeNaturalAbundance{};
    vector<int> numBound{};
    vector<bool> isMaxAbundanceIsotope{};
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
        string molecularFormula,
        string peptideSequence,
        vector<Adduct>& adducts,
        double boundLigandExactMass,
        int minNumBoundLigand,
        int maxNumBoundLigand,
        double peptidePredictedIsotopeRatioThreshold,
        bool ms1UseDaTol,
        double ms1PpmTol,
        double ms1DaTol,
        bool debug=false
        );
};

/**
 * @brief The MaldesiParameters class
 * dedicated class for maldesi parameters, primarily used when
 * encoding/decoding from a C++ context, e.g. mzkitcpp
 */
class MaldesiParameters {
public:
    int minNumBoundLigand = -1;
    int maxNumBoundLigand = -1;
    double boundLigandExactMass = -1.0;
    vector<Adduct> adducts{};

    //single encoded MaldesiParameters object
    //e.g. {minNumBoundLigand=1;maxNumBoundLigand=2;boundLigandExactMass=100;}
    static shared_ptr<MaldesiParameters> decode(const string& encodedParams);

    //series of scan-specific MaldesiParameters objects
    //e.g. {1={minNumBoundLigand=1;maxNumBoundLigand=2;boundLigandExactMass=100;}, 2={minNumBoundLigand=1;adducts={[M+H]+&[M+2H]2+&}}}
    static map<int, shared_ptr<MaldesiParameters>> decodeScanSpecific(const string& encodedParams);
};

#endif // MALDESI_H
