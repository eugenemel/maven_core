#pragma once

#ifndef MALDESI_H
#define MALDESI_H

#include "mzSample.h"

class MaldesiIonList;
class MaldesiIonListGenerator;
class MaldesiLibraryParamsSet;
class MaldesiParameters;

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
 * @brief MaldesiLibraryParamsSet
 * Scan-specific parameters should only be used if they are provided,
 * this structure retains some of the information implicit in how this data is encoded
 * in an mzkitcpp maldesi context.
 *
 * Specifically, if a compound has scan-specific information for at least one scan,
 * the compound should only be searched in scans where scan-specific information is provided.
 * this is recorded by adding the compound to the 'compoundsWithScanSpecificParams' set.
 *
 * If the compound is not provided with any scan-specific information, however, the compound is not
 * included in the 'compoundsWithScanSpecificParams' set, and downstream consumers should assume
 * that the compound should be searched in every scan, using the default parameters.
 *
 * In this way, the 'compoundsWithScanSpecificParams' serves multiple purposes: to indicate which
 * compounds have scan-specific params, and also to indicate which compounds should be searched in
 * a subset of scans or in every scan. This is necessary because absence from the 'compoundScanSpecificParamsMap'
 * is insufficient to indicate if a compound should be searched in every scan.
 */
class MaldesiLibraryParamsSet {
public:
    set<string> compoundsWithScanSpecificParams{};
    map<pair<string, int>, shared_ptr<MaldesiParameters>> compoundScanSpecificParamsMap{};
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

    //dual vectors, related by position, e.g. compoundNameVector[0] matches to encodedParamsVector[0]
    static MaldesiLibraryParamsSet decodeLibraryParamsSet(
        const vector<string>& compoundNameVector,
        const vector<string>& encodedParamsVector);

    //return encoded parameter value or defaults, if the encoded value corresponds to defaults
    const int getMinNumBoundLigand(const int defaultValue);
    const int getMaxNumBoundLigand(const int defaultValue);
    const double getBoundLigandExactMass(const double defaultValue);
    const vector<Adduct>& getAdducts(const vector<Adduct>& defaultValue);
};

#endif // MALDESI_H
