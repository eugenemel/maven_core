#include "maldesi.h"
#include "mzMassCalculator.h"

MaldesiIonList MaldesiIonListGenerator::getMaldesiIonList(
    double compoundMz,
    vector<Adduct>& adducts,
    bool ms1UseDaTol,
    double ms1PpmTol,
    double ms1DaTol,
    bool debug) {

    MaldesiIonList ionList;

    ionList.isSinglePrecursorMz = adducts.empty();

    if (adducts.empty()) {
        // If no adducts supplied, assume that a single precursor m/z is provided (no need for manipulation)
        ionList.searchMz.push_back(compoundMz);
        ionList.ionName.push_back("");
        ionList.chg.push_back(1);

    } else {
        for (Adduct& adduct : adducts) {
            double precursorMz = adduct.computeAdductMass(compoundMz);
            ionList.searchMz.push_back(precursorMz);
            ionList.ionName.push_back(adduct.name);
            ionList.chg.push_back(adduct.charge);

            if (debug) {
                cout << "MaldesiIonListGenerator::getMaldesiIonList(): " << compoundMz << " " << adduct.name << ": mz=" << precursorMz << endl;
            }
        }
    }

    for (double mz : ionList.searchMz) {
        double minMz, maxMz = 0;
        if (ms1UseDaTol) {
            minMz = mz - ms1DaTol;
            maxMz = mz + ms1DaTol;
        } else {
            minMz = mz - mz*ms1PpmTol/1e6;
            maxMz = mz + mz*ms1PpmTol/1e6;
        }
        ionList.minMz.push_back(minMz);
        ionList.maxMz.push_back(maxMz);
    }

    return ionList;
}

MaldesiIonList MaldesiIonListGenerator::getLargePeptideProteinBindingAssayIonList(
    string peptideSequence,
    vector<Adduct>& adducts,
    double boundLigandExactMass,
    int maxNumBoundLigand,
    double peptidePredictedIsotopeRatioThreshold,
    bool ms1UseDaTol,
    double ms1PpmTol,
    double ms1DaTol,
    bool debug
    ) {
    MaldesiIonList ionList;
    ionList.isSinglePrecursorMz = false;

    double peptideExactMass = MassCalculator::computePeptideNeutralMass(peptideSequence);

    map<int, double> peptideIsotopeDist = MassCalculator::peptideC13Distribution(
        peptideSequence,
        0.0107, // natural abundance of C13
        peptidePredictedIsotopeRatioThreshold,
        false); // debug

    int numBoundLigand = 0;
    while (numBoundLigand <= maxNumBoundLigand) {
        for (auto it = peptideIsotopeDist.begin(); it != peptideIsotopeDist.end(); ++it) {
            for (Adduct& adduct : adducts) {

                stringstream s;
                s << std::fixed << setprecision(4);

                int numC13 = it->first;
                double isotopeProbability = it->second;
                string adductName = adduct.name;

                s << adductName << " "
                  << "[M+" << numC13 << "] "
                  << "f=" << isotopeProbability << " "
                  << "L=" << numBoundLigand;

                string ionName = s.str();

                // peptide exact mass + ligand mass shift + C13 mass shift
                double combinedExactMass = peptideExactMass + numBoundLigand * boundLigandExactMass + 1.003354835 * numC13;

                double mz = adduct.computeAdductMass(combinedExactMass);

                if (debug) {
                    cout << "getLargePeptideProteinBindingAssayIonList(): "<< peptideSequence << " " << ionName << ": mz=" << mz << endl;
                }

                double minMz, maxMz = 0;
                if (ms1UseDaTol) {
                    minMz = mz - ms1DaTol;
                    maxMz = mz + ms1DaTol;
                } else {
                    minMz = mz - mz*ms1PpmTol/1e6;
                    maxMz = mz + mz*ms1PpmTol/1e6;
                }

                // record results in ion list
                ionList.searchMz.push_back(mz);
                ionList.minMz.push_back(minMz);
                ionList.maxMz.push_back(maxMz);
                ionList.ionName.push_back(ionName);
                ionList.chg.push_back(adduct.charge);
            }
        }

        numBoundLigand++;
    }

    return ionList;
}
