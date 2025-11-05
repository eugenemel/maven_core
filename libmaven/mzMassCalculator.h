#ifndef MASSCALC_H
#define MASSCALC_H

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <map>
#include "mzSample.h" 
#include "mzUtils.h"
#include "Fragment.h"

class Compound;
class Isotope;
class Adduct;
class Atom;
class IsotopicAbundance;

using namespace std;

// EI    M+ M-
// Cl    M+H M+X
// ACPI  M+H M+X M-X
// FI    M+H  M+X
// ESI    [M+nH]^n+  [M-nX]^n-
// FAB   M+X M+N
// ACPI  M+H  M+X

class NaturalAbundanceData;
class NaturalAbundanceDistribution;

class MassCalculator { 

    public:
        enum IonizationType { ESI=0, EI=1};
        static IonizationType ionizationType;
        static Adduct* PlusHAdduct;
        static Adduct* MinusHAdduct;
        static Adduct* ZeroMassAdduct;

    struct Match {

        Match& operator= (const Match& b) {
            name = b.name;
            formula = b.formula;
            mass = b.mass;
            diff = b.diff;
            rtdiff = b.rtdiff;
            compoundLink = b.compoundLink;
            adductLink =   b.adductLink;
            fragScore = b.fragScore;
            atomCounts=b.atomCounts;
            return *this;
        }

        vector<int> compAtomCounts( Match& b) {
            vector<int>delta = this->atomCounts;
            if ( b.atomCounts.size() < this->atomCounts.size() ) return delta;
            for(int i=0; i< atomCounts.size(); i++ ) delta[i] = atomCounts[i] - b.atomCounts[i];
            return delta;
        }

        std::string name;
        std::string formula;
        double mass=0;
        double diff=0;
        double rtdiff=0;
        vector<int>atomCounts;
        Compound* compoundLink=NULL;
        Adduct*   adductLink=NULL;
        FragmentationMatchScore fragScore;
    };


    MassCalculator(){}
    static double computeNeutralMass(map<string, int>& atoms);
    static double computeNeutralMass(string formula);
    static double computePeptideNeutralMass(string peptideSequence);
    static map<string,int> getComposition(string formula);
    static map<string,int> getComposition(Adduct* adduct);
    static void addAtoms(map<string, int>& reference, map<string, int> toAdd);
    static void subtractAtoms(map<string, int>& reference, map<string, int> toSubtract);
    static void multiplyAtoms(map<string, int>& reference, int factor);

    static void applyMZeroMzOffset(vector<Isotope>& isotopes, double peakGroupMz);

    double computeMass(string formula, int polarity);
    double computeC13IsotopeMass(string formula);
    map<string,double>computeLabeledMasses(string formula, int polarity);
    map<string,double>computeLabeledAbundances(string formula);

    static int getAdductCharge(string adductName);
    static int getAdductMols(string adductName);
    static Adduct parseAdductFromName(string adductName);

    void matchMass(double mass, double ppm);
    string prettyName(int c, int h, int n, int o, int p, int s, int d=0);
    vector<Match> enumerateMasses(double inputMass, double charge, double maxdiff);
    double adjustMass(double mass,int charge);


    //Issue 758: predict molecular formula for an arbitrary m/z
    static vector<pair<Adduct, map<MassAtom, int>>> candidateAtomMaps(
        double mz,
        vector<Adduct> possibleAdducts,
        map<MassAtom, pair<int, int>> legalAtomCounts,
        double ppmDiff,
        long maxNumCandidates = 2e6,
        bool debug=false
        );

    static vector<pair<string, double>> evaluateAtomMapCandidates(
        double mz,
        const vector<pair<Adduct, map<MassAtom, int>>>& atomMapCandidates,
        bool debug = false);

    static vector<Isotope> computeIsotopes(
        string compoundFormula,
        Adduct *adduct,
        double mz,
        vector<Atom> labeledIsotopes,
        LabeledIsotopeRetentionPolicy labeledIsotopeRetentionPolicy,
        NaturalAbundanceData naturalAbundanceData,
        bool isIncludeNaturalAbundance = false,
        int maxNumExtraNeutrons=INT_MAX,
        double minimumProportionMPlusZero = 0,
        bool debug=false
        );

    //Issue 656: Return a complete natural abundance distribution for every possible observable isotope.
    //Considers all atoms in the molecule and adducts.
    //Organizes information in various ways for easy, intuitive access.
    static NaturalAbundanceDistribution getNaturalAbundanceDistribution(
        string compoundFormula,
        Adduct *adduct,
        NaturalAbundanceData& data,
        double minimumAbundance,
        bool debug=false);

    //Issue 684: Compute a set of stubbed IsotopicAbundance for a case where
    // the atomic composition of a species isn't known (unknown compound mz).
    static vector<IsotopicAbundance> getUnknownFormulaIsotopicAbundances(
        double mz,
        vector<Atom> labeledIsotopes,
        NaturalAbundanceData& naturalAbundanceData,
        int maxNumExtraNeutrons,
        bool debug=false
        );

    //Issue 715
    static string getCachedIsotopeKey(
        string formula,
        Adduct* adduct,
        vector<Atom> labeledIsotopes,
        LabeledIsotopeRetentionPolicy labeledIsotopeRetentionPolicy,
        bool isIncludeNaturalAbundance,
        int maxNumExtraNeutrons,
        double minimumProportionMPlusZero);

    //Issue 711: Declare cache
    static map<string, vector<Isotope>> isotopesCache;

    //Issue 813: make static
    static map<string,int> getPeptideComposition(const string& peptideSeq);

    static bool compDiff(const Match& a, const Match& b ) { return a.diff < b.diff; }
    static bool compNumMatches(const Match& a, const Match& b ) { return a.fragScore.numMatches < b.fragScore.numMatches; }
    static double getElementMass(string elmnt);
    static double getNaturalAbundanceCorrectedQuantValue(
        float isotopeObserved,
        float mZeroObserved,
        double isotopeExpectedAbundance,
        double mZeroExpectedAbundance);

    private:
        static void modifyAtoms(map<string, int>& reference, map<string, int> toAdd, bool isAddAtoms);
        static map<string, int> getAdductComponentComposition(string formula);

};

//flexible class to store natural abundance data. In some cases, these values can change - e.g., plant metabolomics.
class NaturalAbundanceData {
    public:
        //Keys:
        //   <atomicSymbol, massNumber>
        //
        // Values:
        //    atomToAbundance: Proportion of a particular atom with given # neutrons.
        //        Sum of all proportions of a particular atom sums to 1.0
        //
        //    atomToMass: Exact mass associated with a particular atom, # neutrons.
        //    Note that these values are always given in terms of the number of total neutrons,
        //    Not the number of "extra" neutrons to the most common case. This handles cases like Selennium,
        //    Where there isn't a single "base" case.
        map<Atom, double> atomToAbundance{};
        map<Atom, double> atomToMass{};
        map<Atom, int> atomToNumExtraNeutrons{};
        map<int, vector<Atom>> extraNeutronToAtoms{};
        map<string, vector<Atom>> symbolToAtoms{};

        //In the case where there are only two possible isotopes for a given atom,
        // retrieve a deltaMz based on the atom string.
        map<string, double> atomToDeltaMz{};

        void setAtomData(string atomicSymbol, int massNumber, double atomicMass, double naturalAbundance, int numExtraNeutrons);
        void print();

        vector<Atom> getAtomsByExtraNeutrons(int numNeutrons);
        vector<Atom> getAtomsBySymbol(string atomicSymbol);
        double getDeltaMzBySymbol(string atomicSymbol);

        static NaturalAbundanceData defaultNaturalAbundanceData;
};

class IsotopicAbundance {
    public:

        map<Atom, int> atomCounts{};

        //represented as an absolute theoretical probability.
        //The sum of all natural abundance values for all possible atom combinations is exactly one.
        double naturalAbundance = 1.0;

        //represented as a proportion of the natural abundance of the [M+0] species.
        //The [M+0] is exactly 1. Computed by this->naturalAbundance/[M+0]->naturalAbundance
        double naturalAbundanceMonoProportion = 1.0;

        //fields computed via IsotopicAbundance::compute()
        double mass = 0.0;
        double mz = 0.0;
        unsigned int numTotalExtraNeutrons = 0;
        set<Atom> labeledForms{};
        set<Atom> unlabeledForms{};

        void compute(NaturalAbundanceData& naturalAbundanceData, unsigned int chgNumber);

        bool isHasAtom(Atom& atom);

        static IsotopicAbundance createMergedAbundance(IsotopicAbundance& one, IsotopicAbundance& two);
        string getFormula();
        string toString();

        Isotope toIsotope();
};

//Issue 656: Implement flexible approach for isotopic correction.
//Takes into account all atoms, all
class NaturalAbundanceDistribution {
public:
    vector<IsotopicAbundance> isotopicAbundances{};

    // <naturalAbundance, naturalAbundanceMonoProportion>
    pair<double, double> getIsotopicAbundance(Isotope& isotope);
};

#endif
