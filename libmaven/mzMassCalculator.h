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

//This enum only applies to labeled isotopes.
//isotopic species that contain labels from natural abundance may optionally be retained
//(handled downstream).
//
// Key case: isotope has label of interest along with some other labeled form, not of
// interest.
//
// Remember, other isotopes can still survive the check even if they are not labeled forms.
//
enum LabeledIsotopeRetentionPolicy {

    //Only one labeled form allowed per isotope.
    //The label must come from the list of valid labeled forms.
    ONLY_ONE_LABEL,

    //Single labeled species are allowed.
    //Double-labeled species are permitted if they include 13C as one of the labels.
    //All labeled forms must come from the list of valid labeled forms.
    //this is the choice for the original Maven 1.0 implementation.
    ONLY_CARBON_TWO_LABELS,

    //Any degree of labeling is supported, as long as
    //All labeled forms come from the list of valid labeled forms.
    ONE_OR_MORE_LABELS
};

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
    static double computeNeutralMass(string formula);
    static map<string,int> getComposition(string formula);
    static map<string,int> getComposition(Adduct* adduct);
    static void addAtoms(map<string, int>& reference, map<string, int> toAdd);
    static void subtractAtoms(map<string, int>& reference, map<string, int> toSubtract);
    static void multiplyAtoms(map<string, int>& reference, int factor);

    double computeMass(string formula, int polarity);
    double computeC13IsotopeMass(string formula);
    map<string,double>computeLabeledMasses(string formula, int polarity);
    map<string,double>computeLabeledAbundances(string formula);

    void matchMass(double mass, double ppm);
    string prettyName(int c, int h, int n, int o, int p, int s, int d=0);
    vector<Match> enumerateMasses(double inputMass, double charge, double maxdiff);
    double adjustMass(double mass,int charge);

    //Used for computation of actual values for labeled data.
    static vector<Isotope> computeIsotopes(string compoundFormula, Adduct* adduct, int maxNumProtons=INT_MAX, bool isUse13C=true, bool isUse15N=true, bool isUse34S=true, bool isUse2H=true);

    static vector<Isotope> computeIsotopes2(string compoundFormula,
                                            Adduct *adduct,
                                            vector<Atom> heavyIsotopes,
                                            LabeledIsotopeRetentionPolicy labeledIsotopeRetentionPolicy,
                                            NaturalAbundanceData naturalAbundanceData,
                                            bool isIncludeNaturalAbundance = false,
                                            int maxNumExtraNeutrons=INT_MAX,
                                            double minimumProportionMPlusZero = 0
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

    map<string,int> getPeptideComposition(const string& peptideSeq);

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

class Atom {
    public:
        string symbol;
        int massNumber;

        Atom() {symbol = "NONE"; massNumber = -1;}

    Atom(string symbol, int massNumber) {
        this->symbol = symbol;
        this->massNumber = massNumber;
    }

    friend bool operator< (const Atom& a, const Atom& b) {
        if (a.symbol == b.symbol) {
            return a.massNumber < b.massNumber;
        } else {
            return a.symbol < b.symbol;
        }
    }
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

        void setAtomData(string atomicSymbol, int massNumber, double atomicMass, double naturalAbundance, int numExtraNeutrons);
        void print();

        vector<Atom> getAtomsByExtraNeutrons(int numNeutrons);
        vector<Atom> getAtomsBySymbol(string atomicSymbol);

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
