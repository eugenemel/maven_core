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


using namespace std;

// EI    M+ M-
// Cl    M+H M+X
// ACPI  M+H M+X M-X
// FI    M+H  M+X
// ESI    [M+nH]^n+  [M-nX]^n-
// FAB   M+X M+N
// ACPI  M+H  M+X



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

    static vector<Isotope> computeIsotopes(string compoundFormula, Adduct* adduct, int maxNumProtons=INT_MAX, bool isUse13C=true, bool isUse15N=true, bool isUse34S=true, bool isUse2H=true);
    map<string,int> getPeptideComposition(const string& peptideSeq);

    static bool compDiff(const Match& a, const Match& b ) { return a.diff < b.diff; }
    static bool compNumMatches(const Match& a, const Match& b ) { return a.fragScore.numMatches < b.fragScore.numMatches; }

    private:
        static double getElementMass(string elmnt);
        static void modifyAtoms(map<string, int>& reference, map<string, int> toAdd, bool isAddAtoms);
        static map<string, int> getAdductComponentComposition(string formula);

};

#endif
