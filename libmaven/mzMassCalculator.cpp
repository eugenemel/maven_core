#include "mzMassCalculator.h"
#include <regex>

//Issue 711: cache must be initialized out-of-line
map<string, vector<Isotope>> MassCalculator::isotopesCache = {};

using namespace mzUtils;
using namespace std;

/*
aaLookup["G"] = 57.021464;
aaLookup["G"] = 57.021464;
aaLookup["A"] = 71.037114;
aaLookup["S"] = 87.032029;
aaLookup["P"] = 97.052764;
aaLookup["V"] = 99.068414;
aaLookup["T"] = 101.04768;
aaLookup["C"] = 103.00919;
aaLookup["L"] = 113.08406;
aaLookup["I"] = 113.08406;
aaLookup["N"] = 114.04293;
aaLookup["D"] = 115.02694;
aaLookup["Q"] = 128.05858;
aaLookup["K"] = 128.09496;
aaLookup["E"] = 129.04259;
aaLookup["M"] = 131.04048;
aaLookup["H"] = 137.05891;
aaLookup["F"] = 147.06841;
aaLookup["R"] = 156.10111;
aaLookup["Y"] = 163.06333;
aaLookup["W"] = 186.07931;
*/

MassCalculator::IonizationType MassCalculator::ionizationType = MassCalculator::ESI;
Adduct* MassCalculator::PlusHAdduct  = new Adduct("[M+H]+",  PROTON , 1, 1);
Adduct* MassCalculator::MinusHAdduct = new Adduct("[M-H]-", -PROTON, -1, 1);
Adduct* MassCalculator::ZeroMassAdduct = new Adduct("[M]",0 ,1, 1);

string MassCalculator::getCachedIsotopeKey(
    string formula,
    Adduct* adduct,
    vector<Atom> labeledIsotopes,
    LabeledIsotopeRetentionPolicy labeledIsotopeRetentionPolicy,
    bool isIncludeNaturalAbundance,
    int maxNumExtraNeutrons,
    double minimumProportionMPlusZero
    ) {

    string key = "formula=" + formula + ";";
    if (adduct) {
        key = key + "adduct=" + adduct->name + ";";
    }

    key = key + "{";
    for (auto& atom : labeledIsotopes) {
        key = key + to_string(atom.massNumber) + atom.symbol + ";";
    }
    key = key + "}";

    key = key + "labeledIsotopeRetentionPolicy=" + to_string(labeledIsotopeRetentionPolicy) + ";";
    key = key + "isIncludeNaturalAbundance=" + to_string(isIncludeNaturalAbundance) + ";";
    key = key + "maxNumExtraNeutrons=" + to_string(maxNumExtraNeutrons) + ";";
    key = key + "minimumProportionMPlusZero=" + to_string(minimumProportionMPlusZero) + ";";

    return key;
}


/*---------- function to get molar weight of an element or group --------*/
double MassCalculator::getElementMass(string elmnt){

    /* default behavior is to ignore string */
    double val_atome(0);

/* Check for atoms */
    if (elmnt == "H")        val_atome = 1.0078250321;
    else if (elmnt == "2H")   val_atome = 2.01410178;
    else if (elmnt == "D")   val_atome = 2.01410178;
    else if (elmnt == "C")   val_atome = 12.00000000;
    else if (elmnt == "13C") val_atome = 13.003354835336;
    else if (elmnt == "N")   val_atome = 14.0030740052;
    else if (elmnt == "15N") val_atome = 15.000108898266;
    else if (elmnt == "O")   val_atome = 15.994914619257;
    else if (elmnt == "18O") val_atome = 17.999159612136;
    else if (elmnt == "S")   val_atome = 31.97207069;
    else if (elmnt == "34S") val_atome = 33.96786701;
    else if (elmnt == "P")   val_atome = 30.97376151;
    else if (elmnt == "F")   val_atome = 18.9984032;
    else if (elmnt == "Na")  val_atome = 22.98976967;
    else if (elmnt == "Mg")  val_atome = 24.98583702;
    else if (elmnt == "Ca")  val_atome = 39.962591;
    else if (elmnt == "Cl")  val_atome = 34.96885271;
    else if (elmnt == "Br")  val_atome = 78.918336;
    else if (elmnt == "I")   val_atome = 126.904477;
    else if (elmnt == "K")   val_atome = 38.96399867;
    else if (elmnt == "Ca")  val_atome = 39.9625912;
    else if (elmnt == "Se")  val_atome = 79.916521;
    else if (elmnt == "As")  val_atome = 74.921596;
    else if (elmnt == "Si")  val_atome = 27.9769265325;
    else if (elmnt == "Fe")  val_atome = 55.934939;
    else if (elmnt == "Li")  val_atome = 7.016003437;
    return(val_atome);
}
/*-----------------------------------------------------------------------*/

map<string,int> MassCalculator::getComposition(Adduct* adduct){

    //cout << "MassCalculator::getComposition(Adduct*)" << endl;

    map<string, int> atoms {};

    if (!adduct) return atoms;

    string name = adduct->name;

    bool isAfterM = false;

    int formulaStart = -1;
    int formulaEnd = -1;

    vector<string> formulasToAdd;
    vector<string> formulasToSubtract;

    bool isAddFormula = false;

    //use string tokenizer approach to build up map
    for (int i = 0; i < name.length(); i++) {

        if (isAfterM) {

            if (name[i] == ']') {

                //write previous entry (if appropriate)
                if (formulaStart > 0) {

                    if (isAddFormula) {
                        formulasToAdd.push_back(name.substr(formulaStart, (formulaEnd-formulaStart+1)));
                    } else {
                        formulasToSubtract.push_back(name.substr(formulaStart, (formulaEnd-formulaStart+1)));
                    }

                    formulaStart = -1;
                    formulaEnd = -1;
                }

                break;

            } else if (name[i] == '-') {

                //write previous entry (if appropriate)
                if (formulaStart > 0) {

                    if (isAddFormula) {
                        formulasToAdd.push_back(name.substr(formulaStart, (formulaEnd-formulaStart+1)));
                    } else {
                        formulasToSubtract.push_back(name.substr(formulaStart, (formulaEnd-formulaStart+1)));
                    }

                    formulaStart = -1;
                    formulaEnd = -1;
                }

                isAddFormula = false;

            } else if (name[i] == '+') {

                //write previous entry (if appropriate)
                if (formulaStart > 0) {

                    if (isAddFormula) {
                        formulasToAdd.push_back(name.substr(formulaStart, (formulaEnd-formulaStart+1)));
                    } else {
                        formulasToSubtract.push_back(name.substr(formulaStart, (formulaEnd-formulaStart+1)));
                    }

                    formulaStart = -1;
                    formulaEnd = -1;
                }

                isAddFormula = true;

            } else { // in the middle of a name

                if (formulaStart == -1) formulaStart = i;
                formulaEnd = i;
            }

        } else if (name[i] == 'M') {
            isAfterM = true;
        }
    }

    for (string posFormula : formulasToAdd) {
        addAtoms(atoms, getAdductComponentComposition(posFormula));
    }

    for (string negFormula : formulasToSubtract) {
        subtractAtoms(atoms, getAdductComponentComposition(negFormula));
    }

    return atoms;
}

int MassCalculator::getAdductCharge(string adductName) {
    int adductCharge = 1;
    regex suffixRe("\\](.*)$");
    regex chgNumRe("(\\d+)[+-]");

    smatch match, match2;
    if (regex_search(adductName, match, suffixRe)) {
        string suffix = match[1].str();
        if (regex_search(suffix, match2, chgNumRe)) {
            string suffixStr = match2[1].str();
            adductCharge = stoi(suffixStr.substr(0, suffixStr.size()-2));
            if (suffix.substr(suffix.size()-1, 1) == "-") {
                adductCharge = -1 * adductCharge;
            }
        } else if (suffix == "-") {
            adductCharge = -1;
        }
    }

    return adductCharge;
}

int MassCalculator::getAdductMols(string adductName) {
    int adductNumMols = 1;
    regex grabNumMols("\\[([^+-]+)[+-]");
    regex grabInts("^([^M]*)");
    smatch match, match2;
    if (regex_search(adductName, match, grabNumMols)) {
        string numMolsStr = match[1].str();
        if (regex_search (numMolsStr, match2, grabInts)) {
            string numInts = match2[1].str();
            if (numInts != "") {
                adductNumMols = stoi(numInts);
            }
        }
    }

    return adductNumMols;
}

Adduct MassCalculator::parseAdductFromName(string adductName) {

    Adduct adduct;
    adduct.name = adductName;

    map<string, int> atomComposition = MassCalculator::getComposition(&(adduct));
    int adductCharge = MassCalculator::getAdductCharge(adductName);
    int adductMols = MassCalculator::getAdductMols(adductName);

    //# of electrons gained or lost
    double massShift = EMASS * -1 * adductCharge;

    for (auto it = atomComposition.begin(); it != atomComposition.end(); ++it) {
        massShift += MassCalculator::getElementMass(it->first) * it->second;
    }
    adduct.mass = massShift;
    adduct.nmol = adductMols;
    adduct.charge = adductCharge;

    return adduct;
}

void MassCalculator::addAtoms(map<string, int>& reference, map<string, int> toAdd) {
    modifyAtoms(reference, toAdd, true);
}

void MassCalculator::subtractAtoms(map<string, int>& reference, map<string, int> toSubtract) {
    modifyAtoms(reference, toSubtract, false);
}

void MassCalculator::multiplyAtoms(map<string, int>& reference, int factor) {
    for (map<string, int>::iterator it = reference.begin(); it != reference.end(); ++it){
        reference[it->first] *= factor;
    }
}

void MassCalculator::modifyAtoms(map<string, int>& reference, map<string, int> toAdjust, bool isAddAtoms) {

    for (map<string, int>::iterator it = toAdjust.begin(); it != toAdjust.end(); ++it) {

        string element = it->first;
        int num = it->second;

        if (reference.find(element) == reference.end()) {

            if (!isAddAtoms) {
                num *= -1;
            }

            reference.insert(make_pair(element, num));

        } else {
            if (isAddAtoms) {
                reference[element] += num;
            } else {
                reference[element] -= num;
            }
        }
    }

}

map<string, int> MassCalculator::getAdductComponentComposition(string formula){
    if (formula == "ACN") { // acetonitrile
        return map<string, int>{{"C", 2}, {"H", 3}, {"N", 1}};
    } else if (formula == "2ACN") { //2 * acetonitrile
        return map<string, int>{{"C", 4}, {"H", 6}, {"N", 2}};
    } else if (formula == "3ACN") { // 3 * acetonitrile
        return map<string, int>{{"C", 6}, {"H", 9}, {"N", 3}};
    } else if (formula == "TFA") { //trifluoroacetate
        return map<string, int>{{"C", 2}, {"H", 1}, {"O", 2}, {"F", 3}};
    } else if (formula == "FA") { //formic acid (HCOOH)
        return map<string, int>{{"C", 1}, {"H", 2}, {"O", 2}};
    } else if (formula == "tributylamine") {
        return map<string,int>{{"C", 12}, {"H", 27}, {"N", 1}};
    } else if (formula == "3H2O") {
        return map<string, int>{{"H", 6}, {"O", 3}};
    } else if (formula == "2H2O") {
        return map<string, int>{{"H", 4}, {"O", 2}};
    } else if (formula == "2H") {
        return map<string, int>{{"H", 2}};
    } else if (formula == "3H") {
        return map<string, int>{{"H", 3}};
    } else if (formula == "AcOH") { //protonated acetate (C2H4O2)
        return map<string, int>{{"C", 2}, {"H", 4}, {"O", 2}};
    } else if (formula == "AcO") { //deprotonated acetate (C2H3O2)
        return map<string, int>{{"C", 2}, {"H", 3}, {"O", 2}};
    } else if (formula == "NaAcOH") { //Na + protonated acetate (C2H4O2)
        return map<string, int>{{"C", 2}, {"H", 4}, {"Na", 1}, {"O", 2}};
    } else if (formula == "NaAcO") { //Na + deprotonated acetate (C2H3O2)
        return map<string, int>{{"C", 2}, {"H", 3}, {"Na", 1}, {"O", 2}};
    } else if (formula == "2Na") {
        return map<string, int>{{"Na", 2}};
    } else if (formula == "3Na") {
        return map<string, int>{{"Na", 3}};
    } else if (formula == "2K") {
        return map<string, int>{{"K", 2}};
    } else if (formula == "3K") {
        return map<string, int>{{"K", 3}};
    } else if (formula == "DMSO") { //Dimethyl sulfoxide (C2H6OS)
        return map<string, int>{{"C", 2}, {"H", 6}, {"O", 1}, {"S", 1}};
    } else if (formula == "IsoProp") { //Isopropylene (C3H8O)
        return map<string, int>{{"C", 3}, {"H", 8}, {"O", 1}};
    } else if (formula == "CH3OH") {
        return map<string, int>{{"C", 1}, {"H", 4}, {"O", 1}};
    } else if (formula == "Na2O2CH") {
        return map<string, int>{{"C", 1}, {"H", 1}, {"Na", 2}, {"O", 2}};
    } else if (formula == "NaHO2CH") {
        return map<string, int>{{"C", 1}, {"H", 2}, {"Na", 1}, {"O", 2}};
    } else if (formula == "NaO2CH") {
        return map<string, int>{{"C", 1}, {"H", 1}, {"Na", 1}, {"O", 2}};
    } else if (formula == "Hac") { //protonated acetate
        return map<string, int>{{"C", 2}, {"H", 4}, {"O", 2}};
    } else if (formula == "37Cl") {
        return map<string, int>{{"Cl", 1}}; //Issue 378: Note that non-monoisotopic aspect is lost here
    } else {
        return getComposition(formula);
    }
}


/*-------------- parsing function ---------------------------------------*/

map<string,int> MassCalculator::getComposition(string formula) {

	/* define allowed characters for formula */
	const string UPP("ABCDEFGHIKLMNOPRSTUVWXYZ");
	const string LOW("abcdefghiklmnoprstuy");
	const string COEFF("0123456789");

	/* define some variable */
	int SIZE = formula.length();
	map<string,int> atoms;

	/* parse the formula */
	for(int i = 0; i < SIZE; i++) {
		string bloc, coeff_txt;
		int coeff;

		/* start of symbol must be uppercase letter */
		if (UPP.find(formula[i]) != string::npos){
			bloc = formula[i];
			if (LOW.find(formula[i+1]) != string::npos){
				bloc += formula[i+1];
				i++;
			}
		}

		while ( COEFF.find(formula[i+1]) != string::npos ) {
			coeff_txt += formula[i+1];
			i++;
		}

		if ( coeff_txt.length() > 0 ) {
				coeff = string2integer(coeff_txt);
		} else {
				coeff = 1;
        }

		/* compute normally if there was no open bracket */
		//cout << bloc <<  " " << coeff << endl;
		atoms[bloc] += coeff;
	}

	return(atoms);
	/* send back value to main.cpp */
}

double MassCalculator::computeNeutralMass(string formula) {
	map<string,int> atoms = getComposition(formula);
	map<string,int>::iterator itr;

	double mass=0;
    for(itr = atoms.begin(); itr != atoms.end(); itr++ ) {
        mass += getElementMass((*itr).first) * (*itr).second;
	}
	return mass;
}


double MassCalculator::adjustMass(double mass,int charge) {

    if (MassCalculator::ionizationType == EI and charge !=0) {
       return ((mass - charge*EMASS)/charge); // loss of electrons
	}

    if (charge == 0 ) return mass;
    else return   (mass+charge*PROTON)/abs(charge);
}


//input neutral formala with all the hydrogens and charge state of molecule
//output expected mass of the molecule after loss/gain of protons
double MassCalculator::computeMass(string formula, int charge) {
    double mass = computeNeutralMass(formula);
    return adjustMass(mass,charge);
}

vector<Isotope> MassCalculator::computeIsotopes(
    string compoundFormula,
    Adduct *adduct,
    double mz,
    vector<Atom> labeledIsotopes,
    LabeledIsotopeRetentionPolicy labeledIsotopeRetentionPolicy,
    NaturalAbundanceData naturalAbundanceData,
    bool isIncludeNaturalAbundance,
    int maxNumExtraNeutrons,
    double minimumProportionMPlusZero,
    bool debug
    ){

    vector<IsotopicAbundance> isotopicAbundances;
    string cacheKey = "";

    bool isKnownAtomicComposition = compoundFormula != "" && adduct;

    // known atomic composition case
    if (isKnownAtomicComposition) {

       if (debug) cout << "MassCalculator::computeIsotopes(): Known atomic composition case." << endl;

       cacheKey = MassCalculator::getCachedIsotopeKey(
           compoundFormula,
           adduct,
           labeledIsotopes,
           labeledIsotopeRetentionPolicy,
           isIncludeNaturalAbundance,
           maxNumExtraNeutrons,
           minimumProportionMPlusZero);

       if (isotopesCache.find(cacheKey) != isotopesCache.end()) {
                if (debug) cout << "Found Cached vector<Isotope> - returning cached value!" << endl;
                return (isotopesCache.at(cacheKey));
       } else if (debug) {
                cout << "vector<Isotope> not found in cache, will re-compute." << endl;
       }

       //First, enumerate all theoretically possible isotopes based on the formula.
       NaturalAbundanceDistribution abundanceDistribution =
           MassCalculator::getNaturalAbundanceDistribution(
               compoundFormula,
               adduct,
               naturalAbundanceData,

               0, //Issue 695: reversion of Issue 690 - retain 0-intensity isotopes

               // Issue 690: speedup
               // Assume that the monoisotope constitutes at least 1% of total abundance
               // 0.01 * minimumProportionMPlusZero,

               false);

       isotopicAbundances = abundanceDistribution.isotopicAbundances;

    // unknown atomic composition case
    } else {

       if (debug) cout << "MassCalculator::computeIsotopes(): Unknown atomic composition case." << endl;

       isotopicAbundances = MassCalculator::getUnknownFormulaIsotopicAbundances(
           mz,
           labeledIsotopes,
           naturalAbundanceData,
           maxNumExtraNeutrons,
           debug
           );
    }

    //Filter the list of isotopes based on search criteria.
    vector<Isotope> isotopes{};

    for (auto isotopicAbundance : isotopicAbundances) {

       //Avoid isotopes with too many extra neutrons
       if (isotopicAbundance.numTotalExtraNeutrons > maxNumExtraNeutrons) continue;

       if (debug) {
            cout << isotopicAbundance.toString()
                 << ": " << isotopicAbundance.numTotalExtraNeutrons
                 << " extra neutrons." << endl;
       }

       //check that the isotope is one of the preferred label types
       bool isLabelType = false;

       set<Atom> labeledIsotopesPresent{};

       for (auto& atom : labeledIsotopes) {
            bool isHasAtom = isotopicAbundance.isHasAtom(atom);
            if (isHasAtom) {
                labeledIsotopesPresent.insert(atom);
            }
       }

       bool isSinglePrespecifiedLabel = labeledIsotopesPresent.size() == 1 && isotopicAbundance.labeledForms.size() == 1;

       bool isCarbon13MultipleLabelCase = labeledIsotopesPresent.size() == 2 && isotopicAbundance.labeledForms.size() == 2 &&
                                        labeledIsotopesPresent.find(Atom("C", 13)) != labeledIsotopesPresent.end();

       bool isMultipleLabeledForms = !labeledIsotopesPresent.empty() &&
                                     labeledIsotopesPresent.size() == isotopicAbundance.labeledForms.size();

       //Assess agreement based on isotope label retention policy
       if (labeledIsotopeRetentionPolicy == LabeledIsotopeRetentionPolicy::ONLY_ONE_LABEL) {
            isLabelType = isSinglePrespecifiedLabel;
       } else if (labeledIsotopeRetentionPolicy == LabeledIsotopeRetentionPolicy::ONLY_CARBON_TWO_LABELS) {
            isLabelType = isSinglePrespecifiedLabel || isCarbon13MultipleLabelCase;
       } else if (labeledIsotopeRetentionPolicy == LabeledIsotopeRetentionPolicy::ONE_OR_MORE_LABELS) {
            isLabelType = isMultipleLabeledForms;
       }

       //isotopes can also be included if they have a high enough natural abundance.
       bool isValidNaturalAbundance = isIncludeNaturalAbundance && isotopicAbundance.naturalAbundanceMonoProportion > minimumProportionMPlusZero;

       if (debug) {
            cout << "isLabelType? " << isLabelType
                 << ", isValidNaturalAbundance? " << isValidNaturalAbundance
                 <<", is [M+0]? " << (isotopicAbundance.numTotalExtraNeutrons == 0)
                 << ", is unknown atomic composition? " << !isKnownAtomicComposition
                 << endl;
       }

       //retain isotopes if they are the [M+0], of the preferred label type, pass natural abundance criteria, or are associated with an unknown formula.
       if (isLabelType || isValidNaturalAbundance || isotopicAbundance.numTotalExtraNeutrons == 0 || !isKnownAtomicComposition) {
           isotopes.push_back(isotopicAbundance.toIsotope());
       }
    }

    //Issue 711: Write value to cache to avoid excessive recomputations.
    if (cacheKey != "") {
        isotopesCache.insert(make_pair(cacheKey, isotopes));
    }

    if (debug) {
        cout << "MassCalculator::computeIsotopes() Created " << isotopes.size() << " isotopes." << endl;
    }

    return isotopes;
}

vector<IsotopicAbundance> MassCalculator::getUnknownFormulaIsotopicAbundances(
    double mz,
    vector<Atom> heavyIsotopes,
    NaturalAbundanceData& naturalAbundanceData,
    int maxNumExtraNeutrons,
    bool debug
    ) {

    //initialize abundances
    vector<IsotopicAbundance> isotopicAbundances{};

    IsotopicAbundance seed;
    seed.mz = mz;
    isotopicAbundances.push_back(seed);

    //Enumerate every possible combination of heavy atoms, store into atomCounts
    for (Atom& atom : heavyIsotopes) {

        vector<IsotopicAbundance> isotopicAbundancesAddCurrentAtom{};

        for (IsotopicAbundance abundance : isotopicAbundances) {
           for (unsigned int i = 0; i <= maxNumExtraNeutrons; i++) {

                //make a copy
                IsotopicAbundance isotopicAbundance = abundance;
                isotopicAbundance.atomCounts.insert(make_pair(atom, i));
                isotopicAbundance.numTotalExtraNeutrons++;

                //adjust mz
                isotopicAbundance.mz += naturalAbundanceData.getDeltaMzBySymbol(atom.symbol) * i;

                if (isotopicAbundance.numTotalExtraNeutrons <= maxNumExtraNeutrons) {
                    isotopicAbundancesAddCurrentAtom.push_back(isotopicAbundance);
                }
           }
        }

        isotopicAbundances = isotopicAbundancesAddCurrentAtom;
    }

    return isotopicAbundances;
}
/**
 * @brief MassCalculator::enumerateMasses
 * @param inputMass
 * @param charge
 * @param maxdiff
 * @return
 *
 * @deprecated This does not appear to have any active uses in maven_core or maven.
 */
vector<MassCalculator::Match> MassCalculator::enumerateMasses(double inputMass, double charge, double maxdiff) {

    vector<MassCalculator::Match>matches;
    if (charge > 0 ) inputMass = inputMass*abs(charge)-HMASS*abs(charge);
    if (charge < 0 ) inputMass = inputMass*abs(charge)+HMASS*abs(charge);

    for(int c=0; c<30; c++) { //C
        if (c*12 > inputMass) break;
        for(int n=0; n<30; n++) {//N
            if (c*12+n*14 > inputMass) break;
            for(int o=0; o<30;o++){//O
                if (c*12+n*14+o*16 > inputMass) break;
                for(int p=0; p<6;p++) { //P
                    for(int s=0; s<6;s++) { //S
                        int hmax = c*4+o*2+n*4+p*3+s*3;
                        for(int h=0; h<hmax;h++) { //H
                            //double du = ((c*2+n+p)+2-h)/2;
                            //if (du < -0.5 ) continue;
                            //if (round(du / 0.5) != (du/0.5) ) continue;

                            double c12 = c*12.0 +  o*15.9949146221 + n*14.0030740052 + p*30.97376151 + h*1.0078250321 + s*31.97207069;
                            double diff = ppmDist(c12,inputMass);

                            if ( diff < maxdiff) {
                                string name = prettyName(c,h,n,o,p,s);
                                MassCalculator::Match m;
                                m.mass=c12;
                                m.diff= diff;
                                m.compoundLink = NULL;
                                m.adductLink=NULL;
                                m.formula=m.name=name;
                                int tmp[6]= { c,h,n,o,p,s };
                                m.atomCounts.assign(tmp,tmp+6);
                                matches.push_back(m);
                            }
                        }
                    }
                }
            }
        }
    }
    std::sort(matches.begin(), matches.end(), compDiff );
    return matches;
}

std::string MassCalculator::prettyName(int c, int h, int n, int o, int p, int s, int d) {
		char buf[1000];
		string name;
		if ( c != 0 ) { name += "C"; if (c>1) { sprintf(buf,"%d",c);  name += buf;} }
		if ( h != 0 ) { name += "H"; if (h>1) { sprintf(buf,"%d",h);  name += buf;} }
        if ( d != 0 ) { name += "D"; if (d>1) { sprintf(buf,"%d",d);  name += buf;} }
		if ( n != 0 ) { name += "N"; if (n>1) { sprintf(buf,"%d",n);  name += buf;} }
		if ( o != 0 ) { name += "O"; if (o>1) { sprintf(buf,"%d",o);  name += buf;} }
		if ( p != 0 ) { name += "P"; if (p>1) { sprintf(buf,"%d",p);  name += buf;} }
		if ( s != 0 ) { name += "S"; if (s>1) { sprintf(buf,"%d",s);  name += buf;} }
		return name;
}

map<string,int> MassCalculator::getPeptideComposition(const string& peptideSeq) {
    int C=0; int H=0; int O=0; int N=0; int S=0; int P=0;

    for(int i=0; i<peptideSeq.length(); i++ ) {
        char aa = peptideSeq[i];
        switch (aa)  {
            case 'A': C+=3; H+=5;  N+=1;  O+=1; break;
            case 'R': C+=6; H+=12; N+=4;  O+=1; break;
            case 'N': C+=4; H+=6;  N+=2;  O+=2; break;
            case 'D': C+=4; H+=5;  N+=1;  O+=3; break;
            case 'C': C+=3; H+=5;  N+=1;  O+=1;  S+=1; break;
            case 'Q': C+=5; H+=8;  N+=2;  O+=2; break;
            case 'E': C+=5; H+=7;  N+=1;  O+=3; break;
            case 'G': C+=2; H+=3;  N+=1;  O+=1; break;
            case 'H': C+=6; H+=7;  N+=3;  O+=1; break;
            case 'I': C+=6; H+=11; N+=1;  O+=1; break;
            case 'L': C+=6; H+=11; N+=1;  O+=1; break;
            case 'K': C+=6; H+=12; N+=2;  O+=1; break;
            case 'M': C+=5; H+=9;  N+=1;  O+=1; S+=1 ;break;
            case 'F': C+=9; H+=9;  N+=1;  O+=1; break;
            case 'P': C+=5; H+=7;  N+=1;  O+=1; break;
            case 'S': C+=3; H+=5;  N+=1;  O+=2; break;
            case 'T': C+=4; H+=7;  N+=1;  O+=2; break;
            case 'W': C+=11; H+=10;  N+=2;  O+=1; break;
            case 'Y': C+=9; H+=9;  N+=1;  O+=2; break;
            case 'V': C+=5; H+=9;  N+=1;  O+=1; break;
        }
    }

    map<string,int> M;
    M["C"]=C;
    M["H"]=H;
    M["N"]=N;
    M["O"]=O;
    M["S"]=S;
    M["P"]=P;

    return M;
}
/*-----------------------------------------------------------------------*/

void NaturalAbundanceData::setAtomData(
    string atomicSymbol,
    int massNumber,
    double atomicMass,
    double naturalAbundance,
    int numExtraNeutrons){

    auto key = Atom(atomicSymbol, massNumber);

    atomToMass.insert(make_pair(key, atomicMass));
    atomToAbundance.insert(make_pair(key, naturalAbundance));
    atomToNumExtraNeutrons.insert(make_pair(key, numExtraNeutrons));

    if (extraNeutronToAtoms.find(numExtraNeutrons) == extraNeutronToAtoms.end()) {
        extraNeutronToAtoms.insert(make_pair(numExtraNeutrons, vector<Atom>{}));
    }
    extraNeutronToAtoms.at(numExtraNeutrons).push_back(key);

    if (symbolToAtoms.find(atomicSymbol) == symbolToAtoms.end()) {
        symbolToAtoms.insert(make_pair(atomicSymbol, vector<Atom>{}));
    }
    symbolToAtoms.at(atomicSymbol).push_back(key);

}


//Issue 656: Initialize default Natural Abundance via C++11 idiom for initialization of static field.
NaturalAbundanceData NaturalAbundanceData::defaultNaturalAbundanceData = []() -> NaturalAbundanceData {
    NaturalAbundanceData abundanceData;

    abundanceData.setAtomData("C", 12, 12.0, 0.9893, 0);
    abundanceData.setAtomData("C", 13, 13.00335483521, 0.0107, 1);

    abundanceData.setAtomData("H", 1, 1.00782503224, 0.999885, 0);
    abundanceData.setAtomData("H", 2, 2.01410177811, 0.000115, 1);

    abundanceData.setAtomData("O", 16, 15.994914619257, 0.99757, 0);
    abundanceData.setAtomData("O", 18, 17.999159612136, 0.00205, 2);

    abundanceData.setAtomData("N", 14, 14.00307400446, 0.99632, 0);
    abundanceData.setAtomData("N", 15, 15.0001088989, 0.00368, 1);

    //TODO: need exactly two species, any more not supported.
    abundanceData.setAtomData("S", 32, 31.9720711744, 0.9493, 0);
//    abundanceData.setAtomData("S", 33, 32.9714589099, 0.0076, 1);
    abundanceData.setAtomData("S", 34, 33.96786701, 0.0368, 2);
//    abundanceData.setAtomData("S", 36, 35.96708070, 0.0002, 4);

    abundanceData.atomToDeltaMz.insert(make_pair("C", 1.00335483521));
    abundanceData.atomToDeltaMz.insert(make_pair("H", 1.00627674587));
    abundanceData.atomToDeltaMz.insert(make_pair("O", 2.004244992879));
    abundanceData.atomToDeltaMz.insert(make_pair("N", 0.99703489444));
    abundanceData.atomToDeltaMz.insert(make_pair("S", 1.9957958356));

    return abundanceData;
}();

void NaturalAbundanceData::print() {
    for (auto it = extraNeutronToAtoms.begin(); it != extraNeutronToAtoms.end(); ++it) {
        int neutronNum = it->first;
        vector<Atom> atoms = it->second;
        cout << "[M+" << neutronNum << "]:" << endl;

        for (auto atom : atoms) {
            auto val = atomToAbundance.at(atom);
            cout << atom.massNumber << atom.symbol << ": " << val << endl;
        }
        cout << endl;
    }
}

vector<Atom> NaturalAbundanceData::getAtomsByExtraNeutrons(int numNeutrons) {
    if (extraNeutronToAtoms.find(numNeutrons) != extraNeutronToAtoms.end()) {
        return extraNeutronToAtoms.at(numNeutrons);
    }

    return vector<Atom>{};
}

vector<Atom> NaturalAbundanceData::getAtomsBySymbol(string atomicSymbol) {
    if (symbolToAtoms.find(atomicSymbol) != symbolToAtoms.end()) {
        return symbolToAtoms.at(atomicSymbol);
    }

    return vector<Atom>{};
}

double NaturalAbundanceData::getDeltaMzBySymbol(string atomicSymbol){
    if (atomToDeltaMz.find(atomicSymbol) != atomToDeltaMz.end()) {
        return atomToDeltaMz.at(atomicSymbol);
    }
    return 0;
}

NaturalAbundanceDistribution MassCalculator::getNaturalAbundanceDistribution(
    string compoundFormula,
    Adduct *adduct,
    NaturalAbundanceData& data,
    double minAbundance,
    bool debug) {

    NaturalAbundanceDistribution naturalAbundanceDistribution;

    map<string, int> atoms = getComposition(compoundFormula);
    multiplyAtoms(atoms, adduct->nmol);
    addAtoms(atoms, getComposition(adduct));

    //Issue 703: If the adduct atoms reduce the compound atoms to negative amounts, remove them from the map
    map<string, int> cleanedAtoms{};
    for (auto it = atoms.begin(); it != atoms.end(); ++it){
        if (it->second > 0) {
            cleanedAtoms.insert(make_pair(it->first, it->second));
        }
    }
    atoms = cleanedAtoms;

    unsigned int chgNum = abs(adduct->charge);

    // Compute partial probabilities
    // atomSymbol, numRare, probability
    map<string, map<vector<pair<Atom, int>>, double>> atomTypeToPartialProbability{};

    for (auto it = atoms.begin(); it != atoms.end(); ++it) {
        string atomSymbol = it->first;
        int atomTotal = it->second;

        if (debug) cout << "atomsMap: " << atomSymbol << " n=" << atomTotal << endl;

        vector<Atom> possibleMasses = data.getAtomsBySymbol(atomSymbol);

        //skip any elements that do not have multiple isotopic natural abundances recorded
        // (assume 100% monoisotopic)
        if (possibleMasses.size() <= 1) continue;

        if (possibleMasses.size() != 2) {
            //TODO: support more than one possible type of isotopic shift per atom.
            cerr << "Currenly only supports a maximum of two natural isotopic abundances. Aborting.";
            abort();
        }

        Atom commonIsotope = possibleMasses[0];
        Atom rareIsotope = possibleMasses[1];

        if (data.atomToAbundance.at(commonIsotope) < data.atomToAbundance.at(rareIsotope)) {
            Atom tmp = commonIsotope;
            commonIsotope = rareIsotope;
            rareIsotope = tmp;
        }

        double commonAbundance = data.atomToAbundance.at(commonIsotope);
        double rareAbundance = data.atomToAbundance.at(rareIsotope);

        map<vector<pair<Atom, int>>, double> atomPartialMap{};

        for (unsigned int i = 0; i <= atomTotal; i++) {
            double partialProbability = mzUtils::nchoosek(atomTotal,i)*pow(commonAbundance,atomTotal-i)*pow(rareAbundance,i);
            if (partialProbability >= minAbundance) {

                pair<Atom, int> firstPair = make_pair(commonIsotope, (atomTotal-i));
                pair<Atom, int> secondPair = make_pair(rareIsotope, i);

                vector<pair<Atom, int>> counts{firstPair, secondPair};

                atomPartialMap.insert(make_pair(counts, partialProbability));
                if (debug) {
                    cout << atomSymbol << ", " << i << " isotopes = " << partialProbability << endl;
                }
            }
        }

        atomTypeToPartialProbability.insert(make_pair(atomSymbol, atomPartialMap));

    }

    vector<IsotopicAbundance> existingAbundances{};

    for (auto it = atoms.begin(); it != atoms.end(); ++it) {
        string atomSymbol = it->first;
        int numAtoms = it->second;

        vector<IsotopicAbundance> atomAbundances{};

        if (atomTypeToPartialProbability.find(atomSymbol) == atomTypeToPartialProbability.end()) {

            int massNumber = round(MassCalculator::getElementMass(atomSymbol));
            Atom at(atomSymbol, massNumber);
            IsotopicAbundance isotopicAbundance;
            isotopicAbundance.atomCounts.insert(make_pair(at, numAtoms));
            atomAbundances.push_back(isotopicAbundance);

            if (debug) cout << "[all mono] IsotopicAbundance: " << isotopicAbundance.toString() << endl;

        } else {

            auto atomPartialMap = atomTypeToPartialProbability.at(atomSymbol);

            for (auto it2 = atomPartialMap.begin(); it2 != atomPartialMap.end(); ++it2) {
                vector<pair<Atom, int>> atomVals = it2->first;
                double partialProb = it2->second;
                IsotopicAbundance isotopicAbundance;
                for (auto atomVal : atomVals) {
                    isotopicAbundance.atomCounts.insert(atomVal);
                }
                isotopicAbundance.naturalAbundance *= partialProb;

                if (isotopicAbundance.naturalAbundance >= minAbundance) {
                    atomAbundances.push_back(isotopicAbundance);

                    if (debug) cout << "[atomAbundances] IsotopicAbundance: " << isotopicAbundance.toString() << endl;
                }
            }
        }


        if (debug) {
            cout << atomSymbol << ", before merge:"
                 << " existingAbundances = " << existingAbundances.size()
                 << ", atomAbundances = " << atomAbundances.size()
                 << endl;
        }

        vector<IsotopicAbundance> updatedAbundances;
        if (existingAbundances.empty()) {
            updatedAbundances = atomAbundances;
        } else {
            updatedAbundances = vector<IsotopicAbundance>();
            for (auto existingAbundance : existingAbundances) {
                for (auto atomAbundance : atomAbundances) {
                    IsotopicAbundance combinedAbundance = IsotopicAbundance::createMergedAbundance(existingAbundance, atomAbundance);

                    if (combinedAbundance.naturalAbundance >= minAbundance) {
                        updatedAbundances.push_back(combinedAbundance);
                        if (debug) cout << "Updated Abundances: " << combinedAbundance.toString() << endl;
                    }

                    if (debug) cout << "[existing + atom combination] IsotopicAbundance: " << combinedAbundance.toString() << endl;
                }
            }
        }

        existingAbundances = updatedAbundances;

        if (debug) {
            cout << atomSymbol << ", after merge:"
                 << " existingAbundances = " << existingAbundances.size()
                 << endl;
        }

    }

    for (auto& isotopeAbundance : existingAbundances) {
        isotopeAbundance.compute(data, chgNum);
    }

    if (!existingAbundances.empty()) {
        sort(existingAbundances.begin(), existingAbundances.end(), [](IsotopicAbundance& lhs, IsotopicAbundance& rhs){
            return lhs.mass < rhs.mass;
        });

        IsotopicAbundance monoSpecies = existingAbundances[0];

        for (auto& isotopicAbundance : existingAbundances) {
            isotopicAbundance.naturalAbundanceMonoProportion = isotopicAbundance.naturalAbundance/monoSpecies.naturalAbundance;
        }
    }

    naturalAbundanceDistribution.isotopicAbundances = existingAbundances;

    return naturalAbundanceDistribution;
}

pair<double, double> NaturalAbundanceDistribution::getIsotopicAbundance(Isotope& isotope){
    double naturalAbundance = 0.0;
    double naturalAbundanceMonoProportion = 0.0;

    int numC13 = isotope.C13;
    int numN15 = isotope.N15;
    int numS34 = isotope.S34;
    int numH2 = isotope.H2;
    int numO18 = isotope.O18;

    Atom C13("C", 13);
    Atom N15("N", 15);
    Atom S34("S", 34);
    Atom H2("H", 2);
    Atom O18("O", 18);

    for (auto& isotopicAbundance : isotopicAbundances) {

        //Check C13
        int observedNumC13 = 0;
        if (isotopicAbundance.atomCounts.find(C13) != isotopicAbundance.atomCounts.end()) {
            observedNumC13 = isotopicAbundance.atomCounts[C13];
        }
        if (observedNumC13 != numC13) continue;

        //Check N15
        int observedNumN15 = 0;
        if (isotopicAbundance.atomCounts.find(N15) != isotopicAbundance.atomCounts.end()) {
            observedNumN15 = isotopicAbundance.atomCounts[N15];
        }
        if (observedNumN15 != numN15) continue;

        //Check S34
        int observedNumS34 = 0;
        if (isotopicAbundance.atomCounts.find(S34) != isotopicAbundance.atomCounts.end()) {
            observedNumS34 = isotopicAbundance.atomCounts[S34];
        }
        if (observedNumS34 != numS34) continue;

        //Check H2 (D)
        int observedNumH2 = 0;
        if (isotopicAbundance.atomCounts.find(H2) != isotopicAbundance.atomCounts.end()) {
            observedNumH2 = isotopicAbundance.atomCounts[H2];
        }
        if (observedNumH2 != numH2) continue;

        //Check O18
        int observedNumO18 = 0;
        if (isotopicAbundance.atomCounts.find(O18) != isotopicAbundance.atomCounts.end()) {
            observedNumO18 = isotopicAbundance.atomCounts[O18];
        }
        if (observedNumO18 != numO18) continue;

        //All checks pass - record values and exit loop
        naturalAbundance = isotopicAbundance.naturalAbundance;
        naturalAbundanceMonoProportion = isotopicAbundance.naturalAbundanceMonoProportion;
        break;

    }

    return(make_pair(naturalAbundance, naturalAbundanceMonoProportion));
}

void IsotopicAbundance::compute(NaturalAbundanceData& naturalAbundanceData, unsigned int chgNumber) {
    mass = 0.0;

    for (auto it = atomCounts.begin(); it != atomCounts.end(); ++it) {
        Atom atom = it->first;
        int count = it->second;

        if (count <= 0) continue;

        if (naturalAbundanceData.atomToMass.find(atom) != naturalAbundanceData.atomToMass.end()) {
            mass += naturalAbundanceData.atomToMass[atom] * count;
        } else {
            mass += MassCalculator::getElementMass(atom.symbol) * count;
        }

        if (naturalAbundanceData.atomToNumExtraNeutrons.find(atom) != naturalAbundanceData.atomToNumExtraNeutrons.end()) {
            int numExtraNeutrons = naturalAbundanceData.atomToNumExtraNeutrons.at(atom) * count;
            numTotalExtraNeutrons += numExtraNeutrons;
            if (numExtraNeutrons > 0) {
                labeledForms.insert(atom);
            } else {
                unlabeledForms.insert(atom);
            }
        }

    }

    mz = mass;
    if (chgNumber != 0){
        mz = mass/chgNumber;
    }
}

IsotopicAbundance IsotopicAbundance::createMergedAbundance(IsotopicAbundance& one, IsotopicAbundance& two){
    IsotopicAbundance merged;
    merged.naturalAbundance = one.naturalAbundance * two.naturalAbundance;

    merged.atomCounts = two.atomCounts;

    for (auto it = one.atomCounts.begin(); it != one.atomCounts.end(); ++it) {
        if (merged.atomCounts.find(it->first) == merged.atomCounts.end()) {
            merged.atomCounts.insert(make_pair(it->first, it->second));
        } else {
            merged.atomCounts[it->first]+= it->second;
        }
    }

    return merged;
}

string IsotopicAbundance::getFormula() {
    stringstream s;

    for (auto it = atomCounts.begin(); it != atomCounts.end(); ++it) {
        s << "[" << it->first.massNumber << it->first.symbol << "]" << it->second;
    }

    return s.str();
}

string IsotopicAbundance::toString() {
    stringstream s;

    s << std::fixed << setprecision(5); // 5 places after the decimal

    s << getFormula() << " (" << mass << "): " << naturalAbundance << " " << 100*naturalAbundanceMonoProportion << "% [M+0]";

    return s.str();
}


/**
 * @brief IsotopicAbundance::toIsotope
 * @return IsotopicAbundance information reformatted as Isotope class.
 * This support backwards compatability to Isotope, which is used throughout the MAVEN GUI.
 */
Isotope IsotopicAbundance::toIsotope() {

    Isotope isotope;

    Atom C13("C", 13);
    Atom N15("N", 15);
    Atom S34("S", 34);
    Atom H2("H", 2);
    Atom O18("O", 18);

    stringstream isotopeName;

    //Check C13
    int observedNumC13 = 0;
    if (atomCounts.find(C13) != atomCounts.end()) {
        observedNumC13 = atomCounts[C13];
        if (observedNumC13 > 0) {
            isotopeName << "C13";
        }
    }
    isotope.C13 = observedNumC13;

    //Check N15
    int observedNumN15 = 0;
    if (atomCounts.find(N15) != atomCounts.end()) {
        observedNumN15 = atomCounts[N15];
        if (observedNumN15 > 0) {
            isotopeName << "N15";
        }
    }
    isotope.N15 = observedNumN15;

    //Check S34
    int observedNumS34 = 0;
    if (atomCounts.find(S34) != atomCounts.end()) {
        observedNumS34 = atomCounts[S34];
        if (observedNumS34 > 0) {
            isotopeName << "S34";
        }
    }
    isotope.S34 = observedNumS34;

    //Check H2 (D)
    int observedNumH2 = 0;
    if (atomCounts.find(H2) != atomCounts.end()) {
        observedNumH2 = atomCounts[H2];
        if (observedNumH2 > 0) {
            isotopeName << "D";
        }
    }
    isotope.H2 = observedNumH2;

    //Check O18
    int observedNumO18 = 0;
    if (atomCounts.find(O18) != atomCounts.end()) {
        observedNumO18 = atomCounts[O18];
        if (observedNumO18 > 0) {
            isotopeName << "O18";
        }
    }
    isotope.O18 = observedNumO18;

    // abundance and m/z information makes its way to here
    isotope.abundance = naturalAbundance;
    isotope.naturalAbundanceMonoProportion = naturalAbundanceMonoProportion;
    isotope.mz = mz;

    //name information
    if (isotope.isParent()) {
        isotope.name = "C12 PARENT";
    } else {
        isotopeName << "-label";
        if (isotope.C13 > 0) isotopeName << "-" << isotope.C13;
        if (isotope.N15 > 0) isotopeName << "-" << isotope.N15;
        if (isotope.S34 > 0) isotopeName << "-" << isotope.S34;
        if (isotope.H2 > 0) isotopeName << "-" << isotope.H2;
        if (isotope.O18 > 0) isotopeName << "-" << isotope.O18;

        isotope.name = isotopeName.str();
    }

    return isotope;
  }

bool IsotopicAbundance::isHasAtom(Atom& atom) {
    for (auto it = atomCounts.begin(); it != atomCounts.end(); ++it) {
        if (it->second > 0 && it->first.massNumber == atom.massNumber && it->first.symbol == atom.symbol) return true;
    }

    return false;
}

 // Issue 656: correct quant, based on appropriate naturalAbundanceMonoProportion.
 // Each Isotope and IsotopicAbundance keeps track of its naturalAbundanceMonoProportion.
 // Returns any left over abundance after removing expected natural abundance.
 // If no abundance left over, return 0.
double MassCalculator::getNaturalAbundanceCorrectedQuantValue(
      float isotopeObserved, float mZeroObserved,
      double isotopeExpectedAbundance, double mZeroExpectedAbundance) {

    //Do not attempt correction if the [M+0] has an expected abundance outside of valid range.
    if (mZeroExpectedAbundance <= 0 || mZeroExpectedAbundance > 1.0){
        return isotopeObserved;
    }

    float correctedVal = isotopeObserved - (isotopeExpectedAbundance/mZeroExpectedAbundance)*mZeroObserved;

    //Only return corrected vals if the value is strictly greater than 0 - negative values are ignored.
    if (correctedVal > 0) {
        return correctedVal;
    }

    return 0;
}

//Issue 758
vector<pair<Adduct, map<MassAtom, int>>> MassCalculator::candidateAtomMaps(
    double mz,
    vector<Adduct> possibleAdducts,
    map<MassAtom, pair<int, int>> legalAtomCounts,
    double ppmDiff,
    long maxNumCandidates,
    bool debug
    ) {

    //initialize output
    vector<pair<Adduct, map<MassAtom, int>>> atomMapCandidates{};

    double minMz = mz - mz*ppmDiff/1e6;
    double maxMz = mz + mz*ppmDiff/1e6;

    if (debug) {
        cout << "Enumerating valid atom maps with mz= [" << minMz << ", " << maxMz << "]" << endl;
    }

    //Determine which atoms can legally stand alone (that is, do not require any other atoms)
    //For atom X, if any atoms in the map occur before atom X with a non-zero requirement, atom X cannot stand on its own.
    vector<MassAtom> cannotStandAlone{};
    for (auto it = legalAtomCounts.begin(); it != legalAtomCounts.end(); ++it) {
        if (it->second.first > 0) {

            auto nextIt = std::next(it);

            if (nextIt != legalAtomCounts.end()) {
                //all subsequent atoms cannot stand on their own.
                std::transform(nextIt, legalAtomCounts.end(), std::back_inserter(cannotStandAlone),
                               [](const pair<MassAtom, pair<int, int>> pair) { return pair.first; });
            }

            break;
        }
    }

    set<string> cannotStandAloneSymbols{};
    for (MassAtom & massAtom : cannotStandAlone) {
        cannotStandAloneSymbols.insert(massAtom.symbol);
    }

    if (debug) {
        cout << "The following atoms cannot stand alone: ";
        for (MassAtom & massAtom : cannotStandAlone) {
            cout << massAtom.symbol << " ";
        }
        cout << endl;
    }

    //Wrap everything in the adducts map - this enumeration is then solving for the "M" part of the
    //adduct, we don't care about the other stuff
    for (auto& adduct : possibleAdducts) {

        double minParentMass = min(adduct.computeParentMass(maxMz), adduct.computeParentMass(minMz));
        double maxParentMass = max(adduct.computeParentMass(maxMz), adduct.computeParentMass(minMz));

        if (debug) {
            cout << "Enumerating candidates for parent mass, assuming adduct="
                 << adduct.name << "; parent range=[" << minParentMass << ", " << maxParentMass << "]" << endl;
        }

        //Enumerate all possible candidates for this "M".
        //While enumerating, immediately stop if the M so far if the mass is too high.
        vector<map<MassAtom, int>> candidatesMap{};

        //iterate through every possible number of non-carbon atoms.
        for (auto it = legalAtomCounts.begin(); it != legalAtomCounts.end(); ++it) {

            const MassAtom& atomType = it->first;

            int minNum = it->second.first;
            int maxNum = it->second.second;
            int N = maxNum-minNum+1;

            //Case: do not add any of the new atoms.
            //Note that this is only legal if minNum is zero.
            vector<map<MassAtom, int>> mapWithNoNewAdditions{};
            if (minNum == 0) {
                mapWithNoNewAdditions = candidatesMap;
            }

            //Case: only add the new atoms, do not keep any old atoms.
            //check to see if this is legal first
            vector<map<MassAtom, int>> newAdditionsOnly(N);
            unsigned long counter = 0;
            bool isStandAlone = true;
            if (cannotStandAloneSymbols.find(atomType.symbol) != cannotStandAloneSymbols.end()) {
                isStandAlone = false;
                newAdditionsOnly.clear();
            }

            vector<map<MassAtom,int>> combinationMap{};
            for (unsigned int i = minNum; i <= maxNum; i++) {

                if (i == 0) continue;

                if (isStandAlone) {
                    map<MassAtom, int> newMap = map<MassAtom, int>();
                    newMap.insert(make_pair(atomType, i));
                    newAdditionsOnly[counter] = newMap;
                }

                //Case: every atom that is already in the map is updated with the new additions.
                vector<map<MassAtom, int>> ithDuplicate = candidatesMap;
                for (auto & candidate : ithDuplicate) {
                    candidate.insert(make_pair(atomType, i));
                    combinationMap.push_back(candidate);
                }

                counter++;
            }

            unsigned int updatedSize = mapWithNoNewAdditions.size() + newAdditionsOnly.size() + combinationMap.size();
            vector<map<MassAtom, int>> updatedCandidatesMap;
            updatedCandidatesMap.reserve(updatedSize);

            updatedCandidatesMap.insert(updatedCandidatesMap.end(), mapWithNoNewAdditions.begin(), mapWithNoNewAdditions.end());
            updatedCandidatesMap.insert(updatedCandidatesMap.end(), newAdditionsOnly.begin(), newAdditionsOnly.end());
            updatedCandidatesMap.insert(updatedCandidatesMap.end(), combinationMap.begin(), combinationMap.end());

            //prepare for next iteration
            candidatesMap.clear();

            //filtering for next iteration - if the m/z is already too high, abort.
            for (auto & candidate : updatedCandidatesMap) {
                double parentMassTotal = 0.0;
                for(auto it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
                    parentMassTotal += it2->second* it2->first.massValue;
                    if (parentMassTotal > maxParentMass) break;
                }
                if (parentMassTotal <= maxParentMass) {
                    candidatesMap.push_back(candidate);
                    if (candidatesMap.size() > maxNumCandidates) {
                        if (debug) {
                            cerr << "Too many candidates to enumerate - reduce m/z or modify legal atom counts argument.\n"
                                 << "Returning an empty vector."
                                 << endl;
                        }
                        return vector<pair<Adduct, map<MassAtom, int>>>();
                    }
                }
            }

            if (debug) {
                cout << "After introduction of " << atomType.symbol << ": "
                     << candidatesMap.size() << " atom map candidates. "
                     << endl;
            }
        }

        for (auto & candidate : candidatesMap) {
            double parentMassTotal = 0.0;
            for(auto it = candidate.begin(); it != candidate.end(); ++it) {
                parentMassTotal += it->second * it->first.massValue;
            }
            if (parentMassTotal <= maxParentMass && parentMassTotal >= minParentMass) {
                atomMapCandidates.push_back(make_pair(adduct, candidate));
            }
        }
    }

    return atomMapCandidates;
}

vector<pair<string, double>> MassCalculator::evaluateAtomMapCandidates(
    double mz,
    const vector<pair<Adduct, map<MassAtom, int>>>& atomMapCandidates,
    bool debug) {

    vector<pair<string, double>> evaluatedCandidates(atomMapCandidates.size());

    for (unsigned int i = 0; i < atomMapCandidates.size(); i++) {

        Adduct adduct = atomMapCandidates[i].first;
        map<MassAtom, int> atomMap = atomMapCandidates[i].second;

        double parentMassTotal = 0;

        stringstream s;
        for(auto it = atomMap.begin(); it != atomMap.end(); ++it) {
            parentMassTotal += it->second * it->first.massValue;

            s << it->first.symbol;
            if (it->second > 1) {
                s << it->second;
            }
        }

        double candidateMz = adduct.computeAdductMass(parentMassTotal);

        double ppmDiff = (candidateMz-mz)/candidateMz*1e6;
        string formulaString = s.str();

        evaluatedCandidates[i] = make_pair(formulaString, ppmDiff);
    }

    sort(evaluatedCandidates.begin(), evaluatedCandidates.end(), [](pair<string, double>& lhs, pair<string, double>& rhs) {
        return abs(lhs.second) < abs(rhs.second);
    });

    return evaluatedCandidates;
}
