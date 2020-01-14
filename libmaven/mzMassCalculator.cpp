#include "mzMassCalculator.h"

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

/*---------- function to get molar weight of an element or group --------*/
double MassCalculator::getElementMass(string elmnt){

    /* default behavior is to ignore string */
    double val_atome(0);

/* Check for atoms */
    if (elmnt == "H")        val_atome = 1.0078250321;
    else if (elmnt == "D")   val_atome = 2.01410178;
    else if (elmnt == "C")   val_atome = 12.00000000;
    else if (elmnt == "N")   val_atome = 14.0030740052;
    else if (elmnt == "O")   val_atome = 15.9949146221;
    else if (elmnt == "S")   val_atome = 31.97207069;
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
    return(val_atome);
}
/*-----------------------------------------------------------------------*/

map<string,int> MassCalculator::getComposition(Adduct* adduct){

    cout << "MassCalculator::getComposition(Adduct*)" << endl;

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

    //debugging
    cout << "FORMULAS TO ADD:" << endl;
    for (string posFormula : formulasToAdd) {
        cout << posFormula << endl;
        addAtoms(atoms, getAdductComponentComposition(posFormula));
    }
    cout << endl;

    cout << "FORMULAS TO SUBTRACT:" << endl;
    for (string negFormula : formulasToSubtract) {
        cout << negFormula << endl;
        subtractAtoms(atoms, getAdductComponentComposition(negFormula));
    }
    cout << endl;

    return atoms;
}

void MassCalculator::addAtoms(map<string, int>& reference, map<string, int> toAdd) {
    modifyAtoms(reference, toAdd, true);
}

void MassCalculator::subtractAtoms(map<string, int>& reference, map<string, int> toSubtract) {
    modifyAtoms(reference, toSubtract, false);
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
    } else if (formula == "AcOH" || formula == "AcO") { //acetate (C2H4O2)
        return map<string, int>{{"C", 2}, {"H", 4}, {"O", 2}};
    } else if (formula == "NaOAc") {
        return map<string, int>{{"C", 2}, {"H", 4}, {"O", 2}, {"Na", 1}};
    } else if (formula == "2Na") {
        return map<string, int>{{"Na", 2}};
    } else if (formula == "3Na") {
        return map<string, int>{{"Na", 3}};
    } else if (formula == "2K") {
        return map<string, int>{{"K", 2}};
    } else if (formula == "3K") {
        return map<string, int>{{"K", 3}};
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
	for( itr = atoms.begin(); itr != atoms.end(); itr++ ) {
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

vector<Isotope> MassCalculator::computeIsotopes(string compoundFormula, Adduct* adduct, int maxNumProtons, bool isUse13C, bool isUse15N, bool isUse34S, bool isUse2H) {

    const double abC12 = 0.9893;
    const double abC13 = 0.0107;
    const double abN14 = 0.9963620;
    const double abN15=  0.0036420;
    const double abS32=  0.949926;
    const double abS34=  0.042524;
    const double abH  =  0.999885;
    const double abH2  =  0.00011570;
    //const double abO16  = 0.9975716;
    //const double abO18  = 0.0020514;

    const double D_Delta = 2.01410177811-1.00782503224;
    const double C_Delta = 13.00335483521-12.0;
    const double N_Delta = 15.0001088989-14.00307400446;
    const double S_Delta = 33.96786701-31.9720711744;

    map<string, int> atoms = getComposition(compoundFormula);
    addAtoms(atoms, getComposition(adduct));

    //note that this already includes any mass adjustment from the # of electrons
    double parentMass = adduct->computeAdductMass(computeNeutralMass(compoundFormula));

    int CatomCount  =  max(atoms["C"], 0);
    int NatomCount  =  max(atoms["N"], 0);
    int SatomCount  =  max(atoms["S"], 0);
    int HatomCount  =  max(atoms["H"], 0);

     vector<Isotope> isotopes;

     Isotope parent("C12 PARENT", parentMass);
     isotopes.push_back(parent);

    if (isUse13C){
        for (int i=1; i <= CatomCount; i++ ) {
                if (i > maxNumProtons) break;
                Isotope x("C13-label-"+integer2string(i), parentMass + (i*C_Delta),i,0,0,0);
                isotopes.push_back(x);
        }
    }

    if (isUse15N) {
        for (int i=1; i <= NatomCount; i++ ) {
                if (i > maxNumProtons) break;
                Isotope x("N15-label-"+integer2string(i), parentMass + (i*N_Delta),0,i,0,0);
                isotopes.push_back(x);
        }
    }

    if (isUse34S) {
        for (int i=1; i <= SatomCount; i++ ) {
                if (i > maxNumProtons) break;
                Isotope x("S34-label-"+integer2string(i), parentMass + (i*S_Delta),0,0,i,0);
                isotopes.push_back(x);
        }
    }

    if (isUse2H) {
        for (int i=1; i <= HatomCount; i++ ) {
                if (i > maxNumProtons) break;
                Isotope x("D-label-"+integer2string(i), parentMass + (i*D_Delta),0,0,0,i);
                isotopes.push_back(x);
        }
    }


    if (isUse13C && isUse15N) {
        for (int i=1; i <= CatomCount; i++ ) {
            for (int j=1; j <= NatomCount; j++ ) {
                if ((i+j) > maxNumProtons) break;
                string name ="C13N15-label-"+integer2string(i)+"-"+integer2string(j);
                double mass = parentMass + (j*N_Delta) + (i*C_Delta);
                    Isotope x(name,mass,i,j,0,0);
                    isotopes.push_back(x);
            }
        }
    }

    if (isUse13C && isUse34S) {
        for (int i=1; i <= CatomCount; i++ ) {
            for (int j=1; j <= SatomCount; j++ ) {
                if ((i+j) > maxNumProtons) break;
                string name ="C13S34-label-"+integer2string(i)+"-"+integer2string(j);
                double mass = parentMass + (j*S_Delta) + (i*C_Delta);
                Isotope x(name,mass,i,0,j,0);
                isotopes.push_back(x);
            }
        }
    }

    if (isUse13C && isUse2H) {
        for (int i=1; i <= CatomCount; i++ ) {
            for (int j=1; j <= HatomCount; j++ ) {
                if ((i+j) > maxNumProtons) break;
                string name ="C13D-label-"+integer2string(i)+"-"+integer2string(j);
                double mass = parentMass + (j*D_Delta) + (i*C_Delta);
                Isotope x(name,mass,i,0,0,j);
                isotopes.push_back(x);
            }
        }
    }

    for(unsigned int i=0; i < isotopes.size(); i++ ) {
           Isotope& x = isotopes[i];
                int c=x.C13;
                int n=x.N15;
                int s=x.S34;
                int d=x.H2;

		isotopes[i].abundance=
                 mzUtils::nchoosek(CatomCount,c)*pow(abC12,CatomCount-c)*pow(abC13,c)
               * mzUtils::nchoosek(NatomCount,n)*pow(abN14,NatomCount-n)*pow(abN15,n)
               * mzUtils::nchoosek(SatomCount,s)*pow(abS32,SatomCount-s)*pow(abS34,s)
               * mzUtils::nchoosek(HatomCount,d)*pow(abH,HatomCount-d)  *pow(abH2,d);
    }

    return isotopes;
}

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

std::string MassCalculator::prettyName(int c, int h, int n, int o, int p, int s) {
		char buf[1000];
		string name;
		if ( c != 0 ) { name += "C"; if (c>1) { sprintf(buf,"%d",c);  name += buf;} }
		if ( h != 0 ) { name += "H"; if (h>1) { sprintf(buf,"%d",h);  name += buf;} }
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
