//dump basic scan informatio
#include "Peptide.h"
#include <iostream>
#include <iomanip>

bool SHUFFLE=false;

void printTheoryFragmentaiton(Peptide* pep) {

	if (SHUFFLE) { 
		pep=pep->shuffle(0.5, 0);
		cout << "DECOY_";	
	}

	cout << "PEPTIDE:\"" << pep->interactStyleWithCharge() 
			     << "\"\tpreMz:" 
	     		     << setprecision(7) 
			     << pep->monoisotopicMZ() << endl;

	vector<FragmentIon*> theoreticalIons;
	pep->generateFragmentIonsCID(theoreticalIons);

	for(FragmentIon* f: theoreticalIons ) {
		if (f->m_prominence < 4 or f->m_charge >= 3) continue;
		cout << "\t" << f->m_mz << "\t" << f->m_ion << endl;
	}

	if(theoreticalIons.size()) { 
		for(int i=0; i < theoreticalIons.size(); i++ ) delete theoreticalIons[i];
	}
}



int main(int argc, char** argv) {

    for(int i=1; i< argc; i++ ) {
        string opt(argv[i]);
	if (opt == "-shuffle") SHUFFLE=true;
    }

    string peptideString;
    while(cin) {
        getline(cin, peptideString);
        Peptide* pep = new Peptide(peptideString,0);
	printTheoryFragmentaiton(pep);
	delete pep;
    };

    
}
