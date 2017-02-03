//dump basic scan informatio
#include "Peptide.h"
#include <iostream>


int main(int argc, char** argv) {

    for(int i=1; i< argc; i++ ) {

        string peptideString(argv[i]);
        Peptide pep(peptideString,0);

        vector<FragmentIon*> theoreticalIons;
        pep.generateFragmentIonsCID(theoreticalIons);

        cerr << "PEPTIDE: " << peptideString << "\tpreMz:" << pep.monoisotopicMZ() << endl;

        for(FragmentIon* f: theoreticalIons ) {
            if (f->m_prominence < 5) continue;
            cerr << f->m_mz << "\t" << f->m_ion << endl;
        }
    }
}
