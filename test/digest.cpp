#include "Peptide.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

//globals
//
int  MINLENGTH=5;
int  MAXLENTH=30;
int  MAXMISSED=3;
bool TMT_LABELED=false;
bool PRINT_IONS=false;
bool BUILD_INDEX=false;
bool TRYPTIC_DIGEST=false;

using namespace std;

class Protein { 

    public:
        string header; 
        string seq; 

        Protein(string h, string s) { 
            header=h; 
            seq=s; 
            seq.erase(std::remove_if(seq.begin(), seq.end(), [](char x){return std::isspace(x);}),seq.end());
        }

        void trypticDigest() { 
            vector<int>cuts;

            //find cut points
            cuts.push_back(0);
            for(int i=0; i<seq.length(); i++ ) {
                if(seq[i] == 'R' or seq[i] == 'K')  cuts.push_back(i); 
            }
            cuts.push_back(seq.length()-1);


            for(int i=0; i < cuts.size()-1; i++ ) {
                int pos1 =  cuts[i];
                for(int j=i+1; j <= cuts.size(); j++ ) {
                    int pos2 =  cuts[j];
                    if (pos2-pos1 > MAXLENTH) break;
                    if (j-i > MAXMISSED) break;

                    if (pos2-pos1 >=MINLENGTH and pos2-pos1 <= MAXLENTH) { 

                        string peptideSeq;
                        if(i>0) peptideSeq = seq.substr(pos1+1,pos2-pos1);
                        else    peptideSeq = seq.substr(pos1+1,pos2-pos1);

                       calculateMZ(peptideSeq);

                    }
                }
            }
        }

	void printTheoryFragmentaiton(Peptide& pep) {

                cout << "PEPTIDE:\"" << pep.interactStyleWithCharge() << "\"\tpreMz:" << setprecision(7) << pep.monoisotopicMZ() << endl;

		vector<FragmentIon*> theoreticalIons;
		pep.generateFragmentIonsCID(theoreticalIons);

		for(FragmentIon* f: theoreticalIons ) {
			if (f->m_prominence < 4 or f->m_charge >= 3) continue;
			cout << "\t" << f->m_mz << "\t" << f->m_ion << endl;
		}

		if(theoreticalIons.size()) { 
			for(int i=0; i < theoreticalIons.size(); i++ ) delete theoreticalIons[i];
		}
	}


        void calculateMZ(string seq) {

            string modPeptideString;
            if (TMT_LABELED) modPeptideString = "n[230]";

            //fixed mods
            for(int i=0; i<seq.size(); i++ ) {
                if      (seq[i] == 'C') modPeptideString += "C[160]";
                else if (TMT_LABELED and seq[i] == 'K')  { modPeptideString += "K[357]"; }
                else  modPeptideString += seq[i];
            }
            //variable mods

            Peptide pep(modPeptideString,1);

            auto SAVE_COPY=pep.mods;
            for(int z=2; z<5; z++) {
                pep.charge = z;
		if(PRINT_IONS) printTheoryFragmentaiton(pep);
		else cout << pep.interactStyleWithCharge() << endl;

                //add oxidation
                for(int i=0; i<seq.size(); i++ ) {
                    if (seq[i] == 'M') { 
                        pep.mods[i] = "Oxidation";
			if(PRINT_IONS) printTheoryFragmentaiton(pep);
			else cout << pep.interactStyleWithCharge() << endl;

                    }
                }
                pep.mods = SAVE_COPY;
            }
        }
};

std::vector<Protein*> proteins;

void loadFastaFile(string filename) { 

    ifstream fastafile(filename);
    if(! fastafile.is_open() ) { cerr << "Can't open file " << filename; return; }

    string line;
    string header;
    string sequence;

    while ( getline(fastafile,line) ) {
        if (line[0] == '>') { 
            if(!header.empty() and !sequence.empty()) { 
                proteins.push_back( new Protein(header,sequence)); 
            } 

            header = line; 
            sequence=string(); 
        } else sequence += line;
    }

    if(!header.empty() and !sequence.empty()) { 
      proteins.push_back( new Protein(header,sequence)); 
    }
}

int main(int argc, char** argv) {

    //
    // process options
    //
    string fastafile;
    for(int i=1; i< argc; i++ ) {
        string optionString(argv[i]);
        if (optionString == "-i" and i+1<argc) fastafile=string(argv[i+1]);
        if (optionString == "-tmt") TMT_LABELED=true;
        if (optionString == "-ions") PRINT_IONS=true;
        if (optionString == "-digest") TRYPTIC_DIGEST=true;
        if (optionString == "-index")  BUILD_INDEX=true;
    }

    //
    //load fasta file
    //
    if(!fastafile.empty()) loadFastaFile(fastafile);
 
    //
    //digest
    //
    if (TRYPTIC_DIGEST) { 
    	for(Protein* p: proteins)  p->trypticDigest();
    }
    
   //
   //build index
   //
   if (BUILD_INDEX) { 
    	for(Protein* p: proteins)  p->trypticDigest();
   }
}
