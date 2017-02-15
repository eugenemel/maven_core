#include "Peptide.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

//globals
//
int MINLENGTH=5;
int MAXLENTH=20;
int MAXMISSED=3;

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

        void calculateMZ(string seq) {

            //charge states
            for(int z=2; z<5; z++) {
                string peptideString(seq);
                Peptide pep(peptideString,z);
                cout << pep.fullWithCharge() << pep.monoisotopicMZ() << endl;
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
    }

    //
    //load fasta file
    //
    if(!fastafile.empty()) loadFastaFile(fastafile);
 

    //
    //digest
    //
    for(Protein* p: proteins) { 
        p->trypticDigest();
    }
}
