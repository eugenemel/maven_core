#include "proteinutils.h"

Protein::Protein(string header, string seq){
    this->header = header;
    this->seq = seq;
}

vector<Protein*> Protein::loadFastaFile(string filename) {

    vector<Protein*> proteins{};

    ifstream fastafile(filename);
    if(! fastafile.is_open() ) {
        cerr << "Can't open file " << filename;
        return proteins;
    }

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
        } else {
            sequence += line;
        }
    }

    //last entry in file
    if(!header.empty() and !sequence.empty()) {
      proteins.push_back(new Protein(header,sequence));
    }

    return proteins;
}

void Protein::writeFastaFile(vector<Protein *> proteins) {
    //TODO
}
