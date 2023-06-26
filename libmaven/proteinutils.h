#ifndef PROTEINUTILS_H
#define PROTEINUTILS_H

#include "mzSample.h"

class Protein {
public:
        string header;
        string seq;

        Protein(string header, string seq);

        static vector<Protein*> loadFastaFile(string fastaFile);
        void writeFastaFile(vector<Protein*> proteins);
};

#endif // PROTEINUTILS_H
