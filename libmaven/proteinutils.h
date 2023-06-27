#ifndef PROTEINUTILS_H
#define PROTEINUTILS_H

#include "mzSample.h"

static map<char, double> aaMasses = {
    {'A', 71.03711},
    {'R', 156.10111},
    {'N', 114.04293},
    {'D', 115.02694},
    {'C', 103.00919},
    {'E', 129.04259},
    {'Q', 128.05858},
    {'G', 57.02146},
    {'H', 137.05891},
    {'I', 113.08406},
    {'L', 113.08406},
    {'K', 128.09496},
    {'M', 131.04049},
    {'F', 147.06841},
    {'P', 97.05276},
    {'S', 87.03203},
    {'T', 101.04768},
    {'W', 186.07931},
    {'Y', 163.06333},
    {'V', 99.06841},
    {'U', 150.95363},
};

class FastaWritable;
class Protein;
class ProteinFragment;

class FastaWritable {
public:
    virtual string getHeader() const = 0;
    virtual string getSequence() const = 0;

    virtual ~FastaWritable() = 0;
};

class Protein : public FastaWritable{
public:
        string header;
        string seq;
        double mw;

        Protein(string header, string seq);

        void printSummary();
        //vector<ProteinFragment*> fragmentProtein(vector<double>& fragMasses, double tolerance, bool debug);

        string getHeader() const {return header;}
        string getSequence() const {return seq;}

        ~Protein();

};

//Use composition instead of inheritance here
class ProteinFragment : public FastaWritable {
public:
    Protein *protein;
    unsigned long start; // 0-indexed
    unsigned long end; // 0-indexed

    double theoreticalMw;
    double observedMw;

    ProteinFragment(Protein* protein, double theoreticalMw, double observedMw, unsigned long start, unsigned long end);

    double deltaMw;

    string getHeader() const;
    string getSequence() const;

    ~ProteinFragment();
};

class ProteinUtils {
public:
    static double getProteinMass(string seq);
    static void writeFastaFile(vector<FastaWritable*>, string outputFile, unsigned int seqLineMax = 87);
    static vector<Protein*> loadFastaFile(string fastaFile);
};

#endif // PROTEINUTILS_H
