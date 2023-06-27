#include "proteinutils.h"

FastaWritable::~FastaWritable(){}

/**
 * @brief Protein::Protein constructor
 * @param header
 * @param seq
 */
Protein::Protein(string header, string seq){
    this->header = header;
    this->seq = seq;

    this->mw = ProteinUtils::getProteinMass(seq);
}

void Protein::printSummary() {
    cout << header << " [length: " << seq.size() << " AA, MW: " << mw << " Da]" << endl;
}

Protein::~Protein(){}

/**
 * @brief ProteinFragment::ProteinFragment constructor
 * @param protein
 * @param theoreticalMw
 * @param observedMw
 * @param start
 * @param end
 */
ProteinFragment::ProteinFragment(Protein* protein, double theoreticalMw, double observedMw, unsigned long start, unsigned long end){
    this->protein = protein;
    this->theoreticalMw = theoreticalMw;
    this->observedMw = observedMw;
    this->start = start;
    this->end = end;
    this->deltaMw = abs(theoreticalMw-observedMw);
}

string ProteinFragment::getSequence() const {
    return protein->seq.substr(start, (end-start+1));
}

string ProteinFragment::getHeader() const {
    stringstream s;
    s << std::fixed << setprecision(3)
      << protein->header
      << " FRAGMENT ["
      << "seq: " << (start+1) << " - " << (end+1)
      << ", theo mw: " << theoreticalMw << " Da"
      << ", observed mw: " << observedMw << " Da"
      << ", delta mw: " << deltaMw << " Da"
      << "]";
    return s.str();
}

ProteinFragment::~ProteinFragment(){}

/**
 * @brief ProteinUtils
 *
 * all static methods
 */

double ProteinUtils::getProteinMass(string seq) {
    double proteinMass = 0.0;

    for (char aa : seq) {

        //silently ignore any weird characters in the sequence, e.g. 'X' or 'N'.
        if (aaMasses.find(aa) != aaMasses.end()) {
            proteinMass += aaMasses[aa];
        }
    }

    return proteinMass;
}

vector<Protein*> ProteinUtils::loadFastaFile(string filename) {

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

            header = line.substr(1); // remove starting '>' character
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

void ProteinUtils::writeFastaFile(vector<FastaWritable*> entries, string outputFile, unsigned int seqLineMax) {
    ofstream outputFileStream;
    outputFileStream.open(outputFile);

    for (FastaWritable *entry : entries) {
        outputFileStream << ">" << entry->getHeader() << "\n";

        string::size_type N = entry->getSequence().size();
        string::size_type currentPos = 0;

        while (currentPos != string::npos) {

            if (currentPos+seqLineMax < N) {

                outputFileStream << entry->getSequence().substr(currentPos, seqLineMax)
                                 << "\n";

                currentPos = currentPos + seqLineMax;
            } else {
                string::size_type remainder = N-currentPos;
                outputFileStream << entry->getSequence().substr(currentPos, remainder) << "\n";
                break;
            }
        }

        outputFileStream << "\n";
    }

    outputFileStream.close();
}
