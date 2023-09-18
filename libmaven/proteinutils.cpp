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

/**
 * @brief Protein::fragmentProtein
 *
 * Progressively examine the pieces of a protein made by cutting at two places,
 * 'cut1', from the C terminus, and 'cut2', from the N terminus.
 *
 *
 *       M1     M2   M3
 *  N ------|------|------------ C
 *          cut1   cut2
 *
 * Return sorted vector of ProteinFragment*, in order by observed mass, then descending on
 * delta in mw from theoretical to observed.
 *
 * @param fragMasses
 * @param tolerance
 * @param debug
 * @return
 */
vector<ProteinFragment*> Protein::fragmentProtein(vector<double>& fragMasses, double tolerance, bool debug) {

    sort(fragMasses.begin(), fragMasses.end(), [](const double& lhs, const double& rhs){
        return lhs < rhs;
    });

    //initialize output
    vector<ProteinFragment*> fragmentProteins{};

    //initial state
    double m1 = 0.0;
    double m2 = this->mw;
    double m2Cut2 = m2;
    double m3 = 0.0;

    for (unsigned long i = 0; i < this->seq.length(); i++) {

        //update temp variables in preparation for iteration
        m2Cut2 = m2;

        for (unsigned long j = this->seq.length()-1; j > i; j--) {

            //second cut affects m2 and m3 values
            m2Cut2 -= aaMasses.at(this->seq[j]);
            m3 += aaMasses.at(this->seq[j]);

            for (double mass : fragMasses) {

                if (abs(m1-mass) < tolerance && i == 0) {
                    ProteinFragment *fragment = new ProteinFragment(this, mass, m1, 0, (i-1));
                    fragmentProteins.push_back(fragment);
                    if (debug) cout << "(" << (i-1) << ", " << 0 << "): mass=" << mass << ", M1=" << m1 << "Da, delta=" << abs(m1-mass) << endl;
                }

                if (abs(m2Cut2-mass) < tolerance) {
                    ProteinFragment *fragment = new ProteinFragment(this, mass, m2Cut2, i, (j-1));
                    fragmentProteins.push_back(fragment);
                    if (debug) cout << "(" << i << ", " << (j-1) << "): mass="<< mass << ", M2=" << m2Cut2 << " Da, delta=" << abs(m2Cut2-mass) << endl;
                }

                if (abs(m3-mass) < tolerance && i == 0) {
                    ProteinFragment *fragment = new ProteinFragment(this, mass, m3, j, this->seq.length()-1);
                    fragmentProteins.push_back(fragment);
                    if (debug) cout << "(" << i << ", " << (this->seq.length()-1) << "): mass="<< mass << ", M3=" << m3 << " Da, delta=" << abs(m3-mass) << endl;
                }
            }
        }

        //first cut affects m1 and m2 values
        m1 += aaMasses.at(this->seq[i]);
        m2 -= aaMasses.at(this->seq[i]);
        m3 = 0.0;
    }

    sort(fragmentProteins.begin(), fragmentProteins.end(), [](ProteinFragment* lhs, ProteinFragment* rhs){
        if (lhs->theoreticalMw == rhs->theoreticalMw) {
            return lhs->deltaMw < rhs->deltaMw;
        } else {
            return lhs->theoreticalMw < rhs->theoreticalMw;
        }
    });

    return fragmentProteins;

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

    //add water
    proteinMass += 18.0105; // adds to the C terminus end

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
