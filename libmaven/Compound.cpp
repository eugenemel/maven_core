#include "Fragment.h"
#include "mzSample.h"

Compound::Compound(string id, string name, string formula, int charge ) {
		this->id = id;
		this->name = name;
		this->formula = formula;
		this->charge = charge;
        this->exactMass =  mcalc->computeNeutralMass(formula);
		this->expectedRt = -1;
        this->logP = 0;
        this->isDecoy=false;

        cid=0;
        ionizationMode=0;
		precursorMz=0;
		productMz=0;
        collisionEnergy=0;
        virtualFragmentation=false;

}

MassCalculator* Compound::mcalc = new MassCalculator();

float Compound::ajustedMass(int charge) { 
	return Compound::mcalc->computeMass(formula,charge); 
}


FragmentationMatchScore Compound::scoreCompoundHit(Fragment* f, float productPpmTolr=20, bool searchProton=true) {
        FragmentationMatchScore s;
        Compound* cpd = this;

        if (!cpd or cpd->fragment_mzs.size() == 0) return s;

        Fragment t;
        t.precursorMz = cpd->precursorMz;
        t.mzs = cpd->fragment_mzs;
        t.intensity_array = cpd->fragment_intensity;
        t.annotations = cpd->fragment_iontype;

        if (searchProton)  { //special case, check for loss or gain of protons
            int N = t.mzs.size();
            for(int i=0; i<N;i++) {
                t.mzs.push_back( t.mzs[i] + PROTON);
                t.intensity_array.push_back( t.intensity_array[i] );
                t.mzs.push_back( t.mzs[i] - PROTON );
                t.intensity_array.push_back( t.intensity_array[i] );
            }
        }

        //theory fragmentation or library fragmentation = t
        //experimental data = f

//        cerr << endl;
//        cerr << "Compound: " << cpd->name << endl;
//        for (float mz : t.mzs) {
//            cerr << mz << endl;
//        }
//        cerr << endl;

//        cerr << "Observed (precursorMz=" << f->precursorMz << ")" << endl;
//        for (float mz : f->mzs) {
//            cerr << mz << endl;
//        }
//        cerr << endl;

        s = t.scoreMatch(f,productPpmTolr);

        return s;
}

vector<Compound*> Compound::getChildren() {return vector<Compound*>(0);}
vector<Compound*> SummarizedCompound::getChildren() {return children;}

/**
 * @brief SummarizedCompound::computeFragments
 *
 * All m/z values from all compounds are used.
 * Intensity values are averaged based on all m/zs with some intensity value.
 *
 */
void SummarizedCompound::computeFragments() {

    map<int, vector<float>> intensitiesByMz = {};

    for (auto compound : getChildren()) {
        for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

            int mzKey = mzUtils::mzToIntKey(compound->fragment_mzs[i], 1000);
            float intensity = compound->fragment_intensity[i];

            if (intensitiesByMz.find(mzKey) == intensitiesByMz.end()){
                vector<float> intensities(1);
                intensities[0] = intensity;
                intensitiesByMz.insert(make_pair(mzKey, intensities));
            } else {
                intensitiesByMz[mzKey].push_back(intensity);
            }
        }
    }

    fragment_mzs = vector<float>(intensitiesByMz.size());
    fragment_intensity=vector<float>(intensitiesByMz.size());

    unsigned int vecCounter = 0;
    for (map<int, vector<float>>::iterator it = intensitiesByMz.begin(); it != intensitiesByMz.end(); ++it) {

        int mzInt = it->first;
        vector<float> intensities = it->second;

        float avgInt = 0.0f;
        for (auto intensity : intensities) {
            avgInt += intensity;
        }
        avgInt /= intensities.size();

        fragment_mzs[vecCounter] = mzUtils::intKeyToMz(mzInt, 1000);
        fragment_intensity[vecCounter] = avgInt;
        vecCounter++;
    }
}
