#include "Fragment.h"
#include "mzSample.h"

Compound::Compound(string id, string name, string formula, int charge) {
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

Compound::Compound(string id, string name, string formula, int charge, float exactMass) {
    this->id = id;
    this->name = name;
    this->formula = formula;
    this->charge = charge;
    this->exactMass =  exactMass;
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
        t.fragment_labels = cpd->fragment_labels;

        if (searchProton)  { //special case, check for loss or gain of protons
            int N = t.mzs.size();
            for(int i=0; i<N;i++) {

                t.mzs.push_back( t.mzs[i] + PROTON);
                t.intensity_array.push_back( t.intensity_array[i] );
                t.fragment_labels.push_back(t.fragment_labels[i] +"+H");

                t.mzs.push_back( t.mzs[i] - PROTON );
                t.intensity_array.push_back( t.intensity_array[i] );
                t.fragment_labels.push_back(t.fragment_labels[i]+"-H");
            }
        }

        //theory fragmentation or library fragmentation = t
        //experimental data = f

//        cerr << endl;
//        cerr << "Compound: " << cpd->name << endl;
//        for (unsigned int i = 0; i < t.mzs.size(); i++) {
//            float mz = t.mzs.at(i);
//            float intensity = t.intensity_array.at(i);
//            string label = t.fragment_labels.at(i);
//            cerr << "mz=" << mz << ", intensity= "<< intensity << ", label=" << label << endl;
//        }
//        cerr << endl;

//        cerr << "Observed (precursorMz=" << f->precursorMz << ")" << endl;
//        for (float mz : f->mzs) {
//            cerr << mz << endl;
//        }
//        cerr << endl;

        s = t.scoreMatch(f,productPpmTolr);

//        cerr << endl;
//        cerr << "t after scoring match " << cpd->name << endl;
//        for (unsigned int i = 0; i < t.mzs.size(); i++) {
//            float mz = t.mzs.at(i);
//            float intensity = t.intensity_array.at(i);
//            string label = t.fragment_labels.at(i);
//            cerr << "mz=" << mz << ", intensity= " << intensity << ", label=" << label << endl;
//        }
//        cerr << endl;

        return s;
}

vector<Compound*> Compound::getChildren() {return vector<Compound*>(0);}
vector<Compound*> SummarizedCompound::getChildren() {return children;}

/**
 * @brief SummarizedCompound::computeSummarizedData
 *
 * fragment information:
 * All m/z values from all compounds are used.
 * Intensity values are averaged based on all m/zs with some intensity value.
 *
 * metadata information:
 * retain all metdata associated with every compound.
 */
void SummarizedCompound::computeSummarizedData() {

    map<int, vector<float>> intensitiesByMz = {};
    map<int, vector<string>> labelsByMz = {};

    map<pair<string, string>, int> summarizedMetaDataMap{};

    for (auto compound : getChildren()) {
        for (unsigned int i = 0; i < compound->fragment_mzs.size(); i++) {

            int mzKey = mzUtils::mzToIntKey(compound->fragment_mzs[i], 1000000);
            float intensity = compound->fragment_intensity[i];
            string label = compound->fragment_labels[i];

            if (intensitiesByMz.find(mzKey) == intensitiesByMz.end()){
                vector<float> intensities(1);
                intensities[0] = intensity;
                intensitiesByMz.insert(make_pair(mzKey, intensities));
            } else {
                intensitiesByMz[mzKey].push_back(intensity);
            }

            if (labelsByMz.find(mzKey) == labelsByMz.end()){
                vector<string> labels(1);
                labels[0] = label;
                labelsByMz.insert(make_pair(mzKey, labels));
            } else {
                labelsByMz[mzKey].push_back(label);
            }
        }
    }

    fragment_mzs = vector<float>(intensitiesByMz.size());
    fragment_intensity=vector<float>(intensitiesByMz.size());
    fragment_labels = vector<string>(intensitiesByMz.size());

    unsigned int vecCounter = 0;
    for (map<int, vector<float>>::iterator it = intensitiesByMz.begin(); it != intensitiesByMz.end(); ++it) {

        int mzInt = it->first;
        vector<float> intensities = it->second;

        float avgInt = 0.0f;
        for (auto intensity : intensities) {
            avgInt += intensity;
        }
        avgInt /= intensities.size();

        fragment_mzs[vecCounter] = mzUtils::intKeyToMz(mzInt, 1000000);
        fragment_intensity[vecCounter] = avgInt;

        vector<string> labels = labelsByMz[mzInt];
        sort(labels.begin(), labels.end());
        if (labels.size() > 0) {

            string lastLabel = labels[0];
            string mergedLabel = lastLabel;

            if (labels.size() > 1) {
                for (unsigned int i = 1; i < labels.size(); i++) {

                    string thisLabel = labels[i];

                    if (lastLabel == thisLabel) continue;

                    mergedLabel = mergedLabel + "; " + thisLabel;

                    lastLabel = thisLabel;
                }
                fragment_labels[vecCounter] = mergedLabel;
            } else {
                fragment_labels[vecCounter] = lastLabel;
            }
        } else {
            fragment_labels[vecCounter] = "";
        }

        vecCounter++;
    }
}
