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

//breadth-first search
vector<Compound*> SummarizedCompound::getChildren() {
    //output
    vector<Compound*> descendants;

    //becomes the next childset (in the next iteration)
    vector<Compound*> nextGeneration;

    //Initial condition
    vector<Compound*> thisGeneration = children;

    //continue iterating until no more children to retrieve
    while (thisGeneration.size() > 0) {

        for (auto compound : thisGeneration) {
            if (compound->db == "summarized") {
                vector<Compound*> compoundChildren = compound->getChildren();
                nextGeneration.insert(nextGeneration.end(), compoundChildren.begin(), compoundChildren.end());
            } else {
                descendants.push_back(compound);
            }
        }

        thisGeneration = nextGeneration;
        nextGeneration.clear();
    }

    return descendants;
}

/**
 * @brief SummarizedCompound::computeSummarizedData
 *
 * fragment information:
 * All m/z values from all compounds are used.
 * Intensity values are averaged when an m/z is found in more than one compound.
 * Labels should be combined for all summarized compounds.
 *
 * metadata information:
 * retain all metdata associated with every compound.
 */
void SummarizedCompound::computeSummarizedData() {

    map<int, vector<float>> intensitiesByMz = {};
    map<int, vector<string>> labelsByMz = {};

    map<pair<string, string>, unsigned long> summarizedMetaDataMap{};

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

        //Issue 216: metadata
        for (auto it = compound->metaDataMap.begin(); it != compound->metaDataMap.end(); ++it){
            pair<string, string> metaDataPair = pair<string, string>(it->first, it->second);
            if (summarizedMetaDataMap.find(metaDataPair) == summarizedMetaDataMap.end()){
                summarizedMetaDataMap.insert(make_pair(metaDataPair, 0));
            }
            summarizedMetaDataMap[metaDataPair]++;
        }
    }

    //Issue 216: metadata
    for (auto it = summarizedMetaDataMap.begin(); it != summarizedMetaDataMap.end(); ++it){
        pair<string, string> metaDataPair = it->first;
        if (it->second == getChildren().size()) {
            this->metaDataMap.insert(metaDataPair);
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

                    mergedLabel = mergedLabel + "/" + thisLabel;

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

/**
 * @brief Ms3Compound::computeMs3Spectra
 *
 * relies on baseCompound, specifically that baseCompound has labels appropriately.
 *
 * eg
 * 325.243 76 ms3-{605.2}-NL Acyl Chain 1
 *
 * need to start with "ms3-"
 * then the ms2 precursor m/z in {brackets}
 * followed by the rest of the name, after a hyphen
 */
void Ms3Compound::computeMs3Spectra() {
    for (unsigned int i = 0; i < baseCompound->fragment_labels.size(); i++){

        string label = baseCompound->fragment_labels[i];
        float fragMz = baseCompound->fragment_mzs[i];
        float fragIntensity = baseCompound->fragment_intensity[i];

        if (label.find("ms3-{") == 0) {
            unsigned long firstCloseBracket = label.find("}");

            if (firstCloseBracket != string::npos && firstCloseBracket > 5 && firstCloseBracket < label.size()-1){

                string ms3precMzStr = label.substr(5,firstCloseBracket-6);
                string ms3precLabel = label.substr(firstCloseBracket+2, label.size());

                double ms3precMz = -1.0;
                try {
                    ms3precMz = stod(ms3precMzStr);
                } catch (...) {}

                if (ms3precMz > -1.0){

                    int mzKey = mzToIntKey(ms3precMz);

                    if (ms3_fragment_mzs.find(mzKey) == ms3_fragment_mzs.end()) {
                        ms3_fragment_mzs.insert(make_pair(mzKey, vector<float>()));
                    }
                    if (ms3_fragment_intensity.find(mzKey) == ms3_fragment_intensity.end()){
                        ms3_fragment_intensity.insert(make_pair(mzKey, vector<float>()));
                    }
                    if (ms3_fragment_labels.find(mzKey) == ms3_fragment_labels.end()){
                        ms3_fragment_labels.insert(make_pair(mzKey, vector<string>()));
                    }

                    ms3_fragment_mzs[mzKey].push_back(fragMz);
                    ms3_fragment_intensity[mzKey].push_back(fragIntensity);
                    ms3_fragment_labels[mzKey].push_back(ms3precLabel);
                }

            }
        } else {
            fragment_labels.push_back(label);
            fragment_mzs.push_back(fragMz);
            fragment_intensity.push_back(fragIntensity);
        }
    }
}

vector<Compound*> Ms3Compound::getChildren(){
    return vector<Compound*>{baseCompound};
}

vector<pair<string, string>> SummarizedCompound::parseCompoundId(string compoundId, bool debug){
    auto adductPart = compoundId.find("}[");

    string compoundIdNoAdduct = compoundId.substr(0, adductPart+1);

    if (debug) cout << "compound Id no adduct: " << compoundIdNoAdduct << endl;

    auto compositionSummaryStart = compoundId.find("}={");

    vector<pair<string, string>> compoundAdductPairs{};

    if (compositionSummaryStart == string::npos) {

        //general
        string generalSummaryCompoundList = compoundIdNoAdduct.substr(1, compoundIdNoAdduct.size()-2) + ";";

        if (debug) cout << "general summary compound list: " << generalSummaryCompoundList << endl;

        vector<string> encodedCompoundAdductPairs{};

        unsigned long posPrevious = 0;
        unsigned long posCurrent = 0;

        while ((posCurrent = generalSummaryCompoundList.find(";", posPrevious)) != string::npos) {
            string encodedCompoundAdductPair = generalSummaryCompoundList.substr(posPrevious, posCurrent-posPrevious);
            posPrevious = posCurrent + 1;

            encodedCompoundAdductPairs.push_back(encodedCompoundAdductPair);

            if (debug) cout << "encoded compound-adduct pair: " << encodedCompoundAdductPair << endl;
        }

        for (auto encodedPair : encodedCompoundAdductPairs) {
            auto splitPos = encodedPair.find("|");
            string compound = encodedPair.substr(0, splitPos);
            string adduct = encodedPair.substr(splitPos+1, encodedPair.size());

            if (debug) cout << "encoded compound: " << compound << ", encoded adduct: " << adduct << endl;

            compoundAdductPairs.push_back(make_pair(compound, adduct));
        }

    } else {
        //summarized
        string compositionSummaryAdduct = compoundId.substr(adductPart+1);
        if (debug) cout << "composition summary adduct: " << compositionSummaryAdduct << endl;

        string compositionSummaryCompoundList = compoundIdNoAdduct.substr(compositionSummaryStart+3);
        compositionSummaryCompoundList = compositionSummaryCompoundList.substr(0, compositionSummaryCompoundList.size()-1) + ";";

        if (debug) cout << "composition summary compound list: " << compositionSummaryCompoundList << endl;

        unsigned long posPrevious = 0;
        unsigned long posCurrent = 0;

        while ((posCurrent = compositionSummaryCompoundList.find(";", posPrevious)) != string::npos) {
            string compoundName = compositionSummaryCompoundList.substr(posPrevious, posCurrent-posPrevious);
            posPrevious = posCurrent + 1;

            compoundAdductPairs.push_back(make_pair(compoundName, compositionSummaryAdduct));
        }
    }

    return compoundAdductPairs;
}
