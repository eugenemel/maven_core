#include "mzSample.h"
#include <regex>

//TODO: test this
//TODO: this should replace MainWindow::getSrmSlices()
pair<vector<mzSlice*>, vector<SRMTransition*>> QQQProcessor::getSRMSlices(
        vector<mzSample*> samples,
        shared_ptr<QQQSearchParameters> params,
        vector<Compound*> compounds,
        vector<Adduct*> adducts,
        bool debug){

    set<string> srms;
    //+118.001@cid34.00 [57.500-58.500]
    //+ c ESI SRM ms2 102.000@cid19.00 [57.500-58.500]
    //-87.000 [42.500-43.500]
    //- c ESI SRM ms2 159.000 [113.500-114.500]

    regex rx1a("[+/-](\\d+\\.\\d+)");
    regex rx1b("ms2\\s*(\\d+\\.\\d+)");
    regex rx2("(\\d+\\.\\d+)-(\\d+\\.\\d+)");
    int countMatches=0;

    float amuQ1 = params->amuQ1;
    float amuQ3 = params->amuQ3;

    vector<mzSlice*>slices;
    map<pair<float, float>, SRMTransition*> srmTransitions{};

    for(unsigned int i=0; i < samples.size(); i++ ) {
        mzSample* sample = samples[i];
        for(unsigned int j=0; j < sample->scans.size(); j++ ) {
            Scan* scan = sample->getScan(j);
            if (!scan) continue;

            string filterLine(scan->filterLine.c_str());
            if (filterLine.empty()) continue;

            if (srms.find(filterLine) != srms.end()) continue;
            //if (srms.contains(filterLine))  continue;
            srms.insert(filterLine);

            mzSlice* s = new mzSlice(0,0,0,0);
            s->srmId = scan->filterLine.c_str();
            slices.push_back(s);

            //match compounds
            Compound* compound = nullptr;
            Adduct* adduct = nullptr;

            float precursorMz = scan->precursorMz;
            float productMz   = scan->productMz;
            int   polarity= scan->getPolarity();
            if (polarity==0) filterLine[0] == '+' ? polarity=1 : polarity =-1;

            if ( precursorMz < 1e-6f ) {

                smatch m1, m2;
                regex_search(filterLine, m1, rx1a);
                regex_search(filterLine, m2, rx1b);

                if (m1.size() > 0) {
                    precursorMz = stof(m1[0]);
                } else if (m2.size() > 0){
                    precursorMz = stof(m2[0]);
                }
            }

            if ( productMz < 1e-6f ) {

                smatch m3;
                regex_search(filterLine, m3, rx2);

                if (m3.size() > 2) {
                    float lb = stof(m3[1]);
                    float ub = stof(m3[2]);
                    productMz = lb+(0.5f*(ub-lb));
                }
            }

            //relies on compounds
            if (precursorMz > 0 && productMz > 0 ) {
                float dist = FLT_MAX;
                for (auto db_compound : compounds) {

                    //only consider SRM compounds
                    if (db_compound->precursorMz <= 0 || db_compound->productMz <= 0) continue;

                    //compounds must be within tolerance
                    auto q1_dist = abs(db_compound->precursorMz-precursorMz) ;
                    if (q1_dist > amuQ1) continue;

                    auto q3_dist = abs(db_compound->productMz-productMz);
                    if (q3_dist > amuQ3) continue;

                    auto d = sqrt(q1_dist*q1_dist + q3_dist*q3_dist);
                    if (d < dist) {
                        compound = db_compound;
                        dist = d;
                    }
                }
            }

            if (compound) {
                compound->srmId=filterLine;
                s->compound=compound;
                s->rt = compound->expectedRt;
                countMatches++;

                if (!compound->adductString.empty()) {
                    for (auto availableAdduct : adducts) {
                        if (availableAdduct->name == compound->adductString) {
                            adduct = availableAdduct;
                            break;
                        }
                    }
                }
            }

            //Issue 347
            pair<float, float> srmKey = make_pair(precursorMz, productMz);
            SRMTransition *srmTransition = new SRMTransition();
            if (srmTransitions.find(srmKey) != srmTransitions.end()) {
                srmTransition = srmTransitions[srmKey];
            } else {
                srmTransition->precursorMz = precursorMz;
                srmTransition->productMz = productMz;
            }

            srmTransition->mzSlices.push_back(make_pair(sample, s));

            if (compound){
                srmTransition->compound = compound;
                if (compound->expectedRt > 0) {
                    srmTransition->rt = compound->expectedRt;
                }
            }
            if (adduct) srmTransition->adduct = adduct;

            srmTransitions[srmKey] = srmTransition;
        }
        //qDebug() << "SRM mapping: " << countMatches << " compounds mapped out of " << srms.size();
    }

    vector<SRMTransition*> srmTransitionsAsVector(srmTransitions.size());

    unsigned int counter = 0;
    for (auto it = srmTransitions.begin(); it != srmTransitions.end(); ++it) {
        srmTransitionsAsVector[counter] = it->second;
        counter++;
    }

    return make_pair(slices, srmTransitionsAsVector);
}
