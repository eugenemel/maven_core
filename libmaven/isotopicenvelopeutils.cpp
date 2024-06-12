#include "isotopicenvelopeutils.h"
#include "mzMassCalculator.h"

double IsotopicEnvelope::getTotalIntensity() {
    if (totalIntensity < 0) {
        totalIntensity = std::accumulate(intensities.begin(), intensities.end(), 0.0);
    }
    return totalIntensity;
}

void IsotopicEnvelopeGroup::print() {
    if (!group) {
        cout << "Peak group is null." << endl;
        return;
    }

    if (!group->compound) {
        cout << "Compound is null." << endl;
        return;
    }

    if (!group->adduct) {
        cout << "Adduct is null." << endl;
        return;
    }

    cout << group->compound->name << " "
         << group->adduct->name << " (m/z, rt) = ("
         << group->compound->precursorMz << ", "
         << group->medianRt()
         << ")"
         << endl;

    cout << "Isotopes: ";
    cout << "{";
    for (unsigned int i = 0; i < isotopes.size(); i++) {
        Isotope isotope = isotopes.at(i);
        if (i > 0) cout << ", ";
        cout << isotope.name;
    }
    cout << "}\n";

    for (auto it = envelopeBySample.begin(); it != envelopeBySample.end(); ++it) {
        mzSample *sample = it->first;
        IsotopicEnvelope envelope = it->second;
        cout << sample->sampleName << ": ";
        envelope.print();
        cout << "\n";
    }

    cout << endl;
}

IsotopicEnvelope& IsotopicEnvelope::none() {
    static IsotopicEnvelope none;
    none.source = "NONE";
    return none;
}

void IsotopicEnvelope::print() {

    stringstream ss;
    ss << std::fixed << setprecision(4);

    ss << "{";
    for (unsigned int i = 0; i < intensities.size(); i++) {
        if (i > 0) ss << ", ";
        ss << intensities.at(i);
    }
    ss << "}; Total=";
    ss << getTotalIntensity();

    cout << ss.str();
}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopes(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<mzSample*>& samples,
    IsotopeParameters params,
    bool debug){

    //Issue 671: Respect isotopes parameter
    int maxNumProtons = INT_MAX;
    if (params.isExtractNIsotopes) {
        maxNumProtons = params.maxIsotopesToExtract;
    }

    string compoundFormula = "";
    if (compound) {
        compoundFormula = compound->formula;
    }

    vector<Isotope> theoreticalIsotopes = MassCalculator::computeIsotopes(
        compoundFormula,
        adduct,
        static_cast<double>(group->meanMz),
        params.getLabeledIsotopes(),
        params.labeledIsotopeRetentionPolicy,
        NaturalAbundanceData::defaultNaturalAbundanceData,
        params.isNatAbundance,
        maxNumProtons,
        params.natAbundanceThreshold,
        debug
        );

    //TODO: currently does not consider possibility of unknown atomic composition isotopes
    vector<Isotope> isotopes = theoreticalIsotopes;
    if (params.isCondenseTheoreticalIsotopes) {
        isotopes = IsotopicEnvelopeAdjuster::condenseTheoreticalIsotopes(
            theoreticalIsotopes,
            params,
            debug
            );
    }

    IsotopicEnvelopeGroup envelopeGroup;

    if (params.isotopicExtractionAlgorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA) {
        envelopeGroup = extractEnvelopesPeakFullRtBounds(compound, adduct, group, isotopes, params, debug);
    } else if (params.isotopicExtractionAlgorithm == IsotopicExtractionAlgorithm::MAVEN_GUI_VERSION_ONE) {
        envelopeGroup = extractEnvelopesVersion1(compound, adduct, group, isotopes, params, debug);
    }

    //Issue 691: combine overlapping isotopes when the peak height for nearby isotopes is identical.
    if (params.isCombineOverlappingIsotopes) {
        envelopeGroup.combineOverlappingIsotopes(params.ppm, debug);
    }

    envelopeGroup.extractionAlgorithmName = IsotopeParameters::getAlgorithmName(params.isotopicExtractionAlgorithm);
    return envelopeGroup;

}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopesPeakFullRtBounds(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<Isotope>& isotopes,
    IsotopeParameters params,
    bool debug) {

    IsotopicEnvelopeGroup envelopeGroup;

    envelopeGroup.compound = compound;
    envelopeGroup.adduct = adduct;
    envelopeGroup.group = group;
    envelopeGroup.isotopes = isotopes;

    //initialize peak groups
    envelopeGroup.isotopePeakGroups = vector<PeakGroup>(isotopes.size());
    for (unsigned int i = 0; i < envelopeGroup.isotopePeakGroups.size(); i++) {
        Isotope isotope = isotopes[i];

        PeakGroup g;
        g.meanMz = isotope.mz;
        g.tagString = isotope.name;
        g.expectedAbundance = isotope.abundance;
        g.isotopeC13count= isotope.C13;
        g.isotopicIndex = i;
        g.compound = compound;
        g.adduct = adduct;
        g.setType(PeakGroup::GroupType::IsotopeType);

        envelopeGroup.isotopePeakGroups[i] = g;
    }

    for (auto & peak : group->peaks) {

        IsotopicEnvelope envelope;
        envelope.intensities = vector<double>(isotopes.size());

        auto sample = peak.sample;

        for (unsigned int i = 0; i < isotopes.size(); i++) {

            Isotope isotope = isotopes.at(i);

            EIC *eic = sample->getEIC(
                static_cast<float>(isotope.mz - isotope.mz/1e6f*params.ppm),
                static_cast<float>(isotope.mz + isotope.mz/1e6f*params.ppm),
                peak.rtmin,
                peak.rtmax,
                1);

            float mzmin = eic->mzmin;
            float mzmax = eic->mzmax;

            double intensity = std::accumulate(eic->intensity.begin(), eic->intensity.end(), 0.0);

            //Issue 695: Collate with raw data for proper rendering in MAVEN GUI.
            float maxIntensity = 0.0f;
            float rtAtMaxIntensity = peak.rt;
            for (unsigned int i = 0; i < eic->size(); i++) {
                if (eic->intensity[i] > maxIntensity) {
                    maxIntensity = eic->intensity[i];
                    rtAtMaxIntensity = eic->rt[i];
                }
            }

            delete(eic);
            eic = nullptr;

            envelope.intensities.at(i) = intensity;

            Peak p;

            p.rt = rtAtMaxIntensity;
            p.peakIntensity = maxIntensity;

            //RT related
            p.scan = peak.scan;
            p.minscan = peak.minscan;
            p.maxscan = peak.maxscan;
            p.pos = peak.scan;
            p.minpos = peak.minpos;
            p.maxpos = peak.maxpos;
            p.width = (peak.maxpos-peak.minpos);
            p.rtmin = peak.rtmin;
            p.rtmax = peak.rtmax;

            //m/z related
            p.mzmin = mzmin;
            p.mzmax = mzmax;
            p.peakMz = isotope.mz;
            p.baseMz = isotope.mz;
            p.medianMz = isotope.mz;

            //intensity related
            //Issue 680: for easier review in MAVEN, just set all quant types to same value.
            p.peakArea = intensity;
            p.peakAreaCorrected = intensity;
            p.peakAreaTop = intensity;
            p.peakAreaFractional = intensity;
            p.signalBaselineRatio = intensity;
            p.smoothedIntensity = intensity;
            p.smoothedPeakArea = intensity;
            p.smoothedPeakAreaCorrected = intensity;
            p.smoothedPeakAreaTop = intensity;
            p.smoothedSignalBaselineRatio = intensity;
            p.peakAreaFWHM = intensity;
            p.smoothedPeakAreaFWHM = intensity;

            //avoid writing junk into mzrollDB
            p.gaussFitR2 = 0;
            p.noNoiseObs = 0;
            p.signalBaselineRatio = 0;

            p.sample = sample;

            envelopeGroup.isotopePeakGroups[i].addPeak(p);
        }

        envelope.getTotalIntensity();

        envelopeGroup.envelopeBySample.insert(make_pair(sample, envelope));
    }

    return envelopeGroup;
}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopesVersion1(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        IsotopeParameters params,
        bool debug) {

    IsotopicEnvelopeGroup envelopeGroup;

    envelopeGroup.compound = compound;
    envelopeGroup.adduct = adduct;
    envelopeGroup.group = group;
    envelopeGroup.isotopes = isotopes;

    vector<PeakGroup> candidateIsotopePeakGroups = vector<PeakGroup>(isotopes.size());

    for (unsigned int i = 0; i < isotopes.size(); i++) {
        Isotope isotope = isotopes[i];

        PeakGroup g;
        g.meanMz = isotope.mz;
        g.tagString = isotope.name;
        g.expectedAbundance = isotope.abundance;
        g.isotopeC13count= isotope.C13;
        g.isotopicIndex = i;
        g.compound = compound;
        g.adduct = adduct;
        g.setType(PeakGroup::GroupType::IsotopeType);

        candidateIsotopePeakGroups[i] = g;
    }

    for (auto parentPeak : group->peaks ) {

        IsotopicEnvelope envelope;
        envelope.intensities = vector<double>(isotopes.size());

        auto sample = parentPeak.sample;

        if (debug) {
            cout << "\n" << compound->name << " " << adduct->name << endl;
        }

        for (unsigned int i = 0; i < isotopes.size(); i++) {

            Isotope isotope = isotopes.at(i);

            if (debug) {
                cout << "Starting " << isotope.name << " " << sample->sampleName << "..." << endl;
            }

            float mzmin = isotope.mz-isotope.mz/1e6f*params.ppm;
            float mzmax = isotope.mz+isotope.mz/1e6f*params.ppm;

            float rt  =   parentPeak.rt;
            float rtmin = parentPeak.rtmin;
            float rtmax = parentPeak.rtmax;

            float parentPeakIntensity= parentPeak.peakIntensity;

            int scannum = parentPeak.getScan()->scannum;

            //TODO: not hard-coded number of scans value
            int minScan = max(0, (scannum - 3));
            int maxScan = min(static_cast<int>(parentPeak.sample->scanCount()), (scannum + 3));

            float isotopePeakIntensity=0;

            if (debug) {
                cout << "Searching for scans in range [" << minScan << " - " << maxScan << "] ..." << endl;
            }

            //Refine RT of isotopic peak, determine corresponding intensity.
            for (int i= minScan; i < maxScan; i++) {

                Scan* s = sample->getScan(static_cast<unsigned int>(i));

                //look for isotopic mass in the same spectrum
                vector<int> matches = s->findMatchingMzs(mzmin, mzmax);

                for(unsigned int j=0; j < matches.size(); j++ ) {
                    int pos = matches[j];
                    if (s->intensity[static_cast<unsigned int>(pos)] > isotopePeakIntensity ) {
                        isotopePeakIntensity = s->intensity[static_cast<unsigned int>(pos)];
                        rt = s->rt;
                        if (debug) {
                            cout << compound->name << " " << adduct->name << " " << isotope.name << " " << sample->sampleName
                                 << " Scan #" << s->scannum << ": m/z=" << s->mz.at(pos) << ", intensity=" << isotopePeakIntensity
                                 << endl;
                    }
                    }
                }
            }

            if (debug && isotopePeakIntensity < 1e-6f) {
                cout << compound->name << " " << adduct->name << " " << isotope.name << " " << sample->sampleName
                     << " isotope not detected. "<< endl;
            }

            //natural abundance check
            if (params.isIgnoreNaturalAbundance) {

                //peaks are too low intensity to observe? TODO: consider removing / restructuring
                if (isotope.abundance < 1e-8) continue;
                if (isotope.abundance * parentPeakIntensity < 1) continue;

                float observedAbundance = isotopePeakIntensity/(parentPeakIntensity+isotopePeakIntensity);
                float naturalAbundanceError = abs(observedAbundance-isotope.abundance)/isotope.abundance*100.0f;

                if (debug) {
                    cout << sample->sampleName << ", "
                         << isotope.name << ": "
                         << "Expected abundance= " << isotope.abundance << " "
                         << "Observed abundance= " << observedAbundance << " "
                         << "Error= "     << naturalAbundanceError << endl;
                }

                //C12 parent peak is always OK, cannot be disqualified on account of observed/expected abundance issues
                if (naturalAbundanceError > static_cast<float>(params.maxNaturalAbundanceErr) && !isotope.isParent())  continue;
            }

            float w = static_cast<float>(params.maxIsotopeScanDiff)*sample->getAverageFullScanTime();

            //Issue 120: Use sample-specific mz value for peaks instead of average mz
            double c = static_cast<double>(
                sample->correlation(static_cast<float>(isotope.mz), parentPeak.peakMz, params.ppm, rtmin-w,rtmax+w, false)
                );

            if (c < params.minIsotopicCorrelation) continue;

            if (debug) {
                cout << compound->name << " " << adduct->name << " " << isotope.name << " " << sample->sampleName
                     << isotope.mz << " [" << (rtmin-w) << " - " <<  (rtmax+w) << "] correlation=" << c << endl;
            }

            EIC* eic=nullptr;

            for( int i=0; i <= params.maxIsotopeScanDiff; i++ ) {
                float window=i*sample->getAverageFullScanTime();

                //TODO: Issue 371: Handle SRMTransitionType isotopes
                eic = sample->getEIC(mzmin,mzmax,rtmin-window,rtmax+window,1);

                if(!eic) continue;

                eic->getPeakPositionsD(params.peakPickingAndGroupingParameters, debug);

                //TODO: is this a good metric for stopping isotope extraction?
                if (eic->peaks.size() >= 1 ) break;

                //clean up
                delete(eic);
                eic=nullptr;
            }

            if (eic) {
                Peak* nearestPeak=nullptr;
                float d=FLT_MAX;
                for(unsigned int i=0; i < eic->peaks.size(); i++ ) {
                    Peak& x = eic->peaks[i];
                    float dist = abs(x.rt - rt);
                    if ( dist > static_cast<float>(params.maxIsotopeScanDiff)*sample->getAverageFullScanTime()) continue;
                    if ( dist < d ) { d=dist; nearestPeak = &x; }
                }

                if (nearestPeak) {
                    candidateIsotopePeakGroups[i].addPeak(*nearestPeak);

                    //quant type used for natural abundance correction, so used here
                    envelope.intensities.at(i) = nearestPeak->peakIntensity;
                }
                delete(eic);
                eic = nullptr;
            }
        }

        envelope.getTotalIntensity();

        envelopeGroup.envelopeBySample.insert(make_pair(sample, envelope));
    }

    //Issue 615: respect option to retain/exclude peak groups that have no peaks
    if (!params.isKeepEmptyIsotopes) {
        for (auto pg : candidateIsotopePeakGroups) {
            if (pg.peakCount() > 0) {
                envelopeGroup.isotopePeakGroups.push_back(pg);
            }
        }
    } else {
        envelopeGroup.isotopePeakGroups = candidateIsotopePeakGroups;
    }

    return envelopeGroup;
}

IsotopicEnvelopeGroup IsotopicEnvelopeExtractor::extractEnvelopesPeakShrinkingRtBounds(
    Compound *compound,
    Adduct *adduct,
    PeakGroup *group,
    vector<Isotope>& isotopes,
    IsotopeParameters params,
    bool debug){

    IsotopicEnvelopeGroup envelopeGroup;

    //TODO

    return envelopeGroup;
}

vector<Isotope> IsotopicEnvelopeAdjuster::condenseTheoreticalIsotopes(
    vector<Isotope> defaultIsotopes,
    IsotopeParameters params,
    bool debug){

    //WARNING:
    //Isotope.abundance stores the predicted theoretical abundance.
    //In cases where isotopes should be combined, but only one atom is labeled,
    //the natural abundance of the unlabeled atoms must be accounted for.
    //For example, if 13C is labeled, but the 13C cannot be resolved from the 15N,
    //the [M+1] relative abundance must take into account both 13C and 15N.

    //TODO: implement
    return defaultIsotopes;
}

void IsotopicEnvelopeGroup::setIsotopesToChildrenPeakGroups(Classifier *classifier){

    this->group->children.clear();

    for (auto & child : isotopePeakGroups) {

        child.metaGroupId = this->group->metaGroupId;
        child.compound = this->compound;
        child.adduct= this->adduct;
        child.parent = this->group;
        child.setType(PeakGroup::IsotopeType);

        child.groupStatistics();
        if (classifier && classifier->hasModel()) {
            classifier->classify(&child);
            child.groupStatistics();
        }
        this->group->addChild(child);
    }
}

void IsotopicEnvelopeGroup::combineOverlappingIsotopes(float ppm, bool debug) {

    //organize isotopes by name
    map<string, Isotope> isotopesByName{};
    for (Isotope isotope : isotopes) {
        isotopesByName.insert(make_pair(isotope.name, isotope));
    }

    //re-organize data, noting overlapping quant values
    map<pair<mzSample*, float>, vector<string>> quantValToIsotopes{};

    //define new destinations for overlapping quant values
    map<Peak, vector<string>, PeakIntensitySampleComparator> peakToUpdatedPeakGroup{};

    //keep mapping betwen duplicated keys, in case we need those peaks again
    map<Peak, vector<Peak>, PeakIntensitySampleComparator> peakToAllPeaks{};

    PeakGroup monoisotope;

    for (PeakGroup child : isotopePeakGroups) {

        string isotopeName = child.tagString;

        if (isotopeName == "C12 PARENT") {
            monoisotope = child;
            continue;
        }

        for (Peak peak : child.getPeaks()) {
            float peakHeight = peak.peakIntensity;
            mzSample *sample = peak.getSample();

            if (sample && peakHeight > 0) {
                auto key = make_pair(sample, peakHeight);
                if (quantValToIsotopes.find(key) == quantValToIsotopes.end()) {
                    quantValToIsotopes.insert(make_pair(key, vector<string>{}));
                }
                quantValToIsotopes[key].push_back(isotopeName);

                //If the peaks cannot be merged together after all, then need to send peaks back to original isotope.
                peak.tempString = isotopeName;

                if (peakToUpdatedPeakGroup.find(peak) == peakToUpdatedPeakGroup.end()) {
                    peakToUpdatedPeakGroup.insert(make_pair(peak, vector<string>{}));
                    peakToAllPeaks.insert(make_pair(peak, vector<Peak>{}));
                }
                peakToUpdatedPeakGroup[peak].push_back(isotopeName);
                peakToAllPeaks[peak].push_back(peak);

            }
        }
    }


    if (debug) {
        cout << "[IsotopicEnvelopeGroup::combineOverlappingIsotopes()] quantValToIsotopes:" << endl;
        for (auto it = quantValToIsotopes.begin(); it != quantValToIsotopes.end(); ++it) {
            cout << it->first.first->sampleName << ", height=" << it->first.second << ": {";
            for (unsigned int i = 0; i < it->second.size(); i++) {
                if (i > 0) cout << ",";
                cout << it->second[i];
            }
            cout << "}\n";
        }
        cout << endl;
    }

    //identify all of the cases where the same quant value is identified by multiple isotopes.
    set<vector<string>> combinations{};
    for (auto it = quantValToIsotopes.begin(); it != quantValToIsotopes.end(); ++it) {
        if (it->second.size() > 1) {
            combinations.insert(it->second);
        }
    }

    if (debug) {
        cout << "[IsotopicEnvelopeGroup::combineOverlappingIsotopes()] combinations:" << endl;
        for (auto it = combinations.begin(); it != combinations.end(); ++it) {
            cout << "{";
            vector<string> combos = (*it);
            for (unsigned int i = 0; i < combos.size(); i++) {
                if (i > 0) cout << ", ";
                cout << combos[i];
            }
            cout << "}\n";
        }
        cout << endl;
    }

    map<string, bool> isIsotopesOverlappingMz{};

    //add new isotopes based on combinations (if applicable).
    for (auto it = combinations.begin(); it != combinations.end(); ++it){
        Isotope combinedIsotope;

        vector<double> isotopeMzs{};

        for (auto & isotopeName : *it) {
            Isotope isotope = isotopesByName[isotopeName];

            if (combinedIsotope.name != "") {
                combinedIsotope.name += " / ";
                combinedIsotope.name += isotope.name;
            } else {
                combinedIsotope.name = isotope.name;
                combinedIsotope.charge = isotope.charge;
            }

            isotopeMzs.push_back(isotope.mz);

            combinedIsotope.mz += isotope.mz;
            combinedIsotope.abundance += isotope.abundance;
            combinedIsotope.naturalAbundanceMonoProportion += isotope.naturalAbundanceMonoProportion;
            if (isotope.C13) combinedIsotope.C13 = true;
            if (isotope.N15) combinedIsotope.N15 = true;
            if (isotope.H2) combinedIsotope.H2 = true;
            if (isotope.S34) combinedIsotope.S34 = true;
            if (isotope.O18) combinedIsotope.O18 = true;

        }

        sort(isotopeMzs.begin(), isotopeMzs.end());

        bool isMzsOverlapping = true;

        for (unsigned int i = 1; i < isotopeMzs.size(); i++) {

            double previousMz = isotopeMzs[i-1];
            double currentMz = isotopeMzs[i];

            float maxLowerMz = previousMz + previousMz/1e6f*ppm;
            float minHigherMz = currentMz - currentMz/1e6f*ppm;

            isMzsOverlapping = isMzsOverlapping && maxLowerMz >= minHigherMz;
        }

        isIsotopesOverlappingMz.insert(make_pair(combinedIsotope.name, isMzsOverlapping));

        isotopesByName.insert(make_pair(combinedIsotope.name, combinedIsotope));
    }

    //Reformat data, mapping to peaks to a single string key
    map<string, vector<Peak>> isotopeToPeaks{};
    for (auto it = peakToUpdatedPeakGroup.begin(); it != peakToUpdatedPeakGroup.end(); ++it) {

        string isotopeNameKey;

        vector<string> isotopeNames = it->second;
        for (string isotopeName : isotopeNames) {
            if (isotopeNameKey != "") {
                isotopeNameKey += " / ";
                isotopeNameKey += isotopeName;
            } else {
                isotopeNameKey = isotopeName;
            }
        }

        Peak peak = it->first;

        bool isValidMerge = true;
        if (isIsotopesOverlappingMz.find(isotopeNameKey) != isIsotopesOverlappingMz.end()) {
            isValidMerge = isIsotopesOverlappingMz.at(isotopeNameKey);
        }

        if (isValidMerge) {
            if (isotopeToPeaks.find(isotopeNameKey) == isotopeToPeaks.end()) {
                isotopeToPeaks.insert(make_pair(isotopeNameKey, vector<Peak>{}));
            }
            isotopeToPeaks[isotopeNameKey].push_back(peak);
        } else {
            vector<Peak> originalPeaks = peakToAllPeaks.at(peak);
            for(Peak p : originalPeaks) {
                string originalIsotopeName = p.tempString;
                if (isotopeToPeaks.find(originalIsotopeName) == isotopeToPeaks.end()) {
                    isotopeToPeaks.insert(make_pair(originalIsotopeName, vector<Peak>{}));
                }
                isotopeToPeaks[originalIsotopeName].push_back(p);
            }
        }

    }

    // Case: peak would have been mapped to a merged isotope, but not all of the original isotopes were detected
    // TODO

    unsigned int counter = 0;

    //create new peak groups based on merged isotopes
    vector<PeakGroup> updatedPeakGroups(isotopeToPeaks.size()+1); //add one for C12 PARENT
    vector<Isotope> updatedIsotopes(isotopeToPeaks.size()+1);

    //add back C12 parent
    if (isotopesByName.find("C12 PARENT") != isotopesByName.end()) {
        updatedIsotopes[counter] = isotopesByName.at("C12 PARENT");
        updatedPeakGroups[counter] = monoisotope;
        counter++;
    }

    for (auto it = isotopeToPeaks.begin(); it != isotopeToPeaks.end(); ++it) {

        string isotopeName = it->first;
        Isotope isotope = isotopesByName[isotopeName];
        updatedIsotopes[counter] = isotope;

        PeakGroup pg;
        pg.meanMz = isotope.mz;
        pg.expectedAbundance = isotope.abundance;
        pg.isotopeC13count = isotope.C13;
        pg.isotopicIndex = counter;
        pg.compound = this->compound;
        pg.adduct = this->adduct;
        pg.tagString = isotopeName;
        pg.setType(PeakGroup::GroupType::IsotopeType);
        pg.peaks = it->second; //peaks

        updatedPeakGroups[counter] = pg;

        counter++;
    }

    this->isotopes = updatedIsotopes;
    this->isotopePeakGroups = updatedPeakGroups;

    //TODO: recompute the isotopic envelopes?
}

float DifferentialIsotopicEnvelopeUtils::compareDifferentialIsotopicEnvelopes(
    vector<PeakGroup>& isotopePeakGroups,
    vector<mzSample*> unlabeledSamples,
    vector<mzSample*> labeledSamples,
    const IsotopeParameters& params,
    bool debug
    ){

    vector<PeakGroup> childrenSortedByMz = isotopePeakGroups;
    sort(childrenSortedByMz.begin(), childrenSortedByMz.end(), [](PeakGroup& lhs, PeakGroup& rhs){
        return lhs.meanMz < rhs.meanMz;
    });

    //agglomerated values
    vector<float> unlabeledIsotopesEnvelope{};
    vector<float> labeledIsotopesEnvelope{};

    //original values
    //(i, j) index matching between quant value and sample.
    vector<vector<float>> unlabeledIsotopeValuesEnvelope{};
    vector<vector<float>> labeledIsotopeValuesEnvelope{};
    vector<vector<mzSample*>> unlabeledIsotopeSamplesEnvelope{};
    vector<vector<mzSample*>> labeledIsotopeSamplesEnvelope{};

    map<mzSample*, float> sampleTotals{};
    for (auto& sample : unlabeledSamples) {
        sampleTotals.insert(make_pair(sample, 0.0f));
    }
    for (auto& sample : labeledSamples) {
        sampleTotals.insert(make_pair(sample, 0.0f));
    }

    for (unsigned int i = 0; i < childrenSortedByMz.size(); i++) {

        PeakGroup pg = childrenSortedByMz.at(i);

        if (debug) {
            cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: "
                 << pg.tagString
                 << " (m/z=" << pg.meanMz << ")"
                 << endl;
        }

        vector<float> unlabeledIsotopeValues{};
        vector<mzSample*> unlabeledIsotopeSamples{};

        vector<float> labeledIsotopeValues{};
        vector<mzSample*> labeledIsotopeSamples{};

        for (Peak p : pg.peaks) {

           if (find(unlabeledSamples.begin(), unlabeledSamples.end(), p.getSample()) != unlabeledSamples.end()){
                float quantVal = p.getQuantByName(params.diffIsoQuantType);
                if (quantVal > 0) {
                    mzSample* sample = p.getSample();
                    unlabeledIsotopeValues.push_back(quantVal);
                    unlabeledIsotopeSamples.push_back(sample);
                    sampleTotals[sample] += quantVal;
                }
                if (debug) {
                    cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: "
                         << pg.tagString
                         << "UNLABELED: '"
                         << p.getSample()->sampleName << "': quant="
                         << quantVal
                         << endl;
                }
           }

           if (find(labeledSamples.begin(), labeledSamples.end(), p.getSample()) != labeledSamples.end()) {
                float quantVal = p.getQuantByName(params.diffIsoQuantType);
                if (quantVal > 0) {
                    mzSample* sample = p.getSample();
                    labeledIsotopeValues.push_back(quantVal);
                    labeledIsotopeSamples.push_back(sample);
                    sampleTotals[sample] += quantVal;
                }
                if (debug) {
                    cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: "
                         << pg.tagString
                         << "  LABELED: '"
                         << p.getSample()->sampleName << "': quant="
                         << quantVal
                         << endl;
                }
           }
        }

        // Issue 725: Insufficient sample case
        // Isotopic extraction starts with peaks identified in the [M+0] peak group.
        // If the [M+0] peak group doesn't have enough measurements, then no other isotopes will,
        // which might produce a vector of all zeros. This can cause strange results in the
        // correlation score. Instead, just return 0 (no incorporation).
        if (i == 0 &&
                ((unlabeledIsotopeValues.size() < params.diffIsoReproducibilityThreshold) ||
                       (labeledIsotopeValues.size() < params.diffIsoReproducibilityThreshold))) {
           return 0;
        }

        unlabeledIsotopeValuesEnvelope.push_back(unlabeledIsotopeValues);
        unlabeledIsotopeSamplesEnvelope.push_back(unlabeledIsotopeSamples);

        labeledIsotopeValuesEnvelope.push_back(labeledIsotopeValues);
        labeledIsotopeSamplesEnvelope.push_back(labeledIsotopeSamples);

        float unlabeledIntensity = 0.0f;

        if (unlabeledIsotopeValues.size() >= params.diffIsoReproducibilityThreshold) {
           if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
                unlabeledIntensity = median(unlabeledIsotopeValues);
           } else if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
                unlabeledIntensity = accumulate(unlabeledIsotopeValues.begin(), unlabeledIsotopeValues.end(), 0.0f) / unlabeledIsotopeValues.size();
           } else if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Sum) {
                unlabeledIntensity = accumulate(unlabeledIsotopeValues.begin(), unlabeledIsotopeValues.end(), 0.0f);
           } else if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Max) {
                unlabeledIntensity = *max_element(unlabeledIsotopeValues.begin(), unlabeledIsotopeValues.end());
           }
        }

        float labeledIntensity = 0.0f;

        if (labeledIsotopeValues.size() >= params.diffIsoReproducibilityThreshold) {
           if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Mean) {
                labeledIntensity = median(labeledIsotopeValues);
           } else if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Median) {
                labeledIntensity = accumulate(labeledIsotopeValues.begin(), labeledIsotopeValues.end(), 0.0f) / labeledIsotopeValues.size();
           } else if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Sum) {
                labeledIntensity = accumulate(labeledIsotopeValues.begin(), labeledIsotopeValues.end(), 0.0f);
           } else if (params.diffIsoAgglomerationType == Fragment::ConsensusIntensityAgglomerationType::Max) {
                labeledIntensity = *max_element(labeledIsotopeValues.begin(), labeledIsotopeValues.end());
           }
        }

        //double zero
        if (unlabeledIntensity <= 0 && labeledIntensity <= 0) {
           if (params.diffIsoIncludeDoubleZero) {
                unlabeledIsotopesEnvelope.push_back(unlabeledIntensity);
                labeledIsotopesEnvelope.push_back(labeledIntensity);
           }

        //single zero
        } else if (unlabeledIntensity <= 0 || labeledIntensity <= 0) {
           if (params.diffIsoIncludeSingleZero) {
                unlabeledIsotopesEnvelope.push_back(unlabeledIntensity);
                labeledIsotopesEnvelope.push_back(labeledIntensity);
           }

        //no zero values
        } else {
           unlabeledIsotopesEnvelope.push_back(unlabeledIntensity);
           labeledIsotopesEnvelope.push_back(labeledIntensity);
        }

        if (debug) {
           cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: "
                << pg.tagString
                << " UNLABELED: {";
           for (unsigned int i = 0; i < unlabeledIsotopeValues.size(); i++) {
                if (i > 0) cout << ", ";
                cout << unlabeledIsotopeValues[i];
           }
           cout << "} ==> " << unlabeledIntensity << "; LABELED: {";
           for (unsigned int i = 0; i < labeledIsotopeValues.size(); i++) {
                if (i > 0) cout << ", ";
                cout << labeledIsotopeValues[i];
           }
           cout << "} ==> " << labeledIntensity << ";" << endl;
        }

    }

    if (debug) {
        cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: "
             << "ENVELOPE UNLABELED: {";
            for (unsigned int i = 0; i < unlabeledIsotopesEnvelope.size(); i++) {
                if (i > 0) cout << ", ";
                cout << unlabeledIsotopesEnvelope[i];
            }
            cout << "}; ENVELOPE LABELED: {";
            for (unsigned int i = 0; i < labeledIsotopesEnvelope.size(); i++) {
                if (i > 0) cout << ", ";
                cout << labeledIsotopesEnvelope[i];
            }
            cout << "};" << endl;

    }

    // If fewer than 2 measurements in the envelopes,
    // the two envelopes are not very interesting, and so
    // just return a score of 0 (i.e., no deviation in envelopes)
    if (unlabeledIsotopesEnvelope.size() < 2) {
        return 0;
    }

    float score = 0.0f;

    //deviation from r^2
    if (params.diffIsoScoringType == DiffIsoScoringType::PEARSON_CORRELATION) {
        float r = mzUtils::correlation(unlabeledIsotopesEnvelope, labeledIsotopesEnvelope);
        score = 1.0f - r*r;
    } else if (params.diffIsoScoringType == DiffIsoScoringType::F_STATISTIC) {

        //Interested in relative proportions in this case. So, need to normalize appropriately.

        vector<vector<float>> unlabeledIsotopeValuesEnvelopeNorm = vector<vector<float>>(unlabeledIsotopeValuesEnvelope.size());

        for (unsigned int i = 0; i < unlabeledIsotopeValuesEnvelope.size(); i++) {

            vector<float> ithValues = unlabeledIsotopeValuesEnvelope[i];
            vector<mzSample*> ithSamples = unlabeledIsotopeSamplesEnvelope[i];

            vector<float> normIthValues = vector<float>(ithValues.size());

            for (unsigned int j = 0; j < ithValues.size(); j++) {
                float jthValue = ithValues[j];
                mzSample* sample = ithSamples[j];

                float normFactor = sampleTotals.at(sample);

                normIthValues[j] = jthValue / normFactor;
            }

            unlabeledIsotopeValuesEnvelopeNorm[i] = normIthValues;
        }

        vector<vector<float>> labeledIsotopeValuesEnvelopeNorm = vector<vector<float>>(labeledIsotopeValuesEnvelope.size());

        for (unsigned int i = 0; i < labeledIsotopeValuesEnvelope.size(); i++) {

            vector<float> ithValues = labeledIsotopeValuesEnvelope[i];
            vector<mzSample*> ithSamples = labeledIsotopeSamplesEnvelope[i];

            vector<float> normIthValues = vector<float>(ithValues.size());

            for (unsigned int j = 0; j < ithValues.size(); j++) {
                float jthValue = ithValues[j];
                mzSample* sample = ithSamples[j];

                float normFactor = sampleTotals.at(sample);

                normIthValues[j] = jthValue / normFactor;
            }

            labeledIsotopeValuesEnvelopeNorm[i] = normIthValues;
        }

        if (debug) {
            cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: "
                 << "SOURCE ENVELOPES UNLABELED:\n";
            for (unsigned int i = 0; i < unlabeledIsotopeValuesEnvelopeNorm.size(); i++) {
                vector<float> unlabeledIsotopeValuesNorm = unlabeledIsotopeValuesEnvelopeNorm.at(i);
                cout << "{";
                for (unsigned int j = 0; j < unlabeledIsotopeValuesNorm.size(); j++) {
                    if (j > 0) cout << ", ";
                    cout << unlabeledIsotopeValuesNorm[i];
                }
                cout << "}; ";
            }
            cout << endl;

            cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: "
                 << "SOURCE ENVELOPES LABELED:\n";
            for (unsigned int i = 0; i < labeledIsotopeValuesEnvelopeNorm.size(); i++) {
                vector<float> labeledIsotopeValuesNorm = labeledIsotopeValuesEnvelopeNorm.at(i);
                cout << "{";
                for (unsigned int j = 0; j < labeledIsotopeValuesNorm.size(); j++) {
                    if (j > 0) cout << ", ";
                    cout << labeledIsotopeValuesNorm[i];
                }
                cout << "}; ";
            }
            cout << endl;
        }

        score = DifferentialIsotopicEnvelopeUtils::scoreByFStatistic(
            unlabeledIsotopeValuesEnvelopeNorm,
            labeledIsotopeValuesEnvelopeNorm,
            params,
            debug);
    }

    if (debug) {
            cout << "[IsotopicEnvelopeEvaluator::differentialIsotopicEnvelopes()]: score="
                 << score
                 << endl;
    }

    return score;
}

float DifferentialIsotopicEnvelopeUtils::scoreByFStatistic(
    vector<vector<float>> unlabeledIsotopeValuesEnvelope,
    vector<vector<float>> labeledIsotopeValuesEnvelope,
    const IsotopeParameters& params,
    bool debug
    ) {

    //Envelopes need to be able to be compared.
    if (unlabeledIsotopeValuesEnvelope.size() != labeledIsotopeValuesEnvelope.size()) {
        return 0.0f;
    }

    float F_sum = 0.0f;

    for (unsigned int i = 0; i < unlabeledIsotopeValuesEnvelope.size(); i++) {

        vector<float> unlabeledIsotope = unlabeledIsotopeValuesEnvelope.at(i);
        vector<float> labeledIsotope = labeledIsotopeValuesEnvelope.at(i);

        // Avoid NaNs
        if (unlabeledIsotope.empty() || labeledIsotope.empty()) {
            if (debug) {
                cout << "DifferentialIsotopicEnvelopeUtils::scoreByFStatistic(): i=" << i << ":\n"
                     << ", unlabeledIsotope.size(): " << unlabeledIsotope.size()
                     << ", labeledIsotope.size(): "<< labeledIsotope.size() << "\n";
            }
            continue;
        }

        vector<float> allObservations = unlabeledIsotope;
        allObservations.insert(allObservations.end(), labeledIsotope.begin(), labeledIsotope.end());

        double sum = std::accumulate(allObservations.begin(), allObservations.end(), 0.0);
        float mean_Total = sum / allObservations.size();

        unsigned int n_A = unlabeledIsotope.size();
        double sumUnlabeled = std::accumulate(unlabeledIsotope.begin(), unlabeledIsotope.end(), 0.0);
        float mean_A = sumUnlabeled / n_A;
        float var_A = mzUtils::variance(unlabeledIsotope);

        if (var_A <= 0) {
            if (debug) {
                cout << "DifferentialIsotopicEnvelopeUtils::scoreByFStatistic(): i=" << i << ":\n"
                     << "mean_total=" << mean_Total
                     << ", n_A=" << n_A << ", mean_A=" << mean_A << ", var_A=" << var_A
                     << endl;
            }
            continue;
        }

        unsigned int n_B = labeledIsotope.size();
        double sumLabeled = std::accumulate(labeledIsotope.begin(), labeledIsotope.end(), 0.0);
        float mean_B = sumLabeled / n_B;
        float var_B = mzUtils::variance(labeledIsotope);

        if (var_B <= 0) {
            if (debug) {
                cout << "DifferentialIsotopicEnvelopeUtils::scoreByFStatistic(): i=" << i << ":\n"
                     << "mean_total=" << mean_Total
                     << ", n_A=" << n_A << ", mean_A=" << mean_A << ", var_A=" << var_A
                     << ", n_B=" << n_B << ", mean_B=" << mean_B << ", var_B=" << var_B << "\n"
                     << endl;
            }
            continue;
        }

        unsigned int k = 2;

        float MSB_i = (n_A * (mean_A - mean_Total))*(n_A * (mean_A - mean_Total)) + (n_B * (mean_B - mean_Total))*(n_B * (mean_B - mean_Total));
        float MSW_i = ((n_A - 1) * var_A + (n_B - 1) * var_B) / (n_A + n_B - k);

        F_sum += MSB_i/MSW_i;

        if (debug) {
                cout << "DifferentialIsotopicEnvelopeUtils::scoreByFStatistic(): i=" << i << ":\n"
                     << "mean_total=" << mean_Total
                     << ", n_A=" << n_A << ", mean_A=" << mean_A << ", var_A=" << var_A
                     << ", n_B=" << n_B << ", mean_B=" << mean_B << ", var_B=" << var_B << "\n"
                     << "MSB=" << MSB_i << ", MSW=" << MSW_i << ", F_sum=" << F_sum
                     << endl;
        }
    }

    return F_sum;
}
