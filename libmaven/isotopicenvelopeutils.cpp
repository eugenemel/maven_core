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

    vector<Isotope> theoreticalIsotopes = MassCalculator::computeIsotopes2(
        compound->formula,
        adduct,
        params.getLabeledIsotopes(),
        params.labeledIsotopeRetentionPolicy,
        NaturalAbundanceData::defaultNaturalAbundanceData,
        params.isNatAbundance,
        maxNumProtons,
        params.natAbundanceThreshold,
        false
        );

    vector<Isotope> isotopes = theoreticalIsotopes;
    if (params.isCondenseTheoreticalIsotopes) {
        isotopes = IsotopicEnvelopeAdjuster::condenseTheoreticalIsotopes(
            theoreticalIsotopes,
            params,
            false
            );
    }

    IsotopicEnvelopeGroup envelopeGroup;

    if (params.isotopicExtractionAlgorithm == IsotopicExtractionAlgorithm::PEAK_FULL_RT_BOUNDS_AREA) {
        envelopeGroup = extractEnvelopesPeakFullRtBounds(compound, adduct, group, isotopes, params, debug);
    } else if (params.isotopicExtractionAlgorithm == IsotopicExtractionAlgorithm::MAVEN_GUI_VERSION_ONE) {
        envelopeGroup = extractEnvelopesVersion1(compound, adduct, group, isotopes, params, debug);
    }

    //Issue 691: combine overlapping isotopes when the data is identical.
    if (params.isCombineOverlappingIsotopes) {
        envelopeGroup.combineOverlappingIsotopes(params.ppm);
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

void IsotopicEnvelopeGroup::combineOverlappingIsotopes(float ppm) {

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

    //identify all of the cases where the same quant value is identified by multiple isotopes.
    set<vector<string>> combinations{};
    for (auto it = quantValToIsotopes.begin(); it != quantValToIsotopes.end(); ++it) {
        if (it->second.size() > 1) {
            combinations.insert(it->second);
        }
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
