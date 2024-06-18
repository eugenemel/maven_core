#pragma once

#ifndef ISOTOPICENVELOPEUTILS_H
#define ISOTOPICENVELOPEUTILS_H

#include "mzSample.h"
#include <numeric>


//Issue 673: makes exporting much easier
class IsotopeMatrix {
public:
    vector<string> sampleNames{};
    vector<string> isotopeNames{};
    MatrixXf isotopesData{}; // (rows = samples, columns = isotope names)

    static IsotopeMatrix getIsotopeMatrix(
        PeakGroup *group,
        PeakGroup::QType quantType,
        vector<mzSample*> samples,
        bool isNaturalAbundanceCorrected,
        bool isFractionOfSampleTotal,
        bool debug = false);
};

/**
 * @brief The IsotopicEnvelope class
 * Measurements for a set of isotopes in a single sample.
 * the corresponding isotope for measurements in the intensities field
 * should be stored in a separate IsotopicEnvelopeGroup class.
 */
class IsotopicEnvelope {
public:

    static IsotopicEnvelope& none();

    string source = "unknown";
    double totalIntensity = -1.0;
    vector<double> intensities;

    double getTotalIntensity();
    void print();
};

/**
 * @brief The IsotopicEnvelopeGroup class
 * All measurements must follow the isotopic pattern defined in the
 * isotopes field.
 * Note that this means that all sample measurements must match this pattern -
 * cannot have some samples with one set of isotopes and other samples with a different set
 * of isotopes.
 *
 * This supports that isotopes can be a compound-specific thing - but not a sample specific thing.
 * In general, all samples should ve generated with the same resolving power,
 * so you shouldn't expect different sets of isotopes for different samples.
 */
class IsotopicEnvelopeGroup {
public:

    /** Input fields **/

    Compound *compound = nullptr;
    Adduct *adduct = nullptr;
    PeakGroup *group = nullptr;
    vector<Isotope> isotopes{};

    /** Computed Values **/

    //View cross-isotope measurements by sample
    //Individual (Sample, Isotope) measurements are stored by elements
    //in the IsotopicEnvelope.intensities vector (same indices as IsotopicEnvelopeGroup.isotopes).
    map<mzSample*, IsotopicEnvelope> envelopeBySample{};

    //View cross-sample measurements by isotope (index matches IsotopicEnvelopeGroup.isotopes)
    //Individual (Isotope, Sample) measurements are stored as peaks in the PeakGroup
    vector<PeakGroup> isotopePeakGroups{};

    //Usually provided by IsotopicExtractionParameters
    string extractionAlgorithmName = "";

    //convenience method that sets children of the IsotopicEnvelope.group field
    //to the IsotopicEnvelopeGroup.children (same indices as IsotopicEnvelopeGroup.isotopes).
    void setIsotopesToChildrenPeakGroups(Classifier *classifier = nullptr);

    // Issue 691: combine isotopes that give the same quant readout.
    void combineOverlappingIsotopes(float ppm=20, bool debug=false);

    void print();
};

class IsotopicEnvelopeExtractor {
public:

    /**
     * @brief extractEnvelope
     * Most common entry point. The specific extraction algorithm will be selected based on the input params.
     *
     * Note that sample information is passed in - this supports cases where a sample of interest has no Peak in
     * inthe input PeakGroup.
     *
     * Isotopes are computed based on the input parameters and compound information.
     *
     * @param compound
     * @param adduct
     * @param group
     * @param samples
     * @param params
     * @param debug
     * @return
     */
    static IsotopicEnvelopeGroup extractEnvelopes(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<mzSample*>& samples,
        IsotopeParameters params,
        bool debug=false);

    //usually called from IsotopicEnvelopeExtractor::extractEnvelope()
    static IsotopicEnvelopeGroup extractEnvelopesFromMPlusZeroPeaks(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        IsotopeParameters params,
        bool debug=false);

    //The original approach implemented in MAVEN up through version 2.0 2023-09-25
    static IsotopicEnvelopeGroup extractEnvelopesVersion1(
        Compound *compound,
        Adduct *adduct,
        PeakGroup *group,
        vector<Isotope>& isotopes,
        IsotopeParameters params,
        bool debug=false);
};

class DifferentialIsotopicEnvelopeUtils {
public:
    /**
     * @brief
     * Given a vector<PeakGroup> corresponding to IsotopicEnvelopeGroup.isotopePeakGroups,
     * along with a vector<mzSample*> corresponding to unlabeled samples and a vector<mzSample*> corresponding to
     * labeled samples, evaluate the disparity between the isotopic envelope measured from unlabeled samples compared to
     * labeled samples.
     *
     */
    static float compareDifferentialIsotopicEnvelopes(
        vector<PeakGroup>& isotopePeakGroups,
        vector<mzSample*> unlabeledSamples,
        vector<mzSample*> labeledSamples,
        const IsotopeParameters& params,
        bool debug=false
        );

    /**
     *
     * @brief
     * Scoring appraoch for isotopic envelopes, comparing the variability within a set of envelopes
     * to the variability between two agglomerated envelopes.
     *
     *
     */
    static float scoreByFStatistic(
        vector<vector<float>> unlabeledIsotopeValuesEnvelope,
        vector<vector<float>> labeledIsotopeValuesEnvelope,
        const IsotopeParameters& params,
        bool debug=false
    );

    static IsotopeMatrix constructDiffIsotopeMatrix(
        PeakGroup *group,
        vector<mzSample*> unlabeledSamples,
        vector<mzSample*> labeledSamples,
        const IsotopeParameters& params,
        bool debug=false
        );

    static float scoreByPearsonCorrelationCoefficient(
        PeakGroup *group,
        vector<mzSample*> unlabeledSamples,
        vector<mzSample*> labeledSamples,
        const IsotopeParameters& params,
        bool debug=false
        );
};

#endif // ISOTOPICENVELOPEUTILS_H
