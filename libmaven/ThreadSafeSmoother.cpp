#include "ThreadSafeSmoother.h"
#include <cmath>
#include <iostream>
#include <map>

using namespace std;
namespace mzUtils {
vector<float> VectorSmoother::smooth(vector<float> data, vector<float> weights){

    int halfWindow = static_cast<int>((weights.size() - 1) / 2);
    int dataSize = static_cast<int>(data.size());

    vector<float> smoothedData = vector<float>(data.size(), 0);

    for (int i = 0; i < dataSize; i++){

        int jMin = i-halfWindow;
        int jMax = i+halfWindow;

        float weightSumm = 0;
        unsigned long weightIndex = 0;

        for (int j = jMin; j <= jMax; j++) {

            if (j >= 0 && j < data.size()) {

                float weight = weights.at(weightIndex);
                weightSumm += weight;

                smoothedData.at(i) += weight * data.at(j);
            }

            weightIndex++;

        }

        smoothedData.at(i) /= weightSumm;
    }

    return smoothedData;
}

vector<float> MovingAverageSmoother::getWeights(unsigned long windowSize){
    float frac = 1 / static_cast<float>(windowSize);
    return vector<float>(windowSize, frac);
}

GaussianSmoother::GaussianSmoother(double zMax, double sigma){
    GaussianSmoother::init(zMax, sigma);
}

GaussianSmoother::GaussianSmoother() {
    GaussianSmoother::init(3, 1);
}

void GaussianSmoother::init(double zMax, double sigma){

    //passed arguments
    this->zMax = zMax;
    this->sigma = sigma;

    //computed constants
    this->k1 = 1 / (static_cast<double>(sigma) * sqrt(2 * M_PI));
    this->k2 = 1 / (2 * static_cast<double>(sigma) * static_cast<double>(sigma));
}

//Here, 'windowSize' refers to 'nsr' - the number of points with > 50% amplitude
vector<float> GaussianSmoother::getWeights(unsigned long windowSize) {

    unsigned long halfWindow = static_cast<unsigned long>(windowSize-1)/2;
    double deltaSigma = GaussianSmoother::FWHM_sigma / halfWindow;

    halfWindow = floor(zMax / deltaSigma) + 1;
    windowSize = 2*halfWindow + 1;

    vector<float> weights = vector<float>(windowSize, 0);

    unsigned long index = 0;

    for (unsigned long i = 0; i < halfWindow; i++) {

        double zScore = static_cast<double>(halfWindow-i)*deltaSigma;
        float gaussianWeight = static_cast<float>(getGaussianWeight(zScore));

        weights.at(index) = gaussianWeight;

        index++;
    }

    double zScore = 0;
    float gaussianWeight = static_cast<float>(getGaussianWeight(zScore));

    weights.at(index) = gaussianWeight;

    index++;

    for (unsigned long i = 0; i < halfWindow; i++) {

        double zScore = static_cast<double>(i+1) * deltaSigma;
        float gaussianWeight = static_cast<float>(getGaussianWeight(zScore));

        weights.at(index) = gaussianWeight;

        index++;
    }

    return weights;

}

double GaussianSmoother::getGaussianWeight(double zScore) {
    return k1 * exp(-1 * k2 * (zScore * zScore));
}
}
/**
 * For Testing
 *
 * To create an executable:
 * rename 'testThreadSafeSmoother' to 'main'
 * cd ~/workspace/maven_core
 * clang-omp++ libmaven/ThreadSafeSmoother.cpp -o libmaven/ThreadSafeSmoother -Wall -I/usr/local/Cellar/llvm/8.0.1/include/c++/v1
 *
 * Execute:
 * ./ThreadSafeSmoother
 *
 * Functions tested
 * -- MovingAverageSmoother::getWeights() [2019-09-03]
 * -- GaussianSmoother::getWeights() [2019-09-03]
 * -- MovingAverageSmoother::smooth() [2019-09-03]
 * -- smooth() speedup [2019-09-04]
 * -- verifyVsKnownData [TODO]
 *
 * @brief testThreadSafeSmoother
 * @param argc
 * @param argv
 * @return
 */
int testThreadSafeSmoother(int argc, char *argv[]) {

    mzUtils::GaussianSmoother gaussianSmoother = mzUtils::GaussianSmoother(3, 1);
    mzUtils::MovingAverageSmoother movingAverageSmoother = mzUtils::MovingAverageSmoother();

    for (unsigned int i = 3; i <= 15; i=i+2){

         vector<float> movingAvgWeights = movingAverageSmoother.getWeights(i);
         cout << "moving avg window=" << i << ": ";
         for (auto weight : movingAvgWeights) {
             cout << weight << " ";
         }

         cout << endl;

         vector<float> gaussianWeights = gaussianSmoother.getWeights(i);
         cout << "Gaussian   window=" << i << ": ";
         for (auto weight : gaussianWeights) {
             cout << weight << " ";
         }

         cout << endl;
    }
    cout << endl;

    cout << "TEST SMOOTHING" << endl << endl;

    cout << "ONE: CONSTANT VECTOR" << endl;

    vector<float> allSevens = vector<float>(10, 7);
    vector<float> smoothedSevens = movingAverageSmoother.smooth(allSevens, 5);
    vector<float> gaussianSevens = gaussianSmoother.smooth(allSevens, 5);

    for (auto f : allSevens) {
        cout << f << " ";
    }
    cout << endl;

    for (auto f : smoothedSevens){
        cout << f << " ";
    }
    cout << endl;

    for (auto f : gaussianSevens){
        cout << f << " ";
    }
    cout << endl;
    cout << endl;

    cout << "TWO: ALTERNATING VECTOR" << endl;

    vector<float> fivesAndSevens = vector<float>(10, 7);
    for (int i = 0; i < 10; i++){
        if (i % 2 == 0){
            fivesAndSevens.at(i) = 5;
        }
    }

    vector<float> smoothedFivesAndSevens = movingAverageSmoother.smooth(fivesAndSevens, 3);
    vector<float> gaussianFivesAndSevens = gaussianSmoother.smooth(fivesAndSevens, 3);

    for (auto f : fivesAndSevens) {
        cout << f << " ";
    }
    cout << endl;

    for (auto f : smoothedFivesAndSevens){
        cout << f << " ";
    }
    cout << endl;

    for (auto f : gaussianFivesAndSevens){
        cout << f << " ";
    }
    cout << endl;
    cout << endl;

    cout << "THREE: SHARP PEAK" << endl;

    vector<float> sharpPeak = vector<float>(21, 0);
    for (int i = 0; i < 20; i++){
        if (i >=5 && i < 10) {
            sharpPeak.at(i) = 10-(10-i)*2;
        } else if (i >= 10 && i < 15) {
            sharpPeak.at(i) = 10-(i-10)*2;
        }
    }

    vector<float> movingAvgSharpPeak = movingAverageSmoother.smooth(sharpPeak, 15);
    vector<float> gaussianSharpPeak = gaussianSmoother.smooth(sharpPeak, 5);

    string separator = " " ;

    for (auto f : sharpPeak) {
        cout << f << separator;
    }
    cout << endl;

    for (auto f : movingAvgSharpPeak){
        cout << f << separator;
    }
    cout << endl;

    for (auto f : gaussianSharpPeak){
        cout << f << separator;
    }
    cout << endl;
    cout << endl;
}
