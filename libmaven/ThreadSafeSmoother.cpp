#include "ThreadSafeSmoother.h"
#include <cmath>
#include <iostream>
#include <map>

using namespace std;
namespace mzUtils {

unsigned long VectorSmoother::adjustWindowSize(unsigned long windowSize) {

    //intercept window size to minimum of 3 (1 before and 1 after)
    if (windowSize <= 0) {
        return 3;
    }

    //even numbers become odd
    if (windowSize % 2 == 0){
        return (windowSize+1);
    }

    return windowSize;
}

vector<float> VectorSmoother::smooth(vector<float>& data){

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

void MovingAverageSmoother::computeWeights(){
    float frac = 1 / static_cast<float>(windowSize);
    weights = vector<float>(windowSize, frac);
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
void GaussianSmoother::computeWeights() {

    unsigned long halfWindow = static_cast<unsigned long>(windowSize-1)/2;
    double deltaSigma = GaussianSmoother::FWHM_sigma / halfWindow;

    halfWindow = floor(zMax / deltaSigma) + 1;
    windowSize = 2*halfWindow + 1;

    weights = vector<float>(windowSize, 0);

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

    for (unsigned int i = 3; i <= 15; i=i+2){

        mzUtils::GaussianSmoother gaussianSmoother = mzUtils::GaussianSmoother(i, 3, 1);
        mzUtils::MovingAverageSmoother movingAverageSmoother = mzUtils::MovingAverageSmoother(i);

         vector<float> movingAvgWeights = movingAverageSmoother.weights;
         cout << "moving avg window=" << i << ": ";
         for (auto weight : movingAvgWeights) {
             cout << weight << " ";
         }

         cout << endl;

         vector<float> gaussianWeights = gaussianSmoother.weights;
         cout << "Gaussian   window=" << i << ": ";
         for (auto weight : gaussianWeights) {
             cout << weight << " ";
         }

         cout << endl;
    }
    cout << endl;

    cout << "TEST SMOOTHING" << endl << endl;

    mzUtils::GaussianSmoother gaussianSmoother = mzUtils::GaussianSmoother(5, 3, 1);
    mzUtils::MovingAverageSmoother movingAverageSmoother = mzUtils::MovingAverageSmoother(15);

    cout << "ONE: CONSTANT VECTOR" << endl;

    vector<float> allSevens = vector<float>(10, 7);
    vector<float> smoothedSevens = movingAverageSmoother.smooth(allSevens);
    vector<float> gaussianSevens = gaussianSmoother.smooth(allSevens);

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

    vector<float> smoothedFivesAndSevens = movingAverageSmoother.smooth(fivesAndSevens);
    vector<float> gaussianFivesAndSevens = gaussianSmoother.smooth(fivesAndSevens);

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

    vector<float> movingAvgSharpPeak = movingAverageSmoother.smooth(sharpPeak);
    vector<float> gaussianSharpPeak = gaussianSmoother.smooth(sharpPeak);

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
