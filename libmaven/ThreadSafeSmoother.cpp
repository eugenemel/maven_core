#include "ThreadSafeSmoother.h"
#include <cmath>
#include <iostream>

using namespace std;

vector<float> VectorSmoother::smooth(vector<float> data, vector<float> weights){
    vector<float> smoothedData = vector<float>(data.size(), 0);

    int halfWindow = (weights.size() - 1) / 2;

    for (int i = 0; i < data.size(); i++){
        int jMin = i-halfWindow;
        int jMax = i+halfWindow;

        if (jMin < 0 || jMax >= data.size()) {
            continue;
            //TODO: handle edges?
        }

        int weightIndex = 0;
        for (int j = jMin; j <= jMax; j++) {
            smoothedData.at(j) = smoothedData.at(j) + (data.at(j) * weights.at(weightIndex));
            weightIndex++;
        }

    }

    return smoothedData;
}

vector<float> MovingAverageSmoother::getWeights(unsigned long windowSize){
    float frac = 1 / static_cast<float>(windowSize);
    return vector<float>(windowSize, frac);
}

GaussianSmoother::GaussianSmoother(unsigned long zMax, unsigned long sigma){
    GaussianSmoother::init(zMax, sigma);
}

GaussianSmoother::GaussianSmoother() {
    GaussianSmoother::init(3, 1);
}

void GaussianSmoother::init(unsigned long zMax, unsigned long sigma){

    //passed arguments
    this->zMax = zMax;
    this->sigma = sigma;

    //computed constants
    this->k1 = 1 / (static_cast<double>(sigma) * sqrt(2 * M_PI));
    this->k2 = 1 / (2 * static_cast<double>(sigma) * static_cast<double>(sigma));
}

vector<float> GaussianSmoother::getWeights(unsigned long windowSize) {

    vector<float> weights = vector<float>(windowSize, 0);
    float weightSum = 0;

    unsigned long halfWindow = static_cast<unsigned long>(windowSize-1)/2;
    double deltaSigma = zMax / static_cast<double>(halfWindow+1); // endpoint at zMax sigma is not directly used

    unsigned long index = 0;

    for (unsigned long i = 0; i < halfWindow; i++) {

        double zScore = static_cast<double>(halfWindow-i)*deltaSigma;
        float gaussianWeight = static_cast<float>(getGaussianWeight(zScore));

        weights.at(index) = gaussianWeight;
        weightSum += gaussianWeight;

        index++;
    }

    double zScore = 0;
    float gaussianWeight = static_cast<float>(getGaussianWeight(zScore));

    weights.at(index) = gaussianWeight;
    weightSum += gaussianWeight;

    index++;

    for (unsigned long i = 0; i < halfWindow; i++) {

        double zScore = static_cast<double>(i+1) * deltaSigma;
        float gaussianWeight = static_cast<float>(getGaussianWeight(zScore));

        weights.at(index) = gaussianWeight;
        weightSum += gaussianWeight;

        index++;
    }

    vector<float> normalizedWeights = vector<float>(windowSize, 0);

    for (unsigned int i = 0; i < weights.size(); i++){
        normalizedWeights.at(i) = weights.at(i) / weightSum;
    }

    return normalizedWeights;

}

double GaussianSmoother::getGaussianWeight(double zScore) {
    return k1 * exp(-1 * k2 * (zScore * zScore));
}

/**
 * For Testing
 *
 * To compile:
 * cd ~/workspace/maven_core/libmaven
 * clang-omp++ ThreadSafeSmoother.cpp -o ThreadSafeSmoother -Wall -I/usr/local/Cellar/llvm/8.0.1/include/c++/v1
 *
 * Execute:
 * ./ThreadSafeSmoother
 *
 * Functions tested
 * -- MovingAverageSmoother::getWeights() [2019-09-03]
 * -- GaussianSmoother::getWeights() [2019-09-03]
 * -- MovingAverageSmoother::smooth() [TODO]
 * @brief main
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {

    GaussianSmoother gaussianSmoother = GaussianSmoother(3, 1);
    MovingAverageSmoother movingAverageSmoother = MovingAverageSmoother();

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
}
