#include "ThreadSafeSmoother.h"
#include <cmath>
#include <iostream>

using namespace std;

vector<float> VectorSmoother::smooth(vector<float> data, vector<float> weights){
    vector<float> smoothedData = vector<float>(data.size(), 0);

    //TODO

    return smoothedData;
}

vector<float> MovingAverageSmoother::getWeights(unsigned long windowSize){
    float frac = 1 / static_cast<float>(windowSize);
    return vector<float>(windowSize, frac);
}

vector<float> GaussianSmoother::getWeights(unsigned long windowSize) {

    double windowSizeDouble = static_cast<double>(windowSize);

    vector<float> weights = vector<float>(windowSize, 0);

    double deltaSigma = 6 / windowSizeDouble; // bounds determined by +/- 3 sigma

    unsigned int index = 0;

    while (index < windowSize) {
        double sigma = -3 + deltaSigma * index;
        weights.at(index) = static_cast<float>(getGaussianWeight(sigma)/windowSizeDouble);
        index++;
    }

    return weights;

}

double GaussianSmoother::getGaussianWeight(double sigma) {
    double expVal = -1 / pow(2 * sigma, 2);
    double divider = sqrt(2 * M_PI * pow(sigma, 2));
    return (1 / divider) * exp(expVal);
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
 * -- GaussianSmoother::getWeights() [TODO]
 * @brief main
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {

    GaussianSmoother gaussianSmoother = GaussianSmoother();
    MovingAverageSmoother movingAverageSmoother = MovingAverageSmoother();

    for (unsigned int i = 3; i <= 10; i++){

         vector<float> movingAvgWeights = movingAverageSmoother.getWeights(i);
         cout << "window=" << i << ": ";
         for (auto weight : movingAvgWeights) {
             cout << weight << " ";
         }

         cout << endl;

         vector<float> gaussianWeights = gaussianSmoother.getWeights(i);
         cout << "window=" << i << ": ";
         for (auto weight : gaussianWeights) {
             cout << weight << " ";
         }

         cout << endl;
    }
}
