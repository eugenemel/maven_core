#include "ThreadSafeSmoother.h"
#include <cmath>

namespace mzUtils{

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

    double getGaussianWeight(double sigma) {
        double expVal = -1 / pow(2 * sigma, 2);
        double divider = sqrt(2 * M_PI * pow(sigma, 2));
        return (1 / divider) * exp(expVal);
    }
}
