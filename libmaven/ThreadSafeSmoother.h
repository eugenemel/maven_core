#pragma once
#include <vector>

/**
* @brief The VectorSmoother class
* Thread-safe implementations of 1D smoothing (std::vector is thread-safe)
*/
class VectorSmoother {
public:
   virtual ~VectorSmoother() = 0;
   std::vector<float> smooth(std::vector<float> data, unsigned long windowSize){

       //intercept window size to minimum of 3 (1 before and 1 after)
       if (windowSize <= 0) {
           windowSize = 3;
       }

       if (windowSize % 2 == 0){
           windowSize++; //err on side of longer half-windows
       }

       return smooth(data, getWeights(windowSize));
   }

   virtual std::vector<float> getWeights(unsigned long windowSize) = 0;
   std::vector<float> smooth(std::vector<float> data, std::vector<float> weights);
};

VectorSmoother::~VectorSmoother() { } //need definition, even if empty

/**
* @brief MovingAverageSmoother class (subclass of VectorSmoother)
* equal weights, all based on window size.
*/
class MovingAverageSmoother : VectorSmoother {
public:
    ~MovingAverageSmoother() { }
   std::vector<float> getWeights(unsigned long windowSize);
};

/**
* @brief GaussianSmoother (subclass of VectorSmoother)
* weights sample Gaussian distribution (approximately),
* with furthest-out points pinned to +/- 3 sigma.
*/
class GaussianSmoother : VectorSmoother {
public:
    ~GaussianSmoother() { }
    std::vector<float> getWeights(unsigned long windowSize);
    double getGaussianWeight(double sigma);
};
