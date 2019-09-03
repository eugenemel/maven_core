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

       if (windowSize % 2 != 0){
           windowSize++; //err on side of longer half-windows
       }

       return smooth(data, getWeights(windowSize));
   }
private:
   virtual std::vector<float> getWeights(unsigned long windowSize) = 0;
   std::vector<float> smooth(std::vector<float> data, std::vector<float> weights);
};

/**
* @brief MovingAverageSmoother class (subclass of VectorSmoother)
* equal weights, all based on window size.
*/
class MovingAverageSmoother : VectorSmoother {
    ~MovingAverageSmoother();
private:
   std::vector<float> getWeights(unsigned long windowSize);
};

/**
* @brief GaussianSmoother (subclass of VectorSmoother)
* weights sample Gaussian distribution (approximately),
* with furthest-out points pinned to +/- 3 sigma.
*/
class GaussianSmoother : VectorSmoother {
    ~GaussianSmoother();
private:
   double getGaussianWeight(double sigma);
   std::vector<float> getWeights(unsigned long windowSize);
};
