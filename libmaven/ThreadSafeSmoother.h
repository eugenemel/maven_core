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

    /**
     * @brief zMax
     * Sample Gaussian intensities from (-zMax, +zMax) evenly
     *
     * Gaussian does not include intensity values from +/- zMax itself, but will asymptotically approach this value
     * as the windowSize increases.
     *
     * Examples:
     * windowSize=3, zMax=3
     * --> z= (-1.5, 0, +1.5)
     *
     * windowSize=5, zMax=3
     * --> z = (2, 1, 0, 1, 2)
     *
     * windowSize=15, zMax=3
     * --> z = (-2.625, -2.25, -1.875, -1.5, -1.125, -0.75, -0.375, 0, 0.375, 0.75, 1.125, 1.5, 1.875, 2.25, 2.625)
     *
     * windowSize=15, zMax=4
     * --> z = (3.5, 3, 2.5, 2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)
     */
    unsigned long zMax;

    /**
     * @brief sigma
     * standard deviation. Gaussian always has mean value of 0.
     */
    unsigned long sigma;

    /**
     * @brief k1
     *
     * Gaussian distribution constant:
     *
     * 1 / (sigma * sqrt(2 * PI))
     */
    double k1;

    /**
     * @brief k2
     *
     * Gaussian distribution constant:
     *
     * 1 / (2 * sigma * sigma)
     */
    double k2;

public:

    GaussianSmoother();
    GaussianSmoother(unsigned long zMaxVal, unsigned long sigma);

    ~GaussianSmoother() { }
    std::vector<float> getWeights(unsigned long windowSize);
    double getGaussianWeight(double sigma);

private:
    void init(unsigned long zMaxVal, unsigned long sigma);
};
