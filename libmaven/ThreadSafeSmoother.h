#pragma once
#include <vector>

namespace mzUtils
{

/**
* @brief The VectorSmoother class
* Thread-safe implementations of 1D smoothing (std::vector is thread-safe)
*/
class VectorSmoother {
public:

   unsigned long windowSize;
   std::vector<float> weights; // computed by getWeights().

   explicit VectorSmoother(unsigned long windowSize) { }
   virtual ~VectorSmoother() = 0;

   std::vector<float> smooth(std::vector<float>& data);
   unsigned long adjustWindowSize(unsigned long windowSize);

   virtual void computeWeights() = 0;

};

//inline is required here: https://stackoverflow.com/questions/23780274/duplicate-symbol-error-with-base-class-when-compiling
inline VectorSmoother::~VectorSmoother() { } //need definition, even if empty

/**
* @brief MovingAverageSmoother class (subclass of VectorSmoother)
* equal weights, all based on window size.
*/
class MovingAverageSmoother : public VectorSmoother {
public:

    explicit MovingAverageSmoother(unsigned long windowSize) : VectorSmoother(windowSize) {
        this->windowSize = adjustWindowSize(windowSize);
        computeWeights();
    }

    ~MovingAverageSmoother() { }
    void computeWeights();
};

/**
* @brief GaussianSmoother (subclass of VectorSmoother)
* weights sample Gaussian distribution (approximately),
* with furthest-out points pinned to +/- 3 sigma.
*/
class GaussianSmoother : public VectorSmoother {

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
    double zMax;

    /**
     * @brief sigma
     * standard deviation. Gaussian always has mean value of 0.
     */
    double sigma;

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

    /**
     * @brief FWHM_sigma
     *
     * Used in conjunction with getWindowSizeFromNsr()
     */
    constexpr static double FWHM_sigma = 1.17741002252; // sqrt(2 * ln(2))

    explicit GaussianSmoother(unsigned long windowSize) : GaussianSmoother(windowSize, 3, 1) { }

    explicit GaussianSmoother(unsigned long windowSize, double zMaxVal, double sigma) : VectorSmoother(windowSize) {

        //Here, 'windowSize' refers to 'nsr' - the number of points with > 50% amplitude
        this->windowSize = adjustWindowSize(windowSize);

        GaussianSmoother::init(zMaxVal, sigma);

        computeWeights();
    }

    ~GaussianSmoother() { }

    /**
     *
     * @brief getWeights
     *
     * Previously, Gaussian smoothing was specified in nsr space
     * (see mzUtils::gaussian1d_smoothing), and so the weights vector
     * is also specified this way.
     *
     * Note that the actual length of the weights vector depends on
     * the zMax parameter (how far out to extend the Gaussian in each direction
     * for convolution)
     *
     * @param windowSize (# of points with amplitude > half max)
     * or, the # of points in with 2 * FWHM_sigma
     *
     * @return weights vector
     * size depends on nsr (windowSize) and zMax parameter
     */
    void computeWeights();

    /**
     * @brief getGaussianWeight
     * get a single weight based on a z-score (number of std dev)
     */
    double getGaussianWeight(double sigma);

private:
    void init(double zMaxVal, double sigma);
};

}
