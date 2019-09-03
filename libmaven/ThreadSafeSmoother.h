#pragma once
#include <vector>

namespace mzUtils {

     /**
     * @brief The VectorSmoother class
     * Thread-safe implementations of 1D smoothing (std::vector is thread-safe)
     */
    class VectorSmoother {
    public:
        std::vector<float> smooth(std::vector<float> data, std::vector<float> weights);
        virtual std::vector<float> getWeights(int windowSize) = 0;
        virtual ~VectorSmoother() = 0;
    };

    /**
     * @brief MovingAverageSmoother class (subclass of VectorSmoother)
     * equal weights, all based on window size.
     */
    class MovingAverageSmoother : VectorSmoother {
        std::vector<float> getWeights(int windowSize);
    };

    /**
     * @brief GaussianSmoother (subclass of VectorSmoother)
     * weights sample Gaussian distribution (approximately),
     * with furthest-out points pinned to +/- 3 sigma.
     */
    class GaussianSmoother : VectorSmoother {
        std::vector<float> getWeights(int windowSize);
    private:
        double getGaussianWeight(double sigma);
    };
}
