#pragma once
#include <vector>

//class
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

    class MovingAverageSmoother : VectorSmoother {
        std::vector<float> getWeights(int windowSize);
    };

    class GaussianSmoother : VectorSmoother {
        std::vector<float> getWeights(int windowSize);
    };
}
