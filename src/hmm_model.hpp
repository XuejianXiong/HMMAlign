#ifndef HMMALIGN_HMM_MODEL_HPP
#define HMMALIGN_HMM_MODEL_HPP

#include <string>

namespace hmmalign {

struct ModelParams {
    double M_to_M = 0.0;
    double M_to_I = -3.0;
    double M_to_D = -3.0;
    double I_to_M = -0.5;
    double I_to_I = -1.0;
    double D_to_M = -0.5;
    double D_to_D = -1.0;
    int bandwidth = 500;
};

struct EmissionParams {
    double match_score = 0.0;
    double mismatch_score = -2.0;
};

} // namespace hmmalign

#endif // HMMALIGN_HMM_MODEL_HPP
