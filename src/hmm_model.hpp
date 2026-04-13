#ifndef HMMALIGN_HMM_MODEL_HPP
#define HMMALIGN_HMM_MODEL_HPP

#include <vector>
#include <string>
#include <limits>

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
    double match_score = 1.0;
    double mismatch_score = -1.0;
};

// --- NEW: Buffer for memory reuse ---
struct AlignmentBuffer {
    std::vector<double> M, I, D;
    std::vector<int8_t> ptrM, ptrI, ptrD;

    void reserve(size_t total) {
        if (M.size() < total) {
            M.assign(total, -std::numeric_limits<double>::infinity());
            I.assign(total, -std::numeric_limits<double>::infinity());
            D.assign(total, -std::numeric_limits<double>::infinity());
            ptrM.assign(total, -1);
            ptrI.assign(total, -1);
            ptrD.assign(total, -1);
        } else {
            // Re-zero/reset existing memory efficiently
            std::fill(M.begin(), M.begin() + total, -std::numeric_limits<double>::infinity());
            std::fill(I.begin(), I.begin() + total, -std::numeric_limits<double>::infinity());
            std::fill(D.begin(), D.begin() + total, -std::numeric_limits<double>::infinity());
            std::fill(ptrM.begin(), ptrM.begin() + total, -1);
            std::fill(ptrI.begin(), ptrI.begin() + total, -1);
            std::fill(ptrD.begin(), ptrD.begin() + total, -1);
        }
    }
};

} // namespace hmmalign
#endif