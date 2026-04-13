#ifndef HMMALIGN_HMM_MODEL_HPP
#define HMMALIGN_HMM_MODEL_HPP

#include <vector>
#include <string>
#include <limits>
#include <algorithm>

namespace hmmalign {

struct ModelParams {
    // Transitions (Log-space)
    float M_to_M = 0.0f;    // Match stays Match
    float M_to_I = -4.5f;   // Gap Open (Insertion): High penalty
    float M_to_D = -4.5f;   // Gap Open (Deletion): High penalty
    
    float I_to_M = -0.6f;   // Gap Close
    float I_to_I = -0.2f;   // Gap Extend: Low penalty (Affine)
    
    float D_to_M = -0.6f;   // Gap Close
    float D_to_D = -0.2f;   // Gap Extend: Low penalty (Affine)
    
    int bandwidth = 500;
};

struct EmissionParams {
    float match_score = 1.0f;
    float mismatch_score = -1.5f; // Punish mismatches more than gap extensions
};

struct AlignmentBuffer {
    std::vector<float> M, I, D;
    std::vector<int8_t> ptrM, ptrI, ptrD;

    void reserve(size_t total) {
        float neg_inf = -std::numeric_limits<float>::infinity();
        if (M.size() < total) {
            M.assign(total, neg_inf);
            I.assign(total, neg_inf);
            D.assign(total, neg_inf);
            ptrM.assign(total, -1);
            ptrI.assign(total, -1);
            ptrD.assign(total, -1);
        } else {
            std::fill(M.begin(), M.begin() + total, neg_inf);
            std::fill(I.begin(), I.begin() + total, neg_inf);
            std::fill(D.begin(), D.begin() + total, neg_inf);
            std::fill(ptrM.begin(), ptrM.begin() + total, -1);
            std::fill(ptrI.begin(), ptrI.begin() + total, -1);
            std::fill(ptrD.begin(), ptrD.begin() + total, -1);
        }
    }
};

} // namespace hmmalign
#endif