#ifndef HMMALIGN_HMM_MODEL_HPP
#define HMMALIGN_HMM_MODEL_HPP

#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <cstdint>

/**
 * @namespace hmmalign
 * @brief Core engine for High-Performance Hidden Markov Model sequence alignment.
 */
namespace hmmalign {

/**
 * @struct ModelParams
 * @brief Defines the transition probabilities for the Affine Gap HMM.
 * * All values are stored in log-space to allow for additive scoring and to prevent 
 * floating-point underflow. Higher values indicate higher probability.
 */
struct ModelParams {
    // --- Transitions (Log-space) ---
    
    /// Probability of staying in a Match/Mismatch state.
    float M_to_M = 0.0f;    
    
    /// Gap Opening Penalty for Insertions.
    float M_to_I = -4.5f;   
    
    /// Gap Opening Penalty for Deletions.
    float M_to_D = -4.5f;   
    
    /// Probability of returning to Match state from an Insertion.
    float I_to_M = -0.6f;   
    
    /// Gap Extension Penalty for Insertions (Affine).
    float I_to_I = -0.2f;   
    
    /// Probability of returning to Match state from a Deletion.
    float D_to_M = -0.6f;   
    
    /// Gap Extension Penalty for Deletions (Affine).
    float D_to_D = -0.2f;   
    
    /** * @brief Constraint for Sakoe-Chiba banding.
     * Limits the DP calculation to a diagonal band of size 2*bandwidth + 1.
     * Reduces complexity from O(N*M) to O(N*bandwidth).
     */
    int bandwidth = 500;
};

/**
 * @struct EmissionParams
 * @brief Scoring matrix constants for biological base matching.
 */
struct EmissionParams {
    float match_score = 1.0f;
    float mismatch_score = -1.5f; 
};

/**
 * @struct AlignmentBuffer
 * @brief Reusable memory pool for Viterbi DP matrices.
 * * To avoid the high cost of frequent heap allocations (malloc/free) in high-throughput 
 * genomics, this buffer is passed by reference and resized only when necessary.
 */
struct AlignmentBuffer {
    /// Score matrices for the three HMM states.
    std::vector<float> M, I, D;
    
    /// Backtrace pointers for path reconstruction.
    std::vector<int8_t> ptrM, ptrI, ptrD;

    /**
     * @brief Efficiently prepares memory for a new alignment task.
     * @param total The total required size (usually read_len * ref_len or band_size).
     */
    void reserve(size_t total) {
        constexpr float neg_inf = -std::numeric_limits<float>::infinity();
        
        if (M.size() < total) {
            // Reallocation required
            M.assign(total, neg_inf);
            I.assign(total, neg_inf);
            D.assign(total, neg_inf);
            ptrM.assign(total, -1);
            ptrI.assign(total, -1);
            ptrD.assign(total, -1);
        } else {
            // Reuse existing memory to minimize heap pressure
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