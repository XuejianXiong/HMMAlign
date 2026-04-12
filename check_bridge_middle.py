import hmmalign._core as core
import time

def test_glocal_middle():
    params = core.ModelParams()
    
    # --- 1. Set Parameters with Positive Reinforcement ---
    # We need match_score to be a reward (e.g., 1.0) so the algorithm 
    # prefers a long 996bp path over a short 1bp path.
    # Note: Ensure your C++ EmissionParams reflect these rewards/penalties.
    params.M_to_M = 0.0    # Log(1.0) - stay in match
    params.M_to_I = -2.0   # Penalty to start gap
    params.I_to_I = -0.5   # Penalty to extend gap
    params.M_to_D = -2.0   # Penalty to start gap
    params.D_to_D = -0.5   # Penalty to extend gap
    params.I_to_M = -0.1   
    params.D_to_M = -0.1
    params.bandwidth = 500 # Wide enough to find the 100bp offset

    # --- 2. Create Sequences ---
    # Reference is 1200bp
    ref = "ACGT" * 300 
    
    # Read is 996bp, matches ref starting at index 100
    # This means read[0] aligns to ref[100]
    read = ref[100 : 100 + 996]
    
    print(f"Aligning Read ({len(read)} bp) to Reference ({len(ref)} bp)")
    print("Scenario: Read matches the MIDDLE (starting at index 100).")
    print("-" * 40)

    # --- 3. Run Wide Band Test ---
    params.bandwidth = 500
    start_time = time.time()
    res_wide = core.viterbi_align(read, ref, params)
    duration = time.time() - start_time
    
    print(f"Wide Band (500)   | Score: {res_wide.score:>8.2f} | Time: {duration:.5f}s")
    print(f"Path Length: {len(res_wide.path)}")

    # --- 4. Run Narrow Band Test ---
    # Must be > 100 to account for the offset between read[0] and ref[100]
    params.bandwidth = 150 
    res_narrow = core.viterbi_align(read, ref, params)
    print(f"Narrow Band (150) | Score: {res_narrow.score:>8.2f}")
    print(f"Path Length: {len(res_narrow.path)}")

    # --- 5. Verification ---
    # If match_score is 1.0, a 996bp perfect match should score 996.0
    # If your match_score is still 0.0, the length is the only way to verify success.
    if len(res_wide.path) == 996:
        print("\n✅ Success! Glocal identified the full middle match.")
        print(f"Path matches: {res_wide.path == res_narrow.path}")
        print(f"Path snippet (first 10): {res_wide.path[:10]}...")
    else:
        print("\n❌ Failure: Alignment was truncated.")
        print(f"Expected Length: 996 | Actual Length: {len(res_wide.path)}")
        print(f"Score: {res_wide.score}")
        print("Hint: Ensure match_score in C++ is a positive reward (e.g., 1.0).")

if __name__ == "__main__":
    test_glocal_middle()