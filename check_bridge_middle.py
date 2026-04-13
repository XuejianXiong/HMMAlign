import hmmalign._core as core
import time

def test_glocal_middle():
    params = core.ModelParams()
    
    # --- 1. Parameters ---
    # Log-probs: Ensure match reward exists in C++ EmissionParams
    params.M_to_M = 0.0
    params.M_to_I = -2.0
    params.I_to_I = -0.5
    params.M_to_D = -2.0
    params.D_to_D = -0.5
    params.I_to_M = -0.1
    params.D_to_M = -0.1

    # --- 2. Sequence Creation ---
    ref = "ACGT" * 300  # 1200bp
    read = ref[100 : 100 + 996] # 996bp matching middle
    
    print(f"Aligning Read ({len(read)} bp) to Reference ({len(ref)} bp)")
    print("Scenario: Read matches the MIDDLE (starting at index 100).")
    print("-" * 40)

    # --- 3. Wide Band Test ---
    params.bandwidth = 500
    start_time = time.time()
    res_wide = core.viterbi_align(read, ref, params)
    duration = time.time() - start_time
    
    print(f"Wide Band (500)   | Score: {res_wide.score:>8.2f} | Time: {duration:.5f}s")
    print(f"CIGAR Output: {res_wide.path}")

    # --- 4. Narrow Band Test ---
    params.bandwidth = 200 
    res_narrow = core.viterbi_align(read, ref, params)
    print(f"Narrow Band (200) | Score: {res_narrow.score:>8.2f}")
    print(f"CIGAR Output: {res_narrow.path}")

    # --- 5. Verification ---
    # In CIGAR format, a perfect 996bp match is the string "996M"
    expected_cigar = f"{len(read)}M"
    
    if res_wide.path == expected_cigar:
        print(f"\n✅ Success! Glocal identified the full match.")
        print(f"Bands match: {res_wide.path == res_narrow.path}")
        print(f"Final CIGAR: {res_wide.path}")
    else:
        print("\n❌ Failure: CIGAR does not match expectation.")
        print(f"Expected: {expected_cigar} | Actual: {res_wide.path}")
        print(f"Path string length: {len(res_wide.path)} (This is the char count of the CIGAR string)")

if __name__ == "__main__":
    test_glocal_middle()