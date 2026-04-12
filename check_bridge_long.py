import hmmalign._core as core
import time

def test_banded_alignment():
    params = core.ModelParams()
    
    # --- 1. Set Log-Probability Parameters ---
    # We use log-space: 0.0 is probability 1.0, negative is a penalty.
    params.M_to_M = 0.0    # High probability to stay in match
    params.M_to_I = -2.0   # Penalty for starting an insertion
    params.I_to_I = -0.5   # Penalty for extending an insertion
    params.M_to_D = -2.0   # Penalty for starting a deletion
    params.D_to_D = -0.5   # Penalty for extending a deletion
    params.I_to_M = -0.1   # Penalty to transition back
    params.D_to_M = -0.1   # Penalty to transition back

    # --- 2. Create Sequences ---
    # Reference is 1200bp
    ref = "ACGT" * 300   
    
    # Read is 996bp, matching the START of the reference
    read_base = ref[:996]
    read_list = list(read_base)
    # Introduce exactly one mismatch at index 500
    read_list[500] = 'T' if read_list[500] != 'T' else 'A' 
    read = "".join(read_list)
    
    print(f"Aligning Read ({len(read)} bp) to Reference ({len(ref)} bp)")
    print("Scenario: Read matches the beginning of the reference perfectly.")
    print("-" * 40)

    # --- 3. Run Wide Band Test ---
    params.bandwidth = 500
    start_wide = time.time()
    res_wide = core.viterbi_align(read, ref, params)
    end_wide = time.time() - start_wide
    print(f"Wide Band (500)  | Score: {res_wide.score:>8.2f} | Time: {end_wide:.5f}s")

    # --- 4. Run Narrow Band Test ---
    # 20 is plenty of room for a single mismatch on a 1000bp sequence
    params.bandwidth = 20 
    start_narrow = time.time()
    res_narrow = core.viterbi_align(read, ref, params)
    end_narrow = time.time() - start_narrow
    print(f"Narrow Band (20) | Score: {res_narrow.score:>8.2f} | Time: {end_narrow:.5f}s")

    # --- 5. Verification ---
    if res_wide.score == res_narrow.score and res_wide.score != float('-inf'):
        print("\n✅ Success! Real numbers achieved and bands match.")
        print(f"Path matches: {res_wide.path == res_narrow.path}")
        # Snippet shows the 'M' states (Match/Mismatch)
        print(f"Path snippet: ...{res_narrow.path[498:503]}...")
    else:
        print("\n❌ Still getting -inf or mismatch.")
        if res_wide.score == float('-inf'):
            print("Hint: Even the Wide Band failed. Check your C++ transition logic.")
        else:
            print(f"Hint: Wide score was {res_wide.score}, but Narrow was {res_narrow.score}.")

    reduction = (1 - (end_narrow / end_wide)) * 100
    print(f"Efficiency Gain: {reduction:.1f}%")

if __name__ == "__main__":
    test_banded_alignment()