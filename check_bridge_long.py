import hmmalign._core as core
import time

def test_banded_alignment():
    params = core.ModelParams()
    
    # 1. Create identical lengths to satisfy Global Alignment
    # 250 * 4 = 1000 bases
    ref = "ACGT" * 250 
    
    # Create a read that is also 1000 bases, but with a mismatch in the middle
    # We change one 'A' to a 'T' at index 500
    read_list = list(ref)
    read_list[500] = 'T' 
    read = "".join(read_list)
    
    print(f"Aligning Read ({len(read)} bp) to Reference ({len(ref)} bp)")
    print("-" * 40)

    # Test 1: Wide Band (Calculation covers almost the whole matrix)
    params.bandwidth = 500
    start_wide = time.time()
    res_wide = core.viterbi_align(read, ref, params)
    end_wide = time.time() - start_wide
    
    # Test 2: Narrow Band (Calculation covers only 1% of the matrix)
    params.bandwidth = 10
    start_narrow = time.time()
    res_narrow = core.viterbi_align(read, ref, params)
    end_narrow = time.time() - start_narrow

    # 2. Results Analysis
    print(f"Wide Band (500)  | Score: {res_wide.score:>8.2f} | Time: {end_wide:.5f}s")
    print(f"Narrow Band (10) | Score: {res_narrow.score:>8.2f} | Time: {end_narrow:.5f}s")

    # 3. Verification
    if res_wide.score == res_narrow.score and res_wide.score != float('-inf'):
        print("\n✅ Success: Scores are real numbers and identical!")
        print(f"Path matches: {res_wide.path == res_narrow.path}")
        # Show the mismatch in the path (should be all 'M')
        print(f"Path snippet around mismatch: ...{res_narrow.path[498:503]}...")
    else:
        print("\n❌ Calculation failed or still returning -inf.")
        print("Check if sequences are exactly the same length for Global Alignment.")

    reduction = (1 - (end_narrow / end_wide)) * 100
    print(f"Efficiency Gain: {reduction:.1f}%")

if __name__ == "__main__":
    test_banded_alignment()