import hmmalign._core as core
import time
import json
import argparse
import os

def load_config(path):
    if path and os.path.exists(path):
        with open(path, 'r') as f:
            return json.load(f)
    return {}

def run_test_case(name, read, ref, expected_op=None, params=None):
    if params is None:
        params = core.ModelParams()
        params.bandwidth = 200

    buf = core.AlignmentBuffer()
    
    # Use perf_counter for microsecond precision
    start = time.perf_counter()
    result = core.viterbi_align_buffered(read, ref, params, buf)
    duration_ms = (time.perf_counter() - start) * 1000

    print(f"🧪 Test: {name}")
    print(f"  Score: {result.score:.2f} | Time: {duration_ms:.4f}ms")
    # Truncate long CIGARs for readability
    display_cigar = result.path if len(result.path) < 60 else f"{result.path[:30]}...{result.path[-20:]}"
    print(f"  CIGAR: {display_cigar}")
    
    if expected_op and expected_op in result.path:
        print("  ✅ Passed")
    elif expected_op:
        print(f"  ❌ Failed (Expected '{expected_op}' in CIGAR)")
    else:
        print("  ✅ Result Logged")
    print("-" * 50)

def test_suite(config_data):
    print("🚀 Starting HMMAlign Suite (Soft-Clipping & HPC Optimized)\n")
    
    # Extract parameters for specialized tests
    p = config_data.get('parameters', {})
    
    # 1. Standard Identity Match
    run_test_case("Perfect Match", "ACGT" * 25, "ACGT" * 25, "100M")

    # 2. Affine Deletion Test
    params_affine = core.ModelParams()
    params_affine.M_to_D = p.get('M_to_D', -5.0)
    params_affine.D_to_D = p.get('D_to_D', -0.1)
    params_affine.bandwidth = p.get('bandwidth', 100)
    
    ref_gap = "A" * 100 + "G" * 50 + "T" * 100
    read_gap = "A" * 100 + "T" * 100
    run_test_case("Affine Deletion (50bp)", read_gap, ref_gap, "50D", params_affine)

    # 3. Soft-Clipping: Leading Adapter
    # Read has 10bp of 'G' that doesn't exist in 'A' reference
    run_test_case(
        "Leading Soft-Clip (Adapter)", 
        "GGGGGGGGGG" + "A" * 50, 
        "A" * 100, 
        "10S"
    )

    # 4. Soft-Clipping: Trailing Adapter
    # Read ends with 15bp of 'C' that doesn't exist in 'T' reference
    run_test_case(
        "Trailing Soft-Clip", 
        "T" * 50 + "CCCCCCCCCCCCCCC", 
        "T" * 100, 
        "15S"
    )

    # 5. Soft-Clipping: Dual Ends (The "Real World" case)
    # Junk at both ends, biological match in the middle
    run_test_case(
        "Dual-End Soft-Clip", 
        "GGGG" + "ACGT" * 10 + "CCCC", 
        "AAAAA" + "ACGT" * 10 + "TTTTT", 
        "4S"
    )

    # 6. Glocal Start Test
    run_test_case("Glocal (Mid-stream)", "ATCG" * 10, "GGGGGGGGGG" + "ATCG" * 10 + "CCCCCC", "40M")

    # 7. High-Noise Stress Test (10kb)
    big_ref = "ACGT" * 2500 
    big_read = "ACGT" * 2490 + "NNNNNNNNNN" 
    run_test_case("10kb Long Read Stress", big_read, big_ref, "M")

    # 8. Scatterned N-Noise
    run_test_case("Scattered N-Noise", "A" * 400 + "NNNN" * 5 + "A" * 400, "A" * 400 + "GCGT" * 5 + "A" * 400, "M")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="Path to config.json")
    args = parser.parse_args()

    # Load parameters
    config_data = load_config(args.config)
    
    test_suite(config_data)