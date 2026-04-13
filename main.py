import hmmalign._core as core
import time
import json
import argparse
import os


def load_config(path):
    if os.path.exists(path):
        with open(path, 'r') as f:
            return json.load(f)
    return {}

def run_test_case(name, read, ref, expected_op=None, params=None):
    if params is None:
        params = core.ModelParams()
        params.bandwidth = 200

    buf = core.AlignmentBuffer()
    start = time.time()
    result = core.viterbi_align_buffered(read, ref, params, buf)
    duration = time.time() - start

    print(f"🧪 Test: {name}")
    print(f"  Score: {result.score:.2f} | Time: {duration*.1000:.3f}ms")
    print(f"  CIGAR: {result.path}")
    
    if expected_op and expected_op in result.path:
        print("  ✅ Passed")
    elif expected_op:
        print(f"  ❌ Failed (Expected '{expected_op}' in CIGAR)")
    else:
        print("  ✅ Result Logged")
    print("-" * 50)

def test_suite(config_data):
    print("🚀 Starting Comprehensive HMMAlign Test Suite\n")

    # 1. Standard Identity Match
    run_test_case("Perfect Match", "ACGT" * 25, "ACGT" * 25, "100M")

    # 2. N-Neutrality Test
    # Read has an 'N', Ref has 'A'. Should score 9.0 (10 bases - 1 neutral)
    # under match=1.0, mismatch=-1.5
    run_test_case(
        "N-Neutrality (Read)", 
        "AAAAANAAAA", 
        "AAAAAAAAAA", 
        "10M"
    )

    # 3. Affine Deletion Test (Structural Variant)
    # Testing if a 50bp gap is consolidated into a single 'D' block
    params = core.ModelParams()
    params.M_to_D = config_data['parameters']['M_to_D']   # -5.0 Expensive Open
    params.D_to_D = config_data['parameters']['D_to_D']   # -0.1 Cheap Extend
    params.bandwidth = config_data['parameters']['bandwidth']   # 100
    
    ref_gap = "A" * 100 + "G" * 50 + "T" * 100
    read_gap = "A" * 100 + "T" * 100
    run_test_case("Affine Deletion (50bp)", read_gap, ref_gap, "50D", params)

    # 4. Glocal Start Test
    # Read should align to the MIDDLE of the reference without penalty
    ref_long = "GGGGGGGGGG" + "ATCG" * 10 + "CCCCCC"
    read_mid = "ATCG" * 10
    run_test_case("Glocal Alignment (Mid-stream)", read_mid, ref_long, "40M")

    # 5. High-Noise Stress Test (10kb)
    # Verifying the Mac Pro performance and buffer stability
    big_ref = "ACGT" * 2500  # 10,000 bp
    big_read = "ACGT" * 2490 + "NNNNNNNNNN" # 10,000 bp with N-tail
    run_test_case("10kb Long Read Stress", big_read, big_ref, "M")

    # 6. Random "N" Noise (Middle & Scattered)
    # 100bp with a "low quality" center
    ref_mid_n = "A" * 400 + "GCGT" * 5 + "A" * 400
    read_mid_n = "A" * 400 + "NNNN" * 5 + "A" * 400
    run_test_case("Scattered N-Noise (Middle)", read_mid_n, ref_mid_n, "820M")

    # 7. High-Density N-Island
    # Testing if the HMM can "bridge" a 20bp gap of total uncertainty
    ref_island = "C" * 1000 + "ATGC" * 5 + "G" * 1000
    read_island = "C" * 1000 + "N" * 20 + "G" * 1000
    run_test_case("N-Island Bridge (20bp)", read_island, ref_island, "2020M")

    # 8. The "N-Start" Edge Case
    # Some basecallers put Ns at the very beginning of a read
    run_test_case("N-Start Alignment", "NNNNNAAAA", "TTTTTAAAA", "9M")


if __name__ == "__main__":
    # Handle the --input argument from rebuild.sh
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    # Load the config.json
    config_data = load_config(args.config)
    #print(config_data)
    
    test_suite(config_data)