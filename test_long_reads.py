import hmmalign._core as core
import time
import random

def generate_noisy_read(reference, length, error_rate=0.05):
    """Simulates a long read with random mismatches/indels."""
    start_pos = 100
    #read_seq = list(reference[start_pos : start_pos + length])
    # Force a 5bp deletion
    read_seq = list(reference[start_pos : start_pos + length//3] 
                    + reference[start_pos + (length//3)+5 : start_pos + length])
    
    for i in range(len(read_seq)):
        if random.random() < error_rate:
            # 80% mismatches, 20% small indels (simplified)
            read_seq[i] = random.choice(['A', 'C', 'G', 'T'])
            
    return "".join(read_seq)

def benchmark_10kb():
    params = core.ModelParams()
    params.bandwidth = 500  # Crucial for long reads
    params.M_to_I = -2.0
    params.M_to_D = -2.0
    
    # 1. Setup 12kb Reference and 10kb Read
    ref = "".join(random.choice(['A', 'C', 'G', 'T']) for _ in range(12000))
    read = generate_noisy_read(ref, 10000, error_rate=0.05)
    
    print(f"🚀 Starting 10kb Long Read Test")
    print(f"Read Length: {len(read)} bp | Ref Length: {len(ref)} bp")
    print("-" * 45)

    # 2. Create Buffer
    buf = core.AlignmentBuffer()

    # 3. Execution
    start_time = time.time()
    result = core.viterbi_align_buffered(read, ref, params, buf)
    end_time = time.time() - start_time

    # 4. Results
    print(f"Alignment Time: {end_time:.4f} seconds")
    print(f"Viterbi Score:  {result.score:.2f}")
    
    # Check CIGAR for length consistency
    # (Extracting total matched/inserted bases from CIGAR)
    #print(f"CIGAR Snippet:  {result.path[:20]}...{result.path[-20:]}")
    print(f"Full CIGAR: {result.path}")
    
    if result.score > 0:
        print("\n✅ Success: 10kb read aligned within the band.")
    else:
        print("\n❌ Failure: Score is too low or -inf. Increase bandwidth?")

def test_sensitivity():
    # 1. Setup: Reference has a 50bp "G" block that the read lacks
    ref_part1 = "A" * 500
    gap_chunk = "G" * 50  # The structural variation (Deletion)
    ref_part2 = "T" * 500
    reference = ref_part1 + gap_chunk + ref_part2
    
    # Read is a continuous match to the reference but missing the 'G' block
    read = ref_part1 + ref_part2 
    
    # 2. Configure "Senior Scientist" Affine Parameters
    params = core.ModelParams()
    params.M_to_D = -5.0   # High Gap Open penalty
    params.D_to_D = -0.1   # Low Gap Extend penalty (favors long gaps)
    params.bandwidth = 200 # Sufficient for a 50bp drift
    
    buf = core.AlignmentBuffer()
    
    print(f"🚀 Running Affine Gap Sensitivity Test")
    print(f"Read Length: {len(read)} | Ref Length: {len(reference)}")
    print(f"Expected: A 50D block in the middle.")
    print("-" * 45)

    # 3. Execution
    start = time.time()
    result = core.viterbi_align_buffered(read, reference, params, buf)
    end = time.time() - start

    # 4. Results
    print(f"Alignment Time: {end:.4f}s")
    print(f"Viterbi Score:  {result.score:.2f}")
    print(f"Full CIGAR:     {result.path}")

    # 5. Validation Logic
    if "50D" in result.path:
        print("\n✅ Success: The HMM correctly identified the 50bp Affine Deletion.")
    else:
        print("\n❌ Failure: The gap was fragmented or incorrectly aligned.")


if __name__ == "__main__":
    #benchmark_10kb()
    test_sensitivity()