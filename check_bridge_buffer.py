import hmmalign._core as core

def test_buffered_alignment():
    # 1. Initialize parameters
    params = core.ModelParams()
    params.M_to_M = 0.0
    params.M_to_I = -2.0
    params.I_to_I = -0.5
    params.M_to_D = -2.0
    params.D_to_D = -0.5
    params.I_to_M = -0.1
    params.D_to_M = -0.1
    params.bandwidth = 100

    # 2. Setup sequences
    ref = "ACGTACGT" # 8bp Reference
    perfect_read = "ACGTACGT"
    noisy_read   = "ACGTTACGT" # 9bp Read (Insertion of 'T' at index 4)

    # 3. Create the Persistent Buffer
    # This memory is allocated once and reused for all subsequent calls
    buf = core.AlignmentBuffer()

    print("Running Buffered HMM Alignment Tests...")
    print("-" * 40)

    # 4. Test Perfect Match
    res_perfect = core.viterbi_align_buffered(perfect_read, ref, params, buf)
    print(f"Perfect Match Score: {res_perfect.score:>8.2f}")
    print(f"Perfect Match CIGAR: {res_perfect.path}")

    # 5. Test Noisy Match (Insertion)
    # The buffer 'buf' will be resized if necessary and reused here
    res_noisy = core.viterbi_align_buffered(noisy_read, ref, params, buf)
    print(f"Noisy Match Score:   {res_noisy.score:>8.2f}")
    print(f"Noisy Match CIGAR:   {res_noisy.path}")

    # 6. Verification
    expected_perfect_cigar = "8M"
    expected_noisy_cigar = "3M1I5M"

    success = True
    if res_perfect.path != expected_perfect_cigar:
        print(f"❌ Perfect match failed. Expected {expected_perfect_cigar}, got {res_perfect.path}")
        success = False
    
    if res_noisy.path != expected_noisy_cigar:
        print(f"❌ Noisy match failed. Expected {expected_noisy_cigar}, got {res_noisy.path}")
        success = False

    if success:
        print("-" * 40)
        print("✅ Verification Complete: Memory-buffered alignment is accurate.")
        print("✅ Success! The HMM correctly identified state transitions using reused memory.")

if __name__ == "__main__":
    test_buffered_alignment()