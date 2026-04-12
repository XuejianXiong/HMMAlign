from hmmalign import _core as core

# 1. Initialize parameters
params = core.ModelParams()
params.M_to_M = 0.0    # Log-space
params.M_to_I = -3.0
params.M_to_D = -3.0

# 2. Define a test case (Exact match vs one with an insertion)
ref = "ACGTACGT"
read_perfect = "ACGTACGT"
read_noisy = "ACGTTACGT" # Extra T at index 4

# 3. Run alignment
result_perfect = core.viterbi_align(read_perfect, ref, params)
result_noisy = core.viterbi_align(read_noisy, ref, params)

print(f"Perfect Match Score: {result_perfect.score}")
print(f"Path:  {result_perfect.path}")

print(f"Noisy Match Score: {result_noisy.score}")
print(f"Path:  {result_noisy.path}")


assert result_perfect.score > result_noisy.score
print("✅ Verification Complete: Perfect match scored higher than noisy match.")


# Verify the extra T was caught as an 'I' (Insertion)
if 'I' in result_noisy.path:
    print("✅ Success! The HMM identified an insertion in the noisy read.")
