from hmmalign import _core as core

# 1. Initialize parameters
params = core.ModelParams()
params.M_to_M = 0.0    # Log-space
params.M_to_I = -3.0
params.M_to_D = -3.0

# 2. Define a test case (Exact match vs one with an insertion)
ref = "ACGTACGT"
read_perfect = "ACGTACGT"
read_noisy = "ACGTTACGT" # Extra T

# 3. Run alignment
score_perfect = core.viterbi_align(read_perfect, ref, params)
score_noisy = core.viterbi_align(read_noisy, ref, params)

print(f"Perfect Match Score: {score_perfect}")
print(f"Noisy Match Score: {score_noisy}")

assert score_perfect > score_noisy
print("✅ Verification Complete: Perfect match scored higher than noisy match.")