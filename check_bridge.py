import hmmalign._core as core

try:
    # This calls the placeholder function we wrote in viterbi.cpp
    result = core.align_test("ACGT", "ACGT")
    print(f"✅ Success! C++ returned: {result}")
except Exception as e:
    print(f"❌ Error: {e}")