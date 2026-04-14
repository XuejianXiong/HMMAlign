import json
import time
import argparse
import hmmalign._core as hmm

def load_fasta(path):
    try:
        with open(path, 'r') as f:
            # Skip headers, join lines, force uppercase for the C++ kernel
            return "".join(line.strip() for line in f if not line.startswith(">")).upper()
    except FileNotFoundError:
        print(f"Error: File not found at {path}")
        return None

def main():
    # 1. Setup CLI Argument Parsing
    parser = argparse.ArgumentParser(description="HMMAlign Suite Benchmarking Tool")
    parser.add_argument("--input", required=True, help="Path to input.json")
    parser.add_argument("--config", required=True, help="Path to config.json")
    args = parser.parse_args()

    # 2. Load Configurations
    try:
        with open(args.config, "r") as f:
            config = json.load(f)
        
        with open(args.input, "r") as f:
            inputs = json.load(f)
    except Exception as e:
        print(f"Failed to load JSON files: {e}")
        return

    # 3. Setup Reference and Query from input.json data block
    # Ensure the input.json has "data": {"reference": "...", "sample": "..."}
    ref_path = inputs.get("data", {}).get("reference")
    sample_path = inputs.get("data", {}).get("sample")
    
    if not ref_path or not sample_path:
        print("Error: input.json must contain data.reference and data.sample paths.")
        return

    # Load reference and query reads
    ref_seq = load_fasta(ref_path)
    query_seq = load_fasta(sample_path)
    
    if not ref_seq or not query_seq:
        return

    # 4. Map JSON Parameters to C++ ModelParams
    params = hmm.ModelParams()
    p = config["parameters"]
    params.M_to_M = p.get("M_to_M", 1.0)
    params.M_to_I = p.get("M_to_I", -2.0)
    params.M_to_D = p.get("M_to_D", -5.0)
    params.I_to_M = p.get("I_to_M", -0.5)
    params.I_to_I = p.get("I_to_I", -0.1)
    params.D_to_M = p.get("D_to_M", -0.5)
    params.D_to_D = p.get("D_to_D", -0.1)
    params.bandwidth = p.get("bandwidth", 500)

    # 5. Execution with Buffer Reuse
    buf = hmm.AlignmentBuffer()
    start_time = time.perf_counter()
    
    result = hmm.viterbi_align_buffered(query_seq, ref_seq, params, buf)
    
    end_time = time.perf_counter()

    # 6. Reporting Output
    elapsed_ms = (end_time - start_time) * 1000
    print(f"\n" + "="*50)
    print(f"🚀 HMMAlign Suite | Model: {config.get('model_name', 'Unknown')}")
    print(f"="*50)
    print(f"Ref:    {ref_path} ({len(ref_seq)} bp)")
    print(f"Sample: {sample_path} ({len(query_seq)} bp)")
    print(f"Time:   {elapsed_ms:.2f} ms")
    print(f"Score:  {result.score:.2f}")
    
    # Clean output for the CIGAR
    if result.path:
        display_path = result.path if len(result.path) < 100 else f"{result.path[:97]}..."
        print(f"CIGAR:  {display_path}")
    else:
        print("CIGAR:  No alignment found.")
    print("="*50 + "\n")

if __name__ == "__main__":
    main()