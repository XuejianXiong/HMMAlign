import json
import time
import argparse
import os
import hmmalign._core as hmm

def load_fasta(path):
    """
    Parses a FASTA file. 
    Returns: (Header, Sequence)
    """
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
            if not lines:
                return "Unknown", ""
            header = lines[0].strip().lstrip('>')
            # Join all lines after the header and force uppercase
            seq = "".join(line.strip() for line in lines[1:]).upper()
            return header, seq
    except Exception as e:
        print(f"Error loading FASTA at {path}: {e}")
        return "Unknown", ""

def write_sam(output_path, read_id, ref_id, cigar, seq, score, ref_len):
    """
    Generates a standard SAM file for IGV/Bioinformatics tool compatibility.
    """
    # SAM Header: @SQ is required for IGV to know the reference length
    header = (
        f"@HD\tVN:1.6\tSO:unsorted\n"
        f"@SQ\tSN:{ref_id}\tLN:{ref_len}\n"
    )
    
    # SAM Record
    # FLAG 0 = Mapped; MAPQ 60 = High Confidence; POS 1 = 1-based start
    # AS:f: is the tag for the Alignment Score
    pos = 1 
    mapq = 60
    record = f"{read_id}\t0\t{ref_id}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*\tAS:f:{score}\n"
    
    with open(output_path, "w") as f:
        f.write(header)
        f.write(record)

def main():
    # 1. CLI Argument Parsing
    parser = argparse.ArgumentParser(description="HMMAlign Suite - Main Execution Engine")
    parser.add_argument("--input", required=True, help="Path to input.json")
    parser.add_argument("--config", required=True, help="Path to config.json")
    args = parser.parse_args()

    # 2. Load JSON Configurations
    try:
        with open(args.config, "r") as f:
            config = json.load(f)
        with open(args.input, "r") as f:
            inputs = json.load(f)
    except Exception as e:
        print(f"Failed to load JSON configurations: {e}")
        return

    # 3. Load Genomic Data
    ref_path = inputs.get("data", {}).get("reference")
    sample_path = inputs.get("data", {}).get("sample")
    
    if not ref_path or not sample_path:
        print("Error: input.json must define data.reference and data.sample")
        return

    ref_header, ref_seq = load_fasta(ref_path)
    read_header, read_seq = load_fasta(sample_path)
    
    if not ref_seq or not read_seq:
        print("Error: Could not load sequences. Check file paths.")
        return

    # 4. Initialize C++ Model Parameters
    params = hmm.ModelParams()
    p = config.get("parameters", {})
    
    # Using .get() with defaults to prevent crashes if JSON keys are missing
    params.M_to_M = p.get("M_to_M", 1.0)
    params.M_to_I = p.get("M_to_I", -2.0)
    params.M_to_D = p.get("M_to_D", -5.0)
    params.I_to_M = p.get("I_to_M", -0.5)
    params.I_to_I = p.get("I_to_I", -0.1)
    params.D_to_M = p.get("D_to_M", -0.5)
    params.D_to_D = p.get("D_to_D", -0.1)
    params.bandwidth = p.get("bandwidth", 500)

    # 5. Core Execution (Viterbi with Buffer Reuse)
    buf = hmm.AlignmentBuffer()
    start_time = time.perf_counter()
    
    result = hmm.viterbi_align_buffered(read_seq, ref_seq, params, buf)
    
    end_time = time.perf_counter()
    elapsed_ms = (end_time - start_time) * 1000

    # 6. Output Generation (SAM File)
    out_dir = inputs["execution"].get("output_dir", "./results")
    out_name = inputs["execution"].get("output_file", "alignment.sam")
    os.makedirs(out_dir, exist_ok=True)
    
    sam_path = os.path.join(out_dir, out_name)
    
    # Clean up IDs for SAM compatibility (strip spaces)
    read_id = read_header.split()[0]
    ref_id = ref_header.split()[0]
    
    write_sam(sam_path, read_id, ref_id, result.path, read_seq, result.score, len(ref_seq))

    # 7. Final Console Report
    print(f"\n" + "="*60)
    print(f"🚀 HMMAlign | Engine: {config.get('model_name', 'v1')}")
    print(f"="*60)
    print(f"Reference: {ref_id} ({len(ref_seq)} bp)")
    print(f"Query:     {read_id} ({len(read_seq)} bp)")
    print(f"Status:    Success")
    print(f"Time:      {elapsed_ms:.2f} ms")
    print(f"Score:     {result.score:.2f}")
    print(f"Output:    {sam_path}")
    
    if result.path:
        display_cigar = result.path if len(result.path) < 80 else f"{result.path[:77]}..."
        print(f"CIGAR:     {display_cigar}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()