# ⬢ HMMAlign
## High-Performance Viterbi Engine for Reference-Guided Sequence Alignment

## ◈ Overview:

HMMAlign is a high-performance computational framework designed for the precise mapping of genomic reads to reference sequences. By leveraging a 3-state Hidden Markov Model (HMM) and a Banded Viterbi algorithm, the engine effectively resolves indel noise and translocation speed variability inherent in modern sequencing technologies.

 - Architecture: C++20 Core Inference Kernel with Pybind11 integration for MLOps orchestration.

 - Feature Set: $O(N \times B)$ Banded Alignment, N-Neutrality handling, and automated CIGAR-compliant traceback.

 - Status: Core Inference Engine & Traceback Verified.

-----------------------------------
## ◈ Key Features

### 1. High-Performance Banded Alignment
Optimized for long-read sequencing (Nanopore/PacBio), the engine utilizes a banded dynamic programming approach. This restricts the search space around the main diagonal using a Sakoe-Chiba constraint, reducing complexity from $O(NM)$ to $O(N \times \text{bandwidth})$.

### 2. Cross-Platform Optimization
The engine is designed for dual-target hardware optimization:

- Apple Silicon: Optimized for ARM64 using -mcpu=apple-m1 and loop unrolling for M-series performance cores.

- Linux HPC: Targets x86_64 AVX2 and FMA instruction sets for high-throughput execution on Intel Xeon/AMD EPYC clusters.

### 3. MLOps-Ready Orchestration
Designed for integration into Agentic AI and automated R&D pipelines:

- Hybrid Configuration: Sane C++ defaults in hmm_model.hpp with dynamic overrides via config.json.

- Workflow Automation: Execution is managed via input.json, allowing for reproducible pipeline runs without recompilation.

-----------------------------------
## ◈ Architectural Core

### 1. Log-Space Dynamic Programming: 
To eliminate numerical underflow and increase speed, all probabilistic calculations are performed in log-space, replacing expensive floating-point multiplications with high-speed additions.

### 2. Robustness to Biological Noise
- **Affine Gap Scoring**: Correctly models long structural variants using differentiated Gap-Open ($\gamma$) and Gap-Extend ($\epsilon$) penalties: $Score = \gamma + (n-1)\epsilon$.

- **Soft-Clipping**: Employs local-entry and local-exit logic to identify and clip non-matching adapter sequences (e.g., 3S14223M).

- **Glocal Alignment**: Supports global-local transitions, allowing a read to align to a specific sub-region of a larger reference without flank penalties.

-----------------------------------
## ◈ Project Structure

```
HMMAlign/
├── .venv/                # Python virtual environment (managed by uv)
├── src/                  # C++ Source Code
│   ├── viterbi.cpp       # Core Viterbi algorithm implementation
│   └── hmm_model.hpp     # HMM state and transition parameters
├── hmmalign/             # Python Package Source
│   ├── __init__.py       # Package entry point and bridge logic
│   └── _core*.so         # Compiled C++ binary (generated)
├── data/                 # Sample FASTA files
├── results/              # Output (SAM files)
├── pyproject.toml        # Project metadata and build dependencies
├── CMakeLists.txt        # C++ build configuration
├── uv.lock               # Deterministic dependency lockfile
├── input.json            # Pipeline execution instructions
├── config.json           # Bio-parameters (bandwidth, penalties)
├── rebuild.sh            # MLOps runner
└── main.py               # Comprehensive Test Suite & Benchmarking
```

-----------------------------------
## ◈ Getting Started

### Prerequisites
- Python 3.12+
- uv (fast Python package manager)
- C++ Compiler (Clang 15+ or GCC 12+)
- jq (for JSON parsing in the runner)

### Installation & Execution

1. Clone and Build:
```
git clone https://github.com/XuejianXiong/HMMAlign.git
cd HMMAlign
./rebuild.sh

```

2. Configuration:

- input.json: Define file paths and MLOps settings (e.g., rebuild_cpp: true).

- config.json: Tune biological parameters (Match/Mismatch scores, Gap penalties).

-----------------------------------
## ◈ Output & Visualization
The engine generates standard SAM (Sequence Alignment Map) files, making it compatible with industry-standard visualization tools like IGV (Integrative Genomics Viewer).

### Summary Statistics Output:
```
════════════════════════════════════════════════════════════
 HMMAlign Summary Statistics
 ──────────────────────────────────────────────────────────
 Reference:  MT_human (16569 bp)
 Query:      MT_orang (16499 bp)
 Status:     Success
 ──────────────────────────────────────────────────────────
 Identity:    75.84%
 Matches:     14,223 bp
 Indels:      2,273 Ins / 2,258 Del
 Score:       10259.21
 Latency:     889.37 ms
════════════════════════════════════════════════════════════
```

-----------------------------------
## ◈ License

MIT License
