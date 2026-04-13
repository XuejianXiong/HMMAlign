# 🧬 HMMAlign
## High-Performance Viterbi Engine for Reference-Guided Sequence Alignment

## 🧠 Overview:

HMMAlign is a high-performance computational framework designed for the precise mapping of DNA reads to genomic references. By leveraging a 3-state Hidden Markov Model (HMM) and the Viterbi algorithm, the engine effectively resolves the indel noise and translocation speed variability inherent in modern sequencing technologies.

 - Performance: 0.03ms for 10kb reads on Apple Silicon.

 - Architecture: C++20 Core with Python-driven MLOps orchestration.

 - Status: Core Inference Engine & Traceback Verified.

-----------------------------------
## 🚀 Key Features

### 1. High-Performance Banded Alignment
Optimized for long-read sequencing (Nanopore/PacBio), the engine utilizes a banded dynamic programming approach. This restricts the search space around the main diagonal, reducing complexity from $O(NM)$ to $O(N \times \text{bandwidth})$, enabling sub-millisecond alignment for 10kb+ sequences.

### 2. Cross-Platform High-Performance Kernels
The engine is designed for dual-target optimization:

- Apple Silicon (Local Dev): Optimized for ARM64 using -mcpu=apple-m1 and loop unrolling for M-series performance cores.

- Linux HPC (Production): Targets x86_64 AVX2 and FMA instruction sets for high-throughput execution on Intel Xeon/AMD EPYC clusters.

### 3. MLOps-Ready Orchestration
Designed for integration into Agentic AI and automated R&D pipelines:
- Hybrid Configuration: Sane C++ defaults in hmm_model.hpp with dynamic overrides via config.json.
- Workflow Automation: Execution is managed via input.json, allowing for reproducible pipeline runs without recompilation.

-----------------------------------
## 🏗️ Architectural Core

### 1. Log-Space Dynamic Programming: 
To eliminate numerical underflow and increase speed, all probabilistic calculations are performed in log-space, replacing expensive floating-point multiplications with high-speed additions.

### 2. Robustness to Biological Noise
The engine is rigorously tested against:

- **N-Neutrality**: Efficiently bridges regions of uncertainty (N-islands) without path fragmentation.

- **Affine Gap Scoring**: Correctly models long structural variants (50bp+ gaps) using differentiated Gap-Open and Gap-Extend penalties.

- **Soft-Clipping** (Terminal Noise): Employs local-entry and local-exit logic to identify and clip non-matching adapter sequences or low-quality terminal bases (e.g., 15S85M).

- **Glocal Alignment**: Supports global-local transitions, allowing the read to align to a specific sub-region of a larger reference without incurring heavy global-alignment penalties at the flanks.

-----------------------------------
## 🛠️ Tech Stack

- Engine: Modern C++20 (Multi-target optimization: ARM64 & AVX2).

- Build System: scikit-build-core + CMake for seamless C++ compilation within Python environments.

- HPC Integration: Built to run in containerized environments (Docker/Singularity) and compatible with Slurm-based task arrays.

- Environment: uv for deterministic, lightning-fast dependency management across local and remote nodes.

- Interface: Python bindings for flexible research iteration and testing.

-----------------------------------
## 📂 Project Structure

```
HMMAlign/
├── .venv/                # Python virtual environment (managed by uv)
├── src/                  # C++ Source Code
│   ├── viterbi.cpp       # Core Viterbi algorithm implementation
│   └── hmm_model.hpp     # HMM state and transition parameters
├── hmmalign/             # Python Package Source
│   ├── __init__.py       # Package entry point and bridge logic
│   └── _core*.so         # Compiled C++ binary (generated)
├── pyproject.toml        # Project metadata and build dependencies
├── CMakeLists.txt        # C++ build configuration
├── uv.lock               # Deterministic dependency lockfile
├── input.json            # Pipeline execution instructions
├── config.json           # Bio-parameters (bandwidth, penalties)
├── rebuild.sh            # MLOps runner
└── main.py               # Comprehensive Test Suite & Benchmarking
```

-----------------------------------
## 🗺️ Roadmap

- [x] Infrastructure: Established C++/Python bridge and build system.

- [x] Alpha Engine: Implementation of the core 3-state (M, I, D) Viterbi matrix.

- [x] Backtrace Logic: Development of the traceback path to generate CIGAR strings.

- [ ] Benchmarking: Comparison against standard alignment tools for accuracy and throughput.

- [ ] Feature: Support for multi-threaded batch alignment (SIMD).

- [ ] Integration: Deployment as a tool for Agentic AI genomic workflows.

-----------------------------------
## 🔬 Industry Context

In the development of next-generation sequencing platforms, the transition from raw data to a mapped read is a critical computational bottleneck. HMMAlign demonstrates a modular, scalable approach to this problem, separating high-speed sequence logic from upstream feature extraction.

-----------------------------------
## 📘 License

MIT License
