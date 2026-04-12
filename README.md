# 🧬 HMMAlign
## High-Performance Viterbi Engine for Reference-Guided Sequence Alignment

## 🧠 Overview:

HMMAlign is a high-performance computational framework designed for the precise mapping of DNA reads to genomic references. By leveraging a 3-state Hidden Markov Model (HMM) and the Viterbi algorithm, the engine effectively resolves the indel noise and translocation speed variability inherent in modern sequencing technologies.

 - Status: 🛠️ Architecting Core C++ Inference Engine / Bridge Verified

-----------------------------------
## 🚀 The Mission

HMMAlign provides the "algorithmic brain" to map sequencing data back to a known reference by:

1. Modeling Hidden States: Treating genomic positions as hidden states and basecall posteriors as emissions.

2. Global Optimization: Using the Viterbi Algorithm to find the single most likely path through the reference space, effectively "stretching" or "compressing" the read to fit the biological truth.

3. Probabilistic Error Correction: Automatically resolving ambiguities in high-noise regions by utilizing the known reference context to guide path selection.

-----------------------------------
## 🏗️ Architectural Core

### 1. Signal-to-Reference Mapping
Unlike a standard "blind" basecaller, HMMAlign is reference-aware. It utilizes a state-space model where transitions are guided by the expected genomic sequence, significantly increasing alignment sensitivity and consensus accuracy.

### 2. The Viterbi Implementation
The core engine is developed with a focus on low-latency primary analysis:

- Log-Space Dynamic Programming: Eliminates numerical underflow and replaces expensive multiplications with additions to increase precision and speed.

- Cache-Friendly Data Structures: Uses 1D memory layouts for DP tables to maximize CPU cache hits during matrix updates.

- Cache-Friendly Data Structures: Uses 1D memory layouts for DP tables to maximize CPU cache hits during matrix updates.

-----------------------------------
## 🛠️ Tech Stack

- Engine: Modern C++ (C++17/20) for the high-performance DP matrix.

- Build System: scikit-build-core + CMake for seamless C++ compilation within Python environments.

- Environment: uv for lightning-fast, deterministic dependency management.

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
├── check_bridge.py       # Sanity check for C++/Python integration
└── main.py               # Main application entry point
```

-----------------------------------
## 🗺️ Roadmap

- [x] Infrastructure: Established C++/Python bridge and build system.

- [ ] Alpha Engine: Implementation of the core 3-state (M, I, D) Viterbi matrix.

- [ ] Backtrace Logic: Development of the traceback path to generate CIGAR strings.

- [ ] Benchmarking: Comparison against standard alignment tools for accuracy and throughput.

-----------------------------------
## 🔬 Industry Context

In the development of next-generation sequencing platforms, the transition from raw data to a mapped read is a critical computational bottleneck. HMMAlign demonstrates a modular, scalable approach to this problem, separating high-speed sequence logic from upstream feature extraction.

-----------------------------------
## 📘 License

MIT License
