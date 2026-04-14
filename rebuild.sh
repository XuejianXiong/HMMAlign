#!/bin/bash

# ==============================================================================
# HMMAlign Pipeline Orchestrator
# 
# Description:
#   Automates the compilation of the C++ HMM kernel and executes the 
#   Python-based alignment pipeline. Uses 'jq' for configuration parsing
#   and 'uv' for high-speed dependency and build management.
#
# Usage:
#   ./rebuild.sh
#
# Requirements:
#   - jq (JSON processor)
#   - uv (Python package installer and runner)
#   - C++ compiler (Clang/GCC)
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

# Configuration Constants
INPUT_FILE="input.json"
CONFIG_FILE="config.json"

# Check for required tools
command -v jq >/dev/null 2>&1 || { echo >&2 "❌ Error: 'jq' is not installed. Aborting."; exit 1; }
command -v uv >/dev/null 2>&1 || { echo >&2 "❌ Error: 'uv' is not installed. Aborting."; exit 1; }

# 1. Parse Execution Settings
if [ ! -f "$INPUT_FILE" ]; then
    echo "❌ Error: Configuration file '$INPUT_FILE' not found."
    exit 1
fi

REBUILD=$(jq -r '.execution.rebuild_cpp // "false"' "$INPUT_FILE")
PY_SCRIPT=$(jq -r '.execution.python_script // "main.py"' "$INPUT_FILE")

echo "════════════════════════════════════════════════════════════"
echo "🛠️  MLOps Pipeline: HMMAlign Initialization"
echo "════════════════════════════════════════════════════════════"

# 2. Build Logic
if [ "$REBUILD" == "true" ]; then
    OS="$(uname -s)"
    echo "🔄 Rebuild triggered by configuration..."
    
    # Clean old artifacts to prevent stale binary links
    echo "🧹 Cleaning legacy build artifacts and shared objects..."
    rm -rf build/
    find . -maxdepth 3 -name "*.so" -delete
    
    # Optimization Flag Selection
    if [ "$OS" == "Darwin" ]; then
        echo "🍎 System: macOS detected. Optimizing for Apple Silicon..."
        # Using -mcpu=native for modern M-series chips
        export CXXFLAGS="-O3 -mcpu=native -ffast-math -std=c++17"
    else
        echo "🐧 System: Linux detected. Optimizing for AVX2/FMA..."
        export CXXFLAGS="-O3 -mavx2 -mfma -ffast-math -std=c++17"
    fi

    # Execute In-Place Editable Install
    # The --config-settings ensure the .so is placed correctly for the Python package
    echo "🚀 Recompiling C++ Engine via 'uv'..."
    uv pip install -e . --config-settings=editable.mode=inplace
    echo "✅ Compilation complete."
else
    echo "⏩ Rebuild skipped. Using existing binaries."
fi

# 3. Pipeline Execution
if [ ! -f "$PY_SCRIPT" ]; then
    echo "❌ Error: Python script '$PY_SCRIPT' not found."
    exit 1
fi

echo "🧪 Running Pipeline: $PY_SCRIPT"
echo "📂 Input:  $INPUT_FILE"
echo "⚙️  Config: $CONFIG_FILE"
echo "─"

# Execute the main logic
uv run python "$PY_SCRIPT" --input "$INPUT_FILE" --config "$CONFIG_FILE"

echo "✅ Pipeline execution finished successfully."
echo "════════════════════════════════════════════════════════════"