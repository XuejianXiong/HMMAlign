#!/bin/bash
set -e

INPUT_FILE="input.json"
CONFIG_FILE="config.json"

# Extract execution settings
REBUILD=$(jq -r '.execution.rebuild_cpp' $INPUT_FILE)
PY_SCRIPT=$(jq -r '.execution.python_script' $INPUT_FILE)

echo "🛠️  MLOps Pipeline Init"

if [ "$REBUILD" == "true" ]; then
    OS="$(uname -s)"
    
    # 1. Clean the environment (Crucial step!)
    echo "🧹 Cleaning old artifacts..."
    rm -rf build/
    find . -name "*.so" -delete
    
    # 2. Set high-performance flags
    if [ "$OS" == "Darwin" ]; then
        echo "🍎 Targeting Apple Silicon..."
        export CXXFLAGS="-O3 -mcpu=apple-m1 -ffast-math"
    else
        echo "🐧 Targeting Linux AVX2..."
        export CXXFLAGS="-O3 -mavx2 -mfma -ffast-math"
    fi

    # 3. Force Rebuild (The "In-Place" fix)
    echo "🚀 Recompiling Engine..."
    # We can use --no-cache and --force-reinstall to ensure uv doesn't skip the build
    uv pip install -e . --config-settings=editable.mode=inplace
else
    echo "⏩ Using existing binary."
fi

echo "🧪 Executing $PY_SCRIPT..."
uv run python $PY_SCRIPT --config $CONFIG_FILE
