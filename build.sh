#!/usr/bin/env bash

set -euo pipefail

# Defaults
BUILD_TYPE="Release"
CLEAN=false
COVERAGE=false

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        debug)
            BUILD_TYPE="Debug"
            ;;
        release)
            BUILD_TYPE="Release"
            ;;
        coverage)
            BUILD_TYPE="Debug"
            COVERAGE=true
            ;;
        clean)
            CLEAN=true
            ;;
        *)
            echo "Unknown argument: $arg"
            echo "Usage: $0 [debug|release|coverage] [clean]"
            exit 1
            ;;
    esac
done

# Detect CPU cores
if command -v nproc >/dev/null 2>&1; then
    THREADS=$(nproc)
elif [[ "$OSTYPE" == "darwin"* ]]; then
    THREADS=$(sysctl -n hw.ncpu)
else
    THREADS=4
fi

echo "Build type: $BUILD_TYPE"
echo "Using $THREADS threads"

BUILD_DIR="build-${BUILD_TYPE,,}"

# Optional clean
if [ "$CLEAN" = true ]; then
    echo "Cleaning $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake -S .. -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DSPIDA_COVERAGE="$( [ "$COVERAGE" = true ] && echo ON || echo OFF )"
cmake --build . --parallel "$THREADS"

