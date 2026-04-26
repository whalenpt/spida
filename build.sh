#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

BUILD_TYPE="Release"
CLEAN=false
COVERAGE=false
TESTS=false

for arg in "$@"; do
    case "$arg" in
        debug) BUILD_TYPE="Debug" ;;
        release) BUILD_TYPE="Release" ;;
        coverage)
            BUILD_TYPE="RelWithDebInfo"
            COVERAGE=true
            TESTS=true
            ;;
        clean) CLEAN=true ;;
        *)
            echo "Unknown argument: $arg"
            exit 1
            ;;
    esac
done

THREADS="${THREADS:-$(nproc)}"

BUILD_DIR="${BUILD_DIR:-build/$(echo "$BUILD_TYPE" | tr '[:upper:]' '[:lower:]')}"

echo "Build dir: $BUILD_DIR"

if [ "$CLEAN" = true ]; then
    rm -rf "$BUILD_DIR"
fi

CMAKE_ARGS=(
    "-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
    "-DSPIDA_COVERAGE=$([ "$COVERAGE" = true ] && echo ON || echo OFF)"
    "-DSPIDA_TEST=$([ "$TESTS" = true ] && echo ON || echo OFF)"
)

cmake -S "$PROJECT_ROOT" -B "$BUILD_DIR" "${CMAKE_ARGS[@]}"
cmake --build "$BUILD_DIR" --parallel "$THREADS"

if [ "$TESTS" = true ]; then
    ctest --test-dir "$BUILD_DIR" --output-on-failure
fi

# Coverage step (only in coverage mode)
if [ "$COVERAGE" = true ]; then
    echo "=============================="
    echo "Generating coverage"
    echo "=============================="

    lcov --capture \
        --directory "$BUILD_DIR" \
        --output-file "$BUILD_DIR/coverage.info" \
        --ignore-errors mismatch

    lcov --remove \
        "$BUILD_DIR/coverage.info" \
        '/usr/*' \
        '*/external/*' \
        --output-file "$BUILD_DIR/coverage.filtered.info"

    echo "Coverage file:"
    echo "$BUILD_DIR/coverage.filtered.info"
fi
