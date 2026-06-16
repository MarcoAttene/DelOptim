#!/usr/bin/env bash
# Build the delmesher Python bindings into a self-contained virtual environment
# under build/venv and run the Python test suite.
#
# This is the recommended way to develop and test the bindings locally; the same
# pyproject.toml drives the wheels published to PyPI (see scripts/build_wheels.sh).
#
# Usage:
#   scripts/build_python_venv.sh           # build + install + run tests
#   scripts/build_python_venv.sh --no-test # build + install only
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

VENV_DIR="$REPO_ROOT/build/venv"
PYTHON="${PYTHON:-python3}"

echo "==> Building the CLI (Release) for the parity tests"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDELMESHER_BUILD_TESTS=OFF >/dev/null
cmake --build build --config Release --target delmesher -j

echo "==> Creating virtual environment at $VENV_DIR"
"$PYTHON" -m venv "$VENV_DIR"
# shellcheck disable=SC1091
VPY="$VENV_DIR/bin/python"
[ -x "$VPY" ] || VPY="$VENV_DIR/Scripts/python.exe"   # Windows layout

"$VPY" -m pip install --upgrade pip >/dev/null

echo "==> Building and installing delmesher into the venv (this compiles the extension)"
"$VPY" -m pip install -v ".[test]"

if [ "${1:-}" = "--no-test" ]; then
    echo "==> Skipping tests (--no-test)"
    exit 0
fi

echo "==> Running the Python test suite"
DELMESHER_BINARY="$REPO_ROOT/build/delmesher" \
    "$VPY" -m pytest tests/python -v
