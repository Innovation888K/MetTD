#!/bin/sh
# Build the reference kernel from source and run the unit tests.
set -e
cd "$(dirname "$0")/.."
gcc -O2 -std=c11 -o mettd_core src/mettd_core.c -lm
python tests/test_kernel.py
