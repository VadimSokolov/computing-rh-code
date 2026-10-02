#!/bin/sh
# link the zero lists and data files into code/ so every script runs from there
cd "$(dirname "$0")/code"
for f in ../zeros/*.npy ../data/*.json; do ln -sf "$f" .; done
echo "ready: cd code && python3 riesz_large.py 1e10"
