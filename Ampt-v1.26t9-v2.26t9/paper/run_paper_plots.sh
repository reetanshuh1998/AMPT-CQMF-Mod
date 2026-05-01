#!/bin/bash
echo "============================================="
echo "Generating PRC 96 (2017) Validation Figures"
echo "============================================="

# Ensure we're running inside the paper directory
cd "$(dirname "$0")"

echo "1. Generating Invariant Transverse Mass Spectra (Figure 22 Analog)..."
python3 plot_fig1_mt_spectra.py

echo "2. Generating Mean Transverse Momentum (Figure 18 Analog)..."
python3 plot_fig18_mean_pt.py

echo "3. Generating Particle Yield Ratios (Figure 19 Analog)..."
python3 plot_fig19_particle_ratios.py

echo "============================================="
echo "All validation figures constructed successfully!"
echo "Figures saved locally in paper/ directory."
echo "4. Generating Mixed Particle Yield Ratios (Figure 20 Analog)..."
python3 plot_fig20_ratios.py
