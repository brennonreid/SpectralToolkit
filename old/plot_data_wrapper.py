#!/usr/bin/env python3
"""
plot_data_wrapper.py

Reads the tab-separated Sigma and Energy data output by the C++ engine 
and generates the Zeta Energy Landscape plot using Matplotlib.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def plot_cpp_output(input_file="zeta_landscape_cpp_data.txt", output_file="zeta_landscape_cpp_plot.png"):
    
    if not os.path.exists(input_file):
        print(f"[FATAL] Data file not found: {input_file}")
        print("Please ensure the C++ program was run successfully and that it output the data file.")
        sys.exit(1)

    try:
        # Load data from the text file. 
        # C++ output is tab-separated (\t) and skips the header (skiprows=1).
        data = np.loadtxt(input_file, delimiter='\t', skiprows=1)
        
        # Data structure: data[:, 0] = Sigma, data[:, 1] = Energy
        sigmas = data[:, 0]
        energies_float = data[:, 1]
    
    except Exception as e:
        print(f"[FATAL] Could not load or parse data from {input_file}: {e}")
        return

    # --- 1. Data Analysis (Finding the Lock) ---
    min_energy_float = energies_float.min()
    min_idx = energies_float.argmin()
    min_sigma = sigmas[min_idx]

    # --- 2. PLOTTING ---
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(sigmas, energies_float, linewidth=2, color='darkblue')

    # Critical line marker (The Lock)
    ax.axvline(
        x=0.5,
        linestyle="--",
        color='skyblue',
        alpha=0.8,
        linewidth=1.5,
        label="Critical line (0.5)",
    )

    # Highlight minimum (The Valley Floor)
    ax.scatter(
        [min_sigma],
        [min_energy_float],
        s=150,
        zorder=5,
        edgecolor="black",
        color='red',
        linewidth=2,
        label=f"Minimum σ ≈ {min_sigma:.8f}",
    )

    # Annotation box
    info_text = f"min σ = {min_sigma:.8f}\nE_min = {min_energy_float:.10f}"
    ax.text(
        0.99,
        0.02,
        info_text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="black"),
    )

    ax.set_title("Zeta Energy Landscape (C++ Arb Core)")
    ax.set_xlabel("Sigma (Re(s))")
    ax.set_ylabel("Total magnitude |zeta(s)|")
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')

    # Save the output file
    fig.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"[INFO] Graph saved to: {output_file}")
    print(f"[INFO] True minimum found at sigma = {min_sigma:.8f}")


if __name__ == "__main__":
    input_file_name = "zeta_landscape_cpp_data.txt"
    output_file_name = "zeta_landscape_cpp_plot.png"

    # Check for input file name (argument 1)
    if len(sys.argv) > 1:
        input_file_name = sys.argv[1]
    
    # Check for output file name (argument 2)
    if len(sys.argv) > 2:
        output_file_name = sys.argv[2]
        
    plot_cpp_output(input_file=input_file_name, output_file=output_file_name)