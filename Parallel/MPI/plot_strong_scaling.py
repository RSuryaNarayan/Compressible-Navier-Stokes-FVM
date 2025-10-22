#!/usr/bin/env python3
"""
This script reads 'strong_scaling_results.dat' and generates a
strong scaling speedup plot using Matplotlib and Pandas.
"""

import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
RESULTS_FILE = "strong_scaling_results.dat"
PLOT_FILE = "strong_scaling_plot.png"

# --- 1. Read Data ---
# Use pandas to read the space-separated file, skipping the header line
try:
    data = pd.read_csv(
        RESULTS_FILE,
        sep=r'\s+', # Use regex for one or more spaces
        comment='#',
        header=None,
        names=['Procs', 'TotalTime', 'Speedup'] # 3 columns
    )
except FileNotFoundError:
    print(f"Error: Results file '{RESULTS_FILE}' not found.")
    print("Please run './run_strong_scaling.sh' first.")
    exit(1)
except pd.errors.EmptyDataError:
    print(f"Error: Results file '{RESULTS_FILE}' is empty.")
    print("There may have been an error during the benchmark.")
    exit(1)

if data.empty or len(data) < 1:
    print(f"Error: No data found in '{RESULTS_FILE}'.")
    exit(1)

# --- 2. Calculate Ideal Speedup ---
try:
    # Get base procs from the first line of the data
    BASE_PROCS = data['Procs'].iloc[0]
    # Ideal Speedup S(P) = P / P_base
    data['IdealSpeedup'] = data['Procs'] / BASE_PROCS
except IndexError:
    print(f"Error: Could not parse data from '{RESULTS_FILE}'.")
    exit(1)

# --- 3. Create Plot ---
plt.figure(figsize=(10, 6))

# Plot Measured Speedup
plt.plot(
    data['Procs'],
    data['Speedup'],
    marker='o',
    linestyle='-',
    linewidth=2,
    label='Measured Speedup'
)

# Plot Ideal Speedup
plt.plot(
    data['Procs'],
    data['IdealSpeedup'],
    linestyle='--',
    color='black',
    label=f'Ideal Speedup (S = P/{int(BASE_PROCS)})'
)

# --- 4. Customize Plot ---
plt.title(f'MPI Solver Strong Scaling (N_total is fixed)', fontsize=16)
plt.xlabel('Number of Processes', fontsize=12)
plt.ylabel(f'Speedup (Relative to {int(BASE_PROCS)} Processes)', fontsize=12)

# Set x-axis to log scale, base 2
plt.xscale('log', base=2)

# Set explicit x-ticks to match the processor counts
plt.xticks(data['Procs'], labels=data['Procs'].astype(str))
# Turn off minor ticks for a cleaner look
plt.minorticks_off() 

plt.grid(True, which='both', linestyle=':', linewidth=0.6)
plt.legend(fontsize=12, loc='upper left')

# --- 5. Save and Show ---
plt.savefig(PLOT_FILE)
print(f"Plot saved to {PLOT_FILE}")

# You can uncomment the line below if you want the script to
# open the plot in a window automatically.
# plt.show()
