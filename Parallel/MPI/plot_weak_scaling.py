#!/usr/bin/env python3
"""
This script reads 'weak_scaling_results.dat' and generates a
weak scaling efficiency plot using Matplotlib and Pandas.
"""

import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
RESULTS_FILE = "weak_scaling_results.dat"
PLOT_FILE = "weak_scaling_plot.png"

# --- 1. Read Data ---
# Use pandas to read the space-separated file, skipping the header line
try:
    data = pd.read_csv(
        RESULTS_FILE,
        sep=r'\s+', # Use regex for one or more spaces (replaces deprecated delim_whitespace)
        comment='#',
        header=None,
        names=['Procs', 'N_total', 'TotalTime', 'Efficiency'] # 4 columns
    )
except FileNotFoundError:
    print(f"Error: Results file '{RESULTS_FILE}' not found.")
    print("Please run './run_scaling.sh' first.")
    exit(1)
except pd.errors.EmptyDataError:
    print(f"Error: Results file '{RESULTS_FILE}' is empty.")
    print("There may have been an error during the benchmark.")
    exit(1)

if data.empty or len(data) < 1:
    print(f"Error: No data found in '{RESULTS_FILE}'.")
    exit(1)

# --- 2. Get Base Parameters for Labels ---
try:
    # Get the base processor count for the Y-axis label
    base_procs = data['Procs'].iloc[0]
    # Calculate the work-per-core for the title
    work_per_core = (data['N_total'] / data['Procs']).mean()
except IndexError:
    print(f"Error: Could not parse data from '{RESULTS_FILE}'.")
    exit(1)


# --- 3. Create Plot ---
plt.figure(figsize=(10, 6))

# Plot Measured Efficiency
plt.plot(
    data['Procs'],
    data['Efficiency'], # Plot the 'Efficiency' column
    marker='o',
    linestyle='-',
    linewidth=2,
    label='Measured Efficiency'
)

# Plot Ideal Efficiency (a flat line at 1.0)
plt.axhline(
    y=1.0,
    linestyle='--',
    color='black',
    label='Ideal Efficiency (100%)'
)

# --- 4. Customize Plot ---
# Title is now dynamic based on the data
plt.title(f'MPI Solver Weak Scaling (Work/Core ≈ {int(work_per_core)} cells)', fontsize=16)
plt.xlabel('Number of Processes', fontsize=12)
# Y-label is now "Efficiency"
plt.ylabel(f'Efficiency (T_{int(base_procs)} / T_P)', fontsize=12)

# Set x-axis to log scale, base 2
plt.xscale('log', base=2)

# Set explicit x-ticks to match the processor counts
plt.xticks(data['Procs'], labels=data['Procs'].astype(str))
# Turn off minor ticks for a cleaner look
plt.minorticks_off() 

# Set y-axis limits to make the efficiency drop-off clear
plt.ylim(0, 1.2) 

plt.grid(True, which='both', linestyle=':', linewidth=0.6)
plt.legend(fontsize=12, loc='lower left')

# --- 5. Save and Show ---
plt.savefig(PLOT_FILE)
print(f"Plot saved to {PLOT_FILE}")

# You can uncomment the line below if you want the script to
# open the plot in a window automatically.
# plt.show()