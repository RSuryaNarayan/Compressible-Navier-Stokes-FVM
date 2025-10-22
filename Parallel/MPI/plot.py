import numpy as np
import matplotlib.pyplot as plt
import glob
import re

# Find only files that match the pattern solution_t<number>.dat
files = sorted(glob.glob("solution_t*.dat"))

# Extract time from filename using regex (robust)
def extract_time(fname):
    match = re.search(r"solution_t([0-9.]+)\.dat", fname)
    return float(match.group(1)) if match else None

# Filter and sort files by time
files = [(extract_time(f), f) for f in files if extract_time(f) is not None]
files.sort()

plt.figure(figsize=(8,5))

for t, f in files:
    x, rho, u, p = np.loadtxt(f, unpack=True, comments='#')
    plt.plot(x, rho, label=f't={t:.2f}')

plt.xlabel('x')
plt.ylabel('Density')
plt.legend()
plt.title('1D Euler – Sod Shock Tube (Density Evolution)')
plt.tight_layout()
plt.savefig('sod_shock_tube_density_evolution.png', dpi=300)
plt.show()
