import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Process input files.')
parser.add_argument('inp_file', type=str, help='Input ramp file (e.g., inp.ramp)')
args = parser.parse_args()  # Parse the command line arguments

# Load grid size from inp.2d file
with open(args.inp_file, 'r') as f:  # Use args.inp_file
    lines = f.readlines()
    nx = int(lines[0].strip().split('=')[1])-1  # Extracting the value after 'nx='
    ny = int(lines[1].strip().split('=')[1])-1  # Extracting the value after 'ny='
    csv_file = lines[-1].strip().split('=')[1].strip()  # Read the last line as csv_file and remove leading/trailing spaces

# Load the CSV data
df = pd.read_csv(csv_file.replace('.csv', '_cells.csv'))  # Updated to use csv_file

# Get unique x and y coordinates to define grid lines
x_coords = df['x_c'] #sorted(df['x'].unique())
y_coords = df['y_c']  #sorted(df['y'].unique())
vol = df['volume'] 

x_coords = np.array(x_coords)
y_coords = np.array(y_coords)
vol = np.array(vol)

x_coords = x_coords.reshape((nx,ny))
y_coords = y_coords.reshape((nx,ny))
vol = vol.reshape((nx,ny))

# Plot 
plt.figure(figsize=(6, 6))
plt.contourf(x_coords, y_coords, vol, cmap='jet', levels=1000,extend='both')  # Increased levels for smoother contours
plt.colorbar(label='Volume')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Cell Volumes ({nx+1} x {ny+1})')
# Aspect ratio and grid styling to match the example
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.savefig('cyl_vol_3.png',bbox_inches='tight')
plt.show()
