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
    nx = int(lines[0].strip().split('=')[1])  # Extracting the value after 'nx='
    ny = int(lines[1].strip().split('=')[1])  # Extracting the value after 'ny='
    csv_file = lines[-1].strip().split('=')[1].strip()  # Read the last line as csv_file and remove leading/trailing spaces

# Load the CSV data
df = pd.read_csv('cylinder.csv')

# Get unique x and y coordinates to define grid lines
x_coords = df['x'] 
y_coords = df['y'] 
index = df['index']

x_coords = np.array(x_coords)
y_coords = np.array(y_coords)
index = np.array(index)

x_coords = x_coords.reshape((nx,ny))
y_coords = y_coords.reshape((nx,ny))
index = index.reshape((nx,ny))

# Plot the grid points to form a mesh
plt.figure(figsize=(8, 8))

#plot grid
for i in range(nx):
    plt.plot(x_coords[i,:],y_coords[i,:], color='black', linewidth=1.0)

for i in range(ny):
    plt.plot(x_coords[:,i],y_coords[:,i], color='black', linewidth=1.0)

plt.scatter(x_coords, y_coords, color='r')

# Annotate each point with the corresponding index value
for i in range(nx):
    for j in range(ny):
        plt.annotate(index[i, j], (x_coords[i, j]*1.03, y_coords[i, j]*1.03), fontsize=8, ha='right')

# Set equal aspect ratio for a circular appearance
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Grid around cylinder')
plt.grid(True, linestyle='--', alpha=0.5)
# plt.savefig('cylinder_1.pdf',bbox_inches='tight')
plt.show()
