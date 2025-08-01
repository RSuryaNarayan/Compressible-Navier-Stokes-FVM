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
df_cells = pd.read_csv(csv_file.replace('.csv', '_cells.csv'))  # Updated to use csv_file
df_nodes = pd.read_csv(csv_file)  # Use csv_file

# Get unique x_c and y_c coordinates to define cell centers
x_c = df_cells['x_c'] 
y_c = df_cells['y_c'] 

x_c = np.array(x_c)
y_c = np.array(y_c)

x_c = x_c.reshape((nx-1,ny-1))
y_c = y_c.reshape((nx-1,ny-1))

# Get unique x_c and y_c coordinates to define cell centers
x_n = df_nodes['x'] 
y_n = df_nodes['y'] 

x_n = np.array(x_n)
y_n = np.array(y_n)

x_n = x_n.reshape((nx,ny))
y_n = y_n.reshape((nx,ny))

# Plot the grid lines for the mesh
plt.figure(figsize=(10, 6))

#plot grid
for i in range(nx):
    plt.plot(x_n[i,:],y_n[i,:], color='black', linewidth=1.0)

for i in range(ny):
    plt.plot(x_n[:,i],y_n[:,i], color='black', linewidth=1.0)

plt.scatter(x_c, y_c, s=0.5, c='r',label='cell centers')
plt.legend()
# Set the plot limits and labels
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'{nx} x {ny}')

# Aspect ratio and grid styling to match the example
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.savefig('cylinder_3.png',bbox_inches='tight')
plt.show()
