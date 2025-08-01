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
df = pd.read_csv('solution100000.csv')  # Updated to use csv_file

# Get unique x and y coordinates to define grid lines
x_coords = df['x_c']  # sorted(df['x'].unique())
y_coords = df['y_c']  # sorted(df['y'].unique())
rho = df['Ri']
u = df['u']
v = df['v']
p = df['p']
T = df['T']

x_coords = np.array(x_coords)
y_coords = np.array(y_coords)
rho = np.array(rho)
u = np.array(u)
v = np.array(v)
p = np.array(p)
T = np.array(T)

x_coords = x_coords.reshape((nx, ny))
y_coords = y_coords.reshape((nx, ny))
rho = rho.reshape((nx, ny))
u = u.reshape((nx, ny))
v = v.reshape((nx, ny))
p = p.reshape((nx, ny))
T = T.reshape((nx, ny))

# Plot rho
plt.figure(1)
plt.contourf(x_coords, y_coords, rho, cmap='jet', levels=1000, extend='both')  # Increased levels for smoother contours
plt.colorbar(label='Density residual (rho)')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Density Field ({nx+1} x {ny+1})')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.savefig('rho.png', bbox_inches='tight')
plt.show()

# Plot u
plt.figure(2)
plt.contourf(x_coords, y_coords, u, cmap='jet', levels=1000, extend='both')  # Increased levels for smoother contours
plt.colorbar(label='Velocity (u)')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Velocity Field (u) ({nx+1} x {ny+1})')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.savefig('u.png', bbox_inches='tight')
plt.show()

# Plot v
plt.figure(3)
plt.contourf(x_coords, y_coords, v, cmap='jet', levels=1000, extend='both')  # Increased levels for smoother contours
plt.colorbar(label='Velocity (v)')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Velocity Field (v) ({nx+1} x {ny+1})')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.savefig('v.png', bbox_inches='tight')
plt.show()

# Plot p
plt.figure(4)
plt.contourf(x_coords, y_coords, p, cmap='jet', levels=1000, extend='both')  # Increased levels for smoother contours
plt.colorbar(label='Pressure (p)')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Pressure Field ({nx+1} x {ny+1})')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.savefig('p.png', bbox_inches='tight')
plt.show()

# Plot T
plt.figure(5)
plt.contourf(x_coords, y_coords, T, cmap='jet', levels=500,extend='both')  # Increased levels for smoother contours
plt.colorbar(label='Temperature (T)')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Temperature Field ({nx+1} x {ny+1})')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.savefig('T.png', bbox_inches='tight')
plt.show()



