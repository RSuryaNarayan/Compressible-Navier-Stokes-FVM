import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib.patches as patches  # Import patches for drawing shapes

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
df = pd.read_csv('/Users/surya/Desktop/solution1240000.csv')  # Updated to use csv_file

# Get unique x and y coordinates to define grid lines
x_coords = df['x_c']  # sorted(df['x'].unique())
y_coords = df['y_c']  # sorted(df['y'].unique())
rho = df['rho']
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
Ma = np.sqrt(u**2 + v**2)/np.sqrt(1.4*p/rho)

# Plot rho
# plt.figure(1)
# contour = plt.contourf(x_coords, y_coords, Ma, cmap='RdBu', levels=1000)  # Increased levels for smoother contours
# contour.set_clim(0, 5)  # Set color limits for the contour
# cbar = plt.colorbar(contour, fraction=0.15, pad=0.02, label='Mach Number')  # Set aspect ratio for colorbar
# cbar.set_ticks(np.linspace(0.015,5,10))  # Set ticks from 0 to 5
# # plt.colorbar(label='Mach Number')
# plt.xlim([-0.25,0.0])
# plt.ylim([0.0, 0.35])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Mach number Field ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('Ma.png', bbox_inches='tight')
# plt.show()

# Plot u
# plt.figure(2)
# plt.contourf(x_coords, y_coords, u, cmap='jet', levels=1000, extend='both')  # Increased levels for smoother contours
# plt.colorbar(label='Velocity (u)')
# plt.xlim([-0.25,0.0])
# plt.ylim([0.0, 0.35])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Velocity Field (u) ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('u.png', bbox_inches='tight')
# plt.show()

# Plot v
# plt.figure(3)
# plt.contourf(x_coords, y_coords, v, cmap='jet', levels=1000, extend='both')  # Increased levels for smoother contours
# plt.colorbar(label='Velocity (v)')
# plt.xlim([-0.25,0.0])
# plt.ylim([0.0, 0.35])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Velocity Field (v) ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('v.png', bbox_inches='tight')
# plt.show()

# Plot velocity vectors
# plt.figure(3)
# plt.quiver(x_coords, y_coords, u, v, color='k', scale=30000, width=0.005)  # Increased levels for smoother contours
# # plt.colorbar(label='Velocity (v)')
# plt.xlim([-0.25,0.0])
# plt.ylim([0.0, 0.35])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Velocity Vector Field (v) ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('velocity_vector.png', bbox_inches='tight')
# plt.show()

# Plot p
# plt.figure(4)
# contour = plt.contourf(x_coords, y_coords, p, cmap='jet', levels=1000)  # Increased levels for smoother contours
# # contour.set_clim(300, 1700)  # Set color limits for the contour
# cbar = plt.colorbar(contour, fraction=0.15, pad=0.02, label='Pressure [Pa]')  # Set aspect ratio for colorbar
# # cbar.set_ticks(np.linspace(300,1750,10))  # Set ticks from 300 to 1750
# plt.xlim([-0.25,0.0])
# plt.ylim([0.0, 0.35])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Pressure Field ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('p.png', bbox_inches='tight')
# plt.show()

# Plot T
# plt.figure(5)
# contour = plt.contourf(x_coords, y_coords, T, cmap='hot', levels=1000)  # Increased levels for smoother contours
# contour.set_clim(300, 1700)  # Set color limits for the contour
# cbar = plt.colorbar(contour, fraction=0.15, pad=0.02, label='Temperature [K]')  # Set aspect ratio for colorbar
# cbar.set_ticks(np.linspace(300,1750,10))  # Set ticks from 300 to 1750

# # Set axis limits before adding the patch
# plt.xlim([-0.25, 0.0])
# plt.ylim([0.0, 0.35])

# # Add quarter circle patch with adjusted z-order
# quarter_circle = patches.Wedge(center=(0.0, 0.0), r=0.1, theta1=0, theta2=-90, color='black', zorder=5)  # Create a quarter circle with zorder
# plt.gca().add_patch(quarter_circle)  # Add the patch to the current axes

# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Temperature Field ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('T.png', bbox_inches='tight')
# plt.show()

# for i in range(nx):
#     for j in range(ny):
#         if (T[i,j]<300):
#             print(x_coords[i,j], y_coords[i,j], T[i,j])


plt.plot(x_coords[:,0], T[:,0], 'k', linewidth=2.0,label='Centerline temperature profile')
plt.axvline(x=-0.04653-0.1, color='r', linestyle='dashdot',label='Analytically Predicted')
plt.axvline(x=-0.06-0.1, color='b', linestyle='dashdot',label='CFD Predicted')
plt.text(-0.04653-0.1, plt.ylim()[1]*0.9, f'x = {-0.04653-0.1} m', color='r', fontsize=10, ha='right')  # Add x-value vertically
plt.text(-0.06-0.1, plt.ylim()[1]*0.7, f'x = {-0.06-0.1} m', color='b', fontsize=10, ha='right')  # Add x-value vertically
plt.legend()
plt.title("Shock stand-off distance comparison")
plt.xlabel("Centerline distance [m]")
plt.ylabel("Temperature [K]")
plt.savefig("/Users/surya/Desktop/shock_standoff.png")
plt.show()