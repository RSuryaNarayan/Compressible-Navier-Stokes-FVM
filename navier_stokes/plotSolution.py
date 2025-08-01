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
df = pd.read_csv('solution2100.csv')  # Updated to use csv_file

# Get unique x and y coordinates to define grid lines
x_coords = df['x_c']  # sorted(df['x'].unique())
y_coords = df['y_c']  # sorted(df['y'].unique())
rho = df['rho']
u = df['u']
v = df['v']
p = df['p']
T = df['T']
gradUy = df['gradUy']
gradTy = df['gradTy']

x_coords = np.array(x_coords)
y_coords = np.array(y_coords)
rho = np.array(rho)
u = np.array(u)
v = np.array(v)
p = np.array(p)
T = np.array(T)
gradUy = np.array(gradUy)
gradTy = np.array(gradTy)

x_coords = x_coords.reshape((nx, ny))
y_coords = y_coords.reshape((nx, ny))
rho = rho.reshape((nx, ny))
u = u.reshape((nx, ny))
v = v.reshape((nx, ny))
p = p.reshape((nx, ny))
T = T.reshape((nx, ny))
gradUy = gradUy.reshape((nx, ny))
gradTy = gradTy.reshape((nx, ny))

mu = 1.458 * 1e-6 * T**1.5 / (T + 110.3)
cv = 287/(1.4-1)
cp = cv + 287
kappa = (cp)*mu/0.7
a = np.sqrt(1.4*p/rho)
mag = np.sqrt(u**2 + v**2)
Ma = mag/a
dx = 1.0/(nx)
dy = 0.25/(ny)

Re_cell = rho[:,0]*(np.sqrt(u[:,0]**2+v[:,0]**2))*dx/mu[:,0]
tau_wall = mu[:,0]*gradUy[:,0]
u_tau = np.sqrt(np.abs(tau_wall)/rho[:,0])
y_plus = rho[:,0] * u_tau * dy / mu[:,0]
q_wall = kappa[:,0] * gradTy[:,0]

Pr = 0.7
rho_inf = 0.01
u_inf = 1.5 * np.sqrt(1.4*287*300)
mu_inf = 1.458 * 1e-6 * 300**1.5 / (300 + 110.3)
qsrtx = q_wall * np.sqrt(x_coords[:,0])
tsqrtx = tau_wall * np.sqrt(x_coords[:,0])

chRex = qsrtx / 0.5/Pr**0.5 / rho_inf**0.5/u_inf**2.5/mu_inf**0.5
cfRex = tsqrtx/0.5/rho_inf**0.5/u_inf**1.5/mu_inf**0.5

plt.figure(1)
plt.plot(x_coords[:,0],cfRex, linewidth=2)
plt.title(r"Viscous Flat Plate Problem, $c_f\sqrt{Re_x}$ vs. $x$")
plt.xlabel("$x[m]$")
plt.ylabel(r"$c_f\sqrt{Re_x}$")
plt.savefig("cfRex.png",bbox_inches='tight')
plt.show()

# # Plot rho
# plt.figure(1)
# contour = plt.contourf(x_coords, y_coords, rho, cmap='jet', levels=1000)  # Increased levels for smoother contours
# cbar = plt.colorbar(contour, fraction=0.012, pad=0.04, label='Density residual (rho)')  # Set aspect ratio for colorbar
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Density Residual Field ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('rho_finest.png', bbox_inches='tight')
# plt.show()

# # Plot u
# plt.figure(2)
# contour = plt.contourf(x_coords, y_coords, u, cmap='jet', levels=1000)  # Increased levels for smoother contours
# cbar = plt.colorbar(contour, fraction=0.012, pad=0.04, label='Velocity (u)')  # Set aspect ratio for colorbar
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Velocity Field (u) ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('u_finest.png', bbox_inches='tight')
# plt.show()

# # Plot v
# plt.figure(3)
# contour = plt.contourf(x_coords, y_coords, v, cmap='jet', levels=1000)  # Increased levels for smoother contours
# cbar = plt.colorbar(contour, fraction=0.012, pad=0.04, label='Velocity (v)')  # Set aspect ratio for colorbar
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Velocity Field (v) ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('v_finest.png', bbox_inches='tight')
# plt.show()

# # # Plot p
# plt.figure(4)
# contour = plt.contourf(x_coords, y_coords, p, cmap='jet', levels=1000)  # Increased levels for smoother contours
# cbar = plt.colorbar(contour, fraction=0.012, pad=0.04)  # Set aspect ratio for colorbar
# cbar.set_label('Pressure (p)')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Pressure Field ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('pressure_finest.png', bbox_inches='tight')
# plt.show()

# # # Plot T
# plt.figure(5)
# contour = plt.contourf(x_coords, y_coords, T, cmap='jet', levels=500)  # Increased levels for smoother contours
# cbar = plt.colorbar(contour, fraction=0.012, pad=0.04)  # Set aspect ratio for colorbar
# cbar.set_label('Temperature (T)')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Temperature Field ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('Temp_finest.png', bbox_inches='tight')
# plt.show()

# # Plot Mach Number
# plt.figure(6)
# contour = plt.contourf(x_coords, y_coords, Ma, cmap='jet', levels=500)  # Increased levels for smoother contours
# cbar = plt.colorbar(contour, fraction=0.012, pad=0.04)  # Set aspect ratio for colorbar
# cbar.set_label('Mach Number (Ma)')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title(f'Mach Number Field ({nx+1} x {ny+1})')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(False)  # Disable default grid to match the style
# plt.savefig('Ma_finest.png', bbox_inches='tight')
# plt.show()


# ## centerline values
# # plt.figure(7)
# # fig, ax1 = plt.subplots()  # Create a figure and a set of subplots

# # # Plot rho on the first y-axis
# # ax1.plot(x_coords[:, ny//2], rho[:, ny//2], 'b', linestyle='-', linewidth=2, label=r'$\rho_{centerline,CFD}$')
# # ax1.set_ylabel(r'$\rho$', color='b')  # Set the y-axis label for rho
# # ax1.plot(x_coords[:, ny//2],  0.21678845882238504*np.ones((nx,1)), 'b', linestyle='--', linewidth=1, label=r'$\rho_{\theta-\beta-M}$')
# # ax1.tick_params(axis='y', labelcolor='b')  # Set the color of the y-axis ticks

# # # Create a second y-axis for T and Ma
# # ax2 = ax1.twinx()  
# # ax2.plot(x_coords[:, ny//2], T[:, ny//2], 'r', linestyle='-', linewidth=2, label='$T_{centerline,CFD}$')
# # ax2.plot(x_coords[:, ny//2], 396.5875982399144*np.ones((nx,1)), 'r', linestyle='--', linewidth=1, label=r'$T_{\theta-\beta-M}$')
# # ax2.set_ylabel('T', color='r')  # Set the y-axis label for T
# # ax2.tick_params(axis='y', labelcolor='r')  # Set the color of the y-axis ticks

# # ax3 = ax1.twinx()  # Create a third y-axis
# # ax3.spines['right'].set_position(('outward', 40))  # Offset the third y-axis
# # ax3.plot(x_coords[:, ny//2], Ma[:, ny//2], 'k', linestyle='-', linewidth=2, label='$Ma_{centerline,CFD}$')
# # ax3.plot(x_coords[:, ny//2], 1.873526006750317*np.ones((nx,1)), 'k', linestyle='--', linewidth=1, label=r'$Ma_{\theta-\beta-M}$')
# # ax3.set_ylabel('$Ma$', color='k')  # Set the y-axis label for Ma
# # ax3.tick_params(axis='y', labelcolor='k')  # Set the color of the y-axis ticks

# # # Combine legends from all three axes
# # lines, labels = ax1.get_legend_handles_labels()
# # lines2, labels2 = ax2.get_legend_handles_labels()
# # lines3, labels3 = ax3.get_legend_handles_labels()
# # ax1.legend(lines + lines2 + lines3, labels + labels2 + labels3, loc='center left')
# # plt.title(r"Centerline values of $\rho,T,Ma$ vs. $x$")
# # ax1.set_xlabel(r"$x$ [m]")
# # plt.savefig("validation.png",bbox_inches='tight')
# # plt.show()

# # plt.figure(9)
# # contour = plt.contourf(x_coords, y_coords, gradVy, cmap='jet', levels=500)  # Increased levels for smoother contours
# # cbar = plt.colorbar(contour, fraction=0.012, pad=0.04)  # Set aspect ratio for colorbar
# # cbar.set_label('gradVy (s-1)')
# # plt.xlabel('x')
# # plt.ylabel('y')
# # plt.title(f'gradVy Field ({nx+1} x {ny+1})')
# # plt.gca().set_aspect('equal', adjustable='box')
# # plt.grid(False)  # Disable default grid to match the style
# # plt.savefig('gradVy.png', bbox_inches='tight')
# # plt.show()