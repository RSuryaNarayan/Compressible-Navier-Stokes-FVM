import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# grid info
nx = 80
ny = 80
# Load the CSV data
df = pd.read_csv('cylinder.csv')

# Get unique x and y coordinates to define grid lines
x_coords = df['x'] 
y_coords = df['y'] 

x_coords = np.array(x_coords)
y_coords = np.array(y_coords)

x_coords = x_coords.reshape((nx,ny))
y_coords = y_coords.reshape((nx,ny))

# Plot the grid points to form a mesh
plt.figure(figsize=(8, 8))

#plot grid
for i in range(nx):
    plt.plot(x_coords[i,:],y_coords[i,:], color='black', linewidth=1.0)

for i in range(ny):
    plt.plot(x_coords[:,i],y_coords[:,i], color='black', linewidth=1.0)

# Set equal aspect ratio for a circular appearance
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Grid around cylinder')
plt.grid(True, linestyle='--', alpha=0.5)
plt.savefig('cylinder_1.pdf',bbox_inches='tight')
plt.show()
