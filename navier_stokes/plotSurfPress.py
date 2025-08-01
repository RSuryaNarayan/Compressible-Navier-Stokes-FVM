import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def get_surf_press(datafile, nx, ny):

    df = pd.read_csv(datafile) 
    x_coords = df['x_c']  
    p = df['p']

    x_coords = np.array(x_coords)
    p = np.array(p)

    x_coords = x_coords.reshape((nx, ny))
    p = p.reshape((nx, ny))

    return x_coords[:,0], p[:,0]


x50, p50 = get_surf_press('solution3000_50.csv',50,25)
x100, p100 = get_surf_press('solution3000_100.csv',100,50)
x200, p200 = get_surf_press('solution3000_200.csv',200,100)
x400, p400 = get_surf_press('solution1000_400_200.csv', 400, 200)
plt.plot(x50, p50, 'r', linestyle='-', linewidth=2, label=r'$G_0(50\times25)$')
plt.plot(x100, p100, 'k', linestyle='-', linewidth=2, label=r'$G_1(100\times50)-baseline$')
plt.plot(x200, p200, 'b', linestyle='-', linewidth=2, label=r'$G_2(200\times100)$')
plt.plot(x400, p400, 'g', linestyle='-', linewidth=2, label=r'$G_3(400\times200)$')
plt.title(r"Surface pressure $(p_{wall})$ vs. $x_{wall}$")
plt.xlabel(r"$x[m]$")
plt.ylabel(r"$p_{wall} [Pa]$")
plt.legend()
plt.savefig("grid_conv.png",bbox_inches='tight')
plt.show()
