import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('wallPressure.txt', delimiter=',')
x = data[:,0]
Pwall = data[:,1]

plt.plot(x, Pwall, linewidth=3)

plt.title(r'Wall Pressure vs. $x$')
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$p$')
plt.rcParams.update({'font.size': 20})
plt.savefig("wallPressure.png")