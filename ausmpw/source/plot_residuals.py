import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('residuals.txt', delimiter=',')

iters = data[:, 0]
res_rho = data[:, 1]
res_rhou = data[:, 2]
res_rhov = data[:, 3]
res_Et = data[:, 4]

plt.semilogy(iters, res_rho, label=r'$\rho$')
plt.semilogy(iters, res_rhou, label=r'$\rho u$')
plt.semilogy(iters, res_rhov, label=r'$\rho v$')
plt.semilogy(iters, res_Et, label=r'$E_t$')

plt.title('Residuals vs. Iterations')
plt.xlabel('Iterations')
plt.legend()
plt.rcParams.update({'font.size': 20})
plt.savefig("residuals.png")