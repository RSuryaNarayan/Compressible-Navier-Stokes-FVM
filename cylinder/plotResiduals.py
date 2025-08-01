import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the CSV file
data100 = pd.read_csv('residuals.csv')
# data200 = pd.read_csv('residuals_200.csv')
# data300 = pd.read_csv('residuals_300.csv')
# data400 = pd.read_csv('residuals_400.csv')
# data600 = pd.read_csv('residuals_600.csv')
# data983 = pd.read_csv('residuals_983.csv')

# Create a log-log plot
plt.figure(figsize=(10, 6))
plt.loglog(data100['Index'], data100['Residual'], 'b', linestyle='-', linewidth=2, label='CFL=7, $dt$=3.35133e-06')
# plt.loglog(data200['Index'], data200['Residual'], 'purple', linestyle='-', linewidth=2, label='CFL=200, $dt$=6.70267e-06')
# plt.loglog(data300['Index'], data300['Residual'], 'r', linestyle='-', linewidth=2, label='CFL=300, $dt$=1.0054e-05')
# plt.loglog(data400['Index'], data400['Residual'], 'g', linestyle='-', linewidth=2, label='CFL=400, $dt$=1.34053e-05')
# plt.loglog(data600['Index'], data600['Residual'], 'orange', linestyle='-', linewidth=2, label='CFL=600, $dt$=2.0108e-05')
# plt.loglog(data983['Index'], data983['Residual'], 'k', linestyle='-', linewidth=5, label='CFL=983, $dt$=3.29436e-05')
plt.legend()

# Set labels and title
plt.xlabel('Time step')
plt.ylabel('Residual')
plt.title('Log-Log Plot of Residuals')
plt.grid(True)
plt.savefig("residuals_2_1.png")

# Show the plot
plt.show()