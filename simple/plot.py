import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Lee datos
data = np.loadtxt('forces.dat', skiprows=1)  # Salta header
i, rx, ry, rz, fx, fy, fz = data.T

# Plot 3D: posiciones como puntos, fuerzas como flechas
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rx, ry, rz, c='b', marker='o', label='Posiciones')
ax.quiver(rx, ry, rz, fx, fy, fz, color='r', length=0.1, normalize=True, label='Fuerzas')  # length ajusta tamaño flechas

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Fuerzas en partículas (Lennard-Jones)')
ax.legend()
plt.show()