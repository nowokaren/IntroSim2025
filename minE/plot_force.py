
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import glob

# ================================
# CONFIGURACIÓN
# ================================
estados_dir = "estados"
L = 1.0
particula_id = 0  # ← cambia para otra partícula


files = sorted(os.listdir(estados_dir))
r = []; f = []
for i in files:
    sample = np.loadtxt("estados/" + i, comments='#').T
    rx = sample[1]
    ry = sample[2]
    rz = sample[3]
    fx = sample[4]
    fy = sample[5]
    fz = sample[6]
    r.append(np.sqrt((rx[1]-rx[0])**2 + (ry[1]-ry[0])**2 + (rz[1]-rz[0])**2 ))
    f.append(np.sqrt((fx[1]-fx[0])**2 + (fy[1]-fy[0])**2 + (fz[1]-fz[0])**2 ))
    break
#plt.plot(f, ".-", alpha = 0.5)
#plt.plot(r, ".-", alpha = 0.5, label="distancia")

#plt.legend()
#plt.show()

#plt.figure()
#energy = np.loadtxt("energy.dat", comments='#').T
#plt.plot(r[:len(energy[1])], energy[1])
#plt.xlabel("distancia")
#plt.show()
print(i)
print(rx)
print(ry)
print(rz)
dt = 1
d = np.sqrt((rx[1]-rx[0])**2 + (ry[1]-ry[0])**2 + (rz[1]-rz[0])**2 )
#print()
f_mag = 24 * 1 * (-2*(1**12 / d**13) + (1**6 / d**7))
fx = f_mag * (rx[1]-rx[0])/d
fy = f_mag * (ry[1]-ry[0])/d
fz = f_mag * (rz[1]-rz[0])/d
R = (rx + fx/2,ry + fy/2,rz + fz/2)
print(R)

