
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
dx = rx[1]-rx[0]
dy = ry[1]-ry[0]
dz = rz[1]-rz[0]
dr = [dx, dy, dz]
for i in range(3):
    if dr[i]>L/2:
        print("antes", dr[i])
        dr[i]-=L
        print("dps", dr[i])
    elif dr[i]<-L/2:
        print("antes", dr[i])
        dr[i]+=L
        print("dsp", dr[i])
#dx, dy, dz =  3.1803218040094228E-002,  0.11412444806163841,   0.30798857106300215
print("dr final" + str(dr))
dt = 1
d = np.sqrt((dr[0])**2 + (dr[1])**2 + (dr[2])**2 )
print(d)
f_mag = - 24 * 1.0 * (2*(1.0**12 / d**13) - (1.0**6 / d**7))
#f_mag = 24 * 1 * (-2*(1**12 / d**13) + (1**6 / d**7))
print("f_mag", f_mag)
fx = f_mag * (dx)/d
fy = f_mag * (dy)/d
fz = f_mag * (dz)/d
R = (rx + fx/2,ry + fy/2,rz + fz/2)
print(R)

