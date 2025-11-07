#!/usr/bin/env python3
# plot_trayectoria.py
# Trayectoria 3D de 1 partícula (índice 0) con PBC
# Usa: python plot_trayectoria.py

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
interval_ms = 1

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

files = sorted(glob.glob(os.path.join(estados_dir, "estado_*.dat")))
n_steps = len(files)

# Detectar N_part
sample = np.loadtxt(files[0], comments='#')
N_part = len(sample)

for particula_id in range(N_part):
    # ================================
    # CARGAR TODO EN MEMORIA (rápido)
    # ================================
    print("Cargando trayectorias...")
    files = sorted(glob.glob(os.path.join(estados_dir, "estado_*.dat")))
    n_steps = len(files)

    # Detectar N_part
    sample = np.loadtxt(files[0], comments='#')
    N_part = sample.shape[0]

    # Array: (paso, partícula, xyz)
    pos = np.zeros((n_steps, N_part, 3))

    for i, f in enumerate(files):
        dat = np.loadtxt(f, comments='#')
        pos[i] = dat[:, 1:4]  # rx, ry, rz

    print(f"Listo: {n_steps} pasos, partícula {particula_id}")

    # ================================
    # APLICAR PBC A TRAYECTORIA (continua)
    # ================================
    traj = pos[:, particula_id, :].copy()

    # Desenrollar PBC: evitar saltos bruscos
    for i in range(1, n_steps):
        dr = traj[i] - traj[i-1]
        dr = dr - L * np.round(dr / L)  # imagen mínima
        traj[i] = traj[i-1] + dr

    print("PBC aplicado → trayectoria continua")

    # ================================
    # PLOT 3D
    # ================================


    # Trayectoria
    rx, ry, rz = traj.T
    ax.plot(rx, ry, rz, c='blue', lw=2, label=f'Partícula {particula_id}')

    # Puntos inicial/final
    ax.scatter(rx[0], ry[0], rz[0], c='green', s=100, label='Inicio')
    ax.scatter(rx[-1], ry[-1], rz[-1], c='red', s=100, label='Fin')

    # Caja
    corners = np.array([[0,0,0],[L,0,0],[L,L,0],[0,L,0],[0,0,L],[L,0,L],[L,L,L],[0,L,L]])
    edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    for e in edges:
        ax.plot(*zip(corners[e[0]], corners[e[1]]), c='k', alpha=0.5)

    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title(f'Trayectoria de Partícula {particula_id}\n'
                f'{n_steps} pasos, PBC aplicado')
    #ax.legend()
    ax.grid(True, alpha=0.3)

    # Vista óptima
    ax.view_init(elev=20, azim=45)
plt.tight_layout()
plt.show()