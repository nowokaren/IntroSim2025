#!/usr/bin/env python3
# plot_E_3D_fast.py
# 10 000 pasos → < 5 s, < 30 % CPU
# -------------------------------------------------
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# ------------------- Parámetros -------------------
estados_dir = "estados"
energy_file = "energy.dat"
L           = 5.0               # lado de la caja
interval_ms = 0.0000001                  # 1 ms → ~1000 FPS (máximo práctico)

print("Cargando datos…")
energy = np.loadtxt(energy_file, comments='#')
steps_E, E_tot = energy[:,0].astype(int), energy[:,1]
files = sorted(glob.glob(os.path.join(estados_dir, "estado_*.dat")))
n_steps = len(files)

sample = np.loadtxt(files[0], comments='#')
N_part = sample.shape[0]                     
# columnas esperadas: idx, rx, ry, rz, fx, fy, fz, V
cols = sample.shape[1]

pos   = np.zeros((n_steps, N_part, 3))   # (paso, part, xyz)
force = np.zeros((n_steps, N_part, 3))

for i, f in enumerate(files):
    dat = np.loadtxt(f, comments='#')
    pos[i]   = dat[:,1:4]                # rx,ry,rz
    force[i] = dat[:,4:7]                # fx,fy,fz
print(f"Listo → {n_steps} pasos, {N_part} partículas")

# ------------------- Figura -------------------
fig = plt.figure(figsize=(14,7))
gs  = fig.add_gridspec(1,2, width_ratios=[2,1.5], wspace=0.3)

# ---- Energía ----
axE = fig.add_subplot(gs[0])
axE.plot(steps_E, E_tot, 'b-', lw=2)
axE.set_xlabel('Paso'); axE.set_ylabel('E total')
axE.set_title('Energía potencial')
axE.grid(alpha=0.3)
vline = axE.axvline(0, c='red', lw=2)

# ---- 3D ----
ax3 = fig.add_subplot(gs[1], projection='3d')
scatter = ax3.scatter([], [], [], c='b')
quiver  = ax3.quiver([], [], [], [], [], [], length=0.3,
                     normalize=True, color='red', alpha=0.7)

ax3.set(xlim=(0,L), ylim=(0,L), zlim=(0,L),
        xlabel='X', ylabel='Y', zlabel='Z')
ax3.set_title('Paso 0')

# caja (una sola vez)
corners = np.array([[0,0,0],[L,0,0],[L,L,0],[0,L,0],
                    [0,0,L],[L,0,L],[L,L,L],[0,L,L]])
edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
         (0,4),(1,5),(2,6),(3,7)]
for a,b in edges:
    ax3.plot(*zip(corners[a], corners[b]), c='k', lw=1, alpha=0.5)

# ------------------- Animación -------------------
def animate(frame):
    # energía
    vline.set_xdata([frame, frame])

    # 3D
    rx, ry, rz = pos[frame].T
    fx, fy, fz = force[frame].T
    scatter._offsets3d = (rx, ry, rz)

    # quiver es el único objeto “pesado”; lo borramos y recreamos
    global quiver
    quiver.remove()
    quiver = ax3.quiver(rx, ry, rz, fx, fy, fz,
                        length=0.3, normalize=True,
                        color='red', alpha=0.7)

    ax3.set_title(f'Paso {frame}')
    return scatter, vline

anim = FuncAnimation(fig, animate, frames=n_steps,
                     interval=interval_ms, blit=False, repeat=True)


plt.show()