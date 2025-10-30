import matplotlib.pyplot as plt
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D

# =============================
# Leer L desde argumento
# =============================
if len(sys.argv) < 2:
    print("Uso: python plot.py L")
    sys.exit(1)

L = float(sys.argv[1])
sigma = 1.0
d_c = 2.5 * sigma
r0 = 2**(1.0/6.0) * sigma  # ≈ 1.122
# =============================
# Leer datos
# =============================
data = np.loadtxt('forces.dat', skiprows=1)  # Salta header
i, rx, ry, rz, fx, fy, fz, V, r_contx,r_conty,r_contz = data.T

# =============================
# Configurar figura 3D
# =============================
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Partículas y fuerzas
ax.scatter(rx, ry, rz, c='b', marker='o', label='Posiciones')

ax.quiver(rx, ry, rz, fx, fy, fz, color='r', length=0.1, normalize=True, label='Fuerzas')

for idx in range(len(rx)):
    x, y, z = rx[idx], ry[idx], rz[idx]

    # NÚMERO CERCA DE LA PARTÍCULA
    ax.text(x, y, z + 0.15, str(idx), color='black', fontsize=12, 
            fontweight='bold', ha='center', va='bottom', zorder=10)

# =============================
# Dibujar bordes de la caja
# =============================
# Esquinas de la caja
corners = np.array([[0,0,0],[L,0,0],[L,L,0],[0,L,0],
                    [0,0,L],[L,0,L],[L,L,L],[0,L,L]])

# Conectar las esquinas con líneas
edges = [
    (0,1),(1,2),(2,3),(3,0), # base
    (4,5),(5,6),(6,7),(7,4), # techo
    (0,4),(1,5),(2,6),(3,7)  # verticales
]

for e in edges:
    x = [corners[e[0],0], corners[e[1],0]]
    y = [corners[e[0],1], corners[e[1],1]]
    z = [corners[e[0],2], corners[e[1],2]]
    ax.plot(x, y, z, c='k', lw=1)
    
# =============================

x0, y0, z0 = rx[0], ry[0], rz[0]
for j in range(1, len(rx)):

    # --- PBC: imagen mínima ---
    dx = rx[j] - x0
    dy = ry[j] - y0
    dz = rz[j] - z0
    dx = dx - L * round(dx / L)
    dy = dy - L * round(dy / L)
    dz = dz - L * round(dz / L)
    d = np.sqrt(dx**2 + dy**2 + dz**2)

    # --- Solo dibuja si d < d_c ---
    if d < 2.5:
        # Posición imagen
        xj_img = x0 + dx
        yj_img = y0 + dy
        zj_img = z0 + dz

        # Color según distancia
        if d < r0:
            color = 'red'
            label = 'Repulsiva (d < 1.12σ)' if j == 1 else ""
        else:
            color = 'blue'
            label = 'Atractiva (1.12σ ≤ d < 2.5σ)' if j == 1 else ""

        ax.plot([x0, xj_img], [y0, yj_img], [z0, zj_img],
                c=color, lw=2.5, alpha=0.8, label=label)
    
x0, y0, z0 = rx[0], ry[0], rz[0]
x_cont0, y_cont0, z_cont0 = r_contx[0], r_conty[0], r_contz[0]

for j in range(1, len(rx)):
    # vector diferencia
    dx = rx[j] - x0
    dy = ry[j] - y0
    dz = rz[j] - z0
    rdist = np.sqrt(dx**2 + dy**2 + dz**2)
    #dx = rx[j] - x_cont0
    #dy = ry[j] - y_cont0
    #dz = rz[j] - z_cont0
    #rdist_cont = np.sqrt(dx**2 + dy**2 + dz**2)
    #print(str(rdist)+ "   "+ str(rdist_cont))
    if rdist > 0.5*L:
        print(dx)
        print(dy)
        print(dz)
        # dibujar línea (ajustando para imagen mínima)
        ax.plot([x0, x0+dx], [y0, y0+dy], [z0, z0+dz], c='m', lw=1.5, label='> L/2' if j==1 else "")
#        ax.plot([x0, x0+dx], [y0, y0+dy], [z0, z0+dz], c='y', lw=1.5, label='> L/2' if j==1 else "")

# =============================
# Etiquetas y título
# =============================
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Fuerzas en partículas (Lennard-Jones)')

ax.legend()
plt.show()
