import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# --- Generování dat ---
# Theta: 0 až pi (polární), Phi: 0 až 2*pi (azimutální)
n_theta, n_phi = 101, 101
theta = np.linspace(0, np.pi, n_theta)
phi = np.linspace(0, 2 * np.pi, n_phi)
THETA, PHI = np.meshgrid(theta, phi)

# Převod na kartézské souřadnice
X = np.sin(THETA) * np.cos(PHI)
Y = np.sin(THETA) * np.sin(PHI)
Z = np.cos(THETA)

# Tvoje data – zde jako příklad (sférické harmonické Y_2^0)
values = np.loadtxt(r'C:\results\Rabi\schnellbruder\measurement\negativity_scan_Rabi(2j=1, ω=1, R=20, λ=1.5, δ=0.5)_-0.283_31.5_0.txt')
VALUES = np.transpose(values.reshape(n_phi, n_theta))  # tvar (101, 101)

# --- Matplotlib verze ---
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

ax.set_box_aspect([1, 1, 1])
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)

# Normalizace hodnot na colormapu
norm = plt.Normalize(VALUES.min(), VALUES.max())
colors = cm.turbo(norm(VALUES))

surf = ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1,
                antialiased=False, shade=True, alpha=0.9)

# Colorbar (přes ScalarMappable)
mappable = cm.ScalarMappable(cmap='turbo', norm=norm)
mappable.set_array([])
plt.colorbar(mappable, ax=ax, shrink=0.4, label='Hodnota')

# Bod na povrchu koule
theta_arrow = np.pi / 4   # 45°
phi_arrow = np.pi / 3     # 60°

# Kartézské souřadnice špičky šipky
x = np.sin(theta_arrow) * np.cos(phi_arrow)
y = np.sin(theta_arrow) * np.sin(phi_arrow)
z = np.cos(theta_arrow)

ax.scatter([x], [y], [z], color='black', s=30, zorder=5)

ax.plot([0, 0], [0, 0], [1, 1.3], color='black', linewidth=2, linestyle='dashed')

ax.set_title('Heatmapa na 3D kouli')
ax.set_axis_off()
plt.tight_layout()
plt.show()