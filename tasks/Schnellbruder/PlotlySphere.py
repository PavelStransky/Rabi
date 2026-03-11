import numpy as np
import plotly.graph_objects as go

n_theta, n_phi = 101, 101
theta = np.linspace(0, np.pi, n_theta)
phi = np.linspace(0, 2 * np.pi, n_phi)
THETA, PHI = np.meshgrid(theta, phi)

X = np.sin(THETA) * np.cos(PHI)
Y = np.sin(THETA) * np.sin(PHI)
Z = np.cos(THETA)

# Tvoje hodnoty

values = np.loadtxt(r'C:\results\Rabi\schnellbruder\measurement\negativity_scan_Rabi(2j=1, ω=1, R=20, λ=1.5, δ=0.5)_-0.283_31.5_0.txt')
VALUES = np.transpose(values.reshape(n_phi, n_theta))  # tvar (101, 101)

fig = go.Figure(data=[
    go.Surface(
        x=X, y=Y, z=Z,
        surfacecolor=VALUES,   # <-- zde dosadíš svá data
        colorscale='turbo',
        colorbar=dict(title='Hodnota'),
        cmin=VALUES.min(),
        cmax=VALUES.max(),
    )
])

theta_arrow = np.pi / 4   # 45°
phi_arrow = np.pi / 3     # 60°

# Šipka pomocí Cone (hrot) + Scatter3d (dřík)
x_tip = np.sin(theta_arrow) * np.cos(phi_arrow)
y_tip = np.sin(theta_arrow) * np.sin(phi_arrow)
z_tip = np.cos(theta_arrow)

# Dřík šipky
fig.add_trace(go.Scatter3d(
    x=[0, x_tip], y=[0, y_tip], z=[0, z_tip],
    mode='lines',
    line=dict(color='black', width=5),
    showlegend=False
))

# Hrot šipky
fig.add_trace(go.Cone(
    x=[x_tip], y=[y_tip], z=[z_tip],
    u=[x_tip * 0.1], v=[y_tip * 0.1], w=[z_tip * 0.1],
    colorscale=[[0, 'black'], [1, 'black']],
    showscale=False,
    sizemode='absolute',
    sizeref=0.15,
    anchor='tip'
))

# Osa koule (od jižního pólu k severnímu)
fig.add_trace(go.Scatter3d(
    x=[0, 0], y=[0, 0], z=[-2, 1.5],  # trochu přesahuje povrch
    mode='lines+text',
    line=dict(color='black', width=5, dash='dash'),
    text=['', 'Z'],                        # popisek na severním pólu
    textposition='top center',
    textfont=dict(size=14, color='gray'),
    showlegend=False
))

# Rovnoběžky (konstantní theta)
for t in np.linspace(0, np.pi, 11)[1:-1]:   # vynech póly
    p = np.linspace(0, 2 * np.pi, 100)
    x = np.sin(t) * np.cos(p)
    y = np.sin(t) * np.sin(p)
    z = np.cos(t) * np.ones_like(p)
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines',
        line=dict(color='white', width=1),
        showlegend=False
    ))

# Poledníky (konstantní phi)
for p in np.linspace(0, 2 * np.pi, 21)[:-1]:  # vynech poslední (= první)
    t = np.linspace(0, np.pi, 100)
    x = np.sin(t) * np.cos(p)
    y = np.sin(t) * np.sin(p)
    z = np.cos(t)
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines',
        line=dict(color='white', width=1),
        showlegend=False
    ))

theta_bod = np.pi / 4
phi_bod = np.pi / 3

x_bod = np.sin(theta_bod) * np.cos(phi_bod)
y_bod = np.sin(theta_bod) * np.sin(phi_bod)
z_bod = np.cos(theta_bod)

# Malá kulička na povrchu
u = np.linspace(0, 2 * np.pi, 20)
v = np.linspace(0, np.pi, 20)
r = 0.05  # poloměr kuličky

x_s = x_bod + r * np.outer(np.cos(u), np.sin(v))
y_s = y_bod + r * np.outer(np.sin(u), np.sin(v))
z_s = z_bod + r * np.outer(np.ones_like(u), np.cos(v))

fig.add_trace(go.Surface(
    x=x_s, y=y_s, z=z_s,
    colorscale=[[0, 'black'], [1, 'black']],
    showscale=False,
    opacity=1.0
))

fig.update_layout(
    title='Interaktivní heatmapa na 3D kouli',
    scene=dict(
        camera=dict(
            eye=dict(x=0, y=2, z=2),      # kamera podél osy Y
            up=dict(x=0, y=0, z=1),        # osa Z směřuje nahoru
        ),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False),
    )
)

fig.show()