# Kód od Radima

import matplotlib.pyplot as plt
import numpy as np

from qutip import *

# Frequency
omega = 2 * np.pi
R = 100

# Coupling constant
g = 0.7 * 2 * np.pi

kappa = 0
gamma = 0
N = 50

temp = 0

p = 0.75
tau = 0.01
gd = 0.001

tlist = np.linspace(0, 2, 2000)

a = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
parity = 1j * np.pi * (a.dag() * a)
parity = parity.expm()

H0 = omega * a.dag() * a + R * omega * sm.dag() * sm + gd * sm.dag() * sm * (a + a.dag()) + gdr * (a + a.dag())
Hi = g * np.sqrt(R) * ((1 + p) * (a.dag() * sm + a * sm.dag()) + (1 - p) * (a.dag() * sm.dag() + a * sm))

def Hi_c(tlist, args):
    return 1 - np.exp(-tlist / tau)

H = [H0, [Hi, Hi_c]]

c_op = []
rate = kappa * (1 + temp)
if rate > 0:
    c_op.append(np.sqrt(rate) * a)

rate = kappa * temp
if rate > 0:
    c_op.append(np.sqrt(rate) * a.dag())

rate = gamma
if rate > 0:
    c_op.append(np.sqrt(rate) * sm)

output = mesolve(H, psi0, tlist, c_op, [a.dag() * a, sm.dag() * sm])