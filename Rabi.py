import matplotlib.pyplot as plt
import numpy as np

from matplotlib.animation import FuncAnimation, FFMpegWriter

from qutip import *
import time
from scipy.sparse.linalg.eigen.arpack.arpack import _augmented_orthonormal_cols

plt.ion()
class Rabi:
    def __init__(self, N=100, omega=1, R=100, gamma=0, delta=0, mu=0, nu=0, adaptiveN=True):
        self.N = N
        self.omega = omega
        self.R = R
        self.gamma = gamma
        self.delta = delta
        self.mu = mu
        self.nu = nu

        self.adaptiveN = adaptiveN

        self.RecalculateConstants()

    def RecalculateConstants(self):
        self.gammac = 0.5 * self.omega
        self.gamma0 = self.gammac / np.abs(self.delta)

        if self.adaptiveN:
            self.N = int(12 * self.R * self.gamma**2) + 100
            self.numEV = int(4 * self.R * self.gamma**2) + 100
            print(f"N = {self.N}, numEV = {self.numEV}")
        else:
            self.numEV = self.N + 1

    def Hamiltonian(self):
        self.RecalculateConstants()

        # Pre-compute operators for the hamiltonian
        a  = tensor(destroy(self.N), qeye(2))
        Jm = tensor(qeye(self.N), destroy(2))
        n = a.dag() * a
        Jz = Jm.dag() * Jm

        H0 = self.omega * n + self.R * self.omega * (Jz - 0.5)
        Hgamma =  (1 + self.delta) * (a.dag() * Jm + a * Jm.dag()) + (1 - self.delta) * (a.dag() * Jm.dag() + a * Jm)
        Hmu = 2 * (a.dag() + a) * (Jz - 0.5)
        Hnu = a.dag() + a

        return H0 + np.sqrt(self.R) * (self.gamma * Hgamma + self.mu * Hmu + self.nu * Hnu)

    def LevelDynamics(self, gammas=None, deltas=None, numEV=100, limits=(-2, 2)):
        energies = []
        xs = None

        fig, ax = plt.subplots(figsize=(8,6))

        def CalculateEnergies():
            startTime = time.time()
            evals, evecs = self.Hamiltonian().eigenstates(eigvals=numEV)
            energies.append(np.array(2 * evals / self.R))
            print(f"N={self.N}, lambda={self.gamma}, delta={self.delta}, R={self.R}, time={time.time() - startTime}")

        if gammas is not None:
            xs = gammas
            ax.set_xlabel(r"$\lambda}$")
            ax.set_title(f'Level Dynamics for R={self.R}, delta={self.delta}')
    
            for gamma in gammas:
                self.gamma = gamma
                CalculateEnergies()

        elif deltas is not None:
            xs = deltas
            ax.set_xlabel(r"$\delta$")
            ax.set_title(f'Level Dynamics for R={self.R}, lambda={self.gamma}')

            for delta in deltas:
                self.delta = delta
                evals, evecs = self.Hamiltonian().eigenstates(eigvals=numEV)
                CalculateEnergies()

        ax.plot(xs, energies, color="black", alpha=0.3)
        ax.set_ylabel('E')
        ax.set_ylim(limits)

        plt.show()
    
    def TimeEvolution(self, times, psi=None, operator=None):
        if psi is None:
            psi = ket("00", [self.N, 2])
            
        # https://qutip.org/docs/latest/guide/dynamics/dynamics-options.html
        options = Options()
        options.nsteps = 1000000
        options.num_cpus = 8

        return sesolve(self.Hamiltonian(), psi, times, options=options, e_ops=[operator])

    def Nt(self, times, psi=None):
        a  = tensor(destroy(self.N), qeye(2))
        n = a.dag() * a

        return self.TimeEvolution(times, psi, operator=n)

    def St(self, times, psi=None):
        if psi is None:
            psi = ket("00", [self.N, 2])

rabi = Rabi(R=100, delta=0.5, gamma=1, adaptiveN=True)
times = np.linspace(0, 100, 101)
solution = rabi.Nt(times)

plt.plot(solution.times, solution.expect[0])
plt.show()


#
# Main calculation
#
def LD():
    """ Figure 5 of Machine.pdf """
    rabi = Rabi(R=5, delta=0.5, adaptiveN=True)
    gammas = np.linspace(0, 2.5, 251)

    rabi.LevelDynamics(gammas=gammas, numEV=100, limits=(-4,4))

def Figure5():
    """ Figure 5 of Machine.pdf """
    rabi = Rabi(N=500, R=np.sqrt(20), delta=0.5, adaptiveN=False)
    gammas = np.linspace(0, 2.5, 501)

    rabi.LevelDynamics(gammas=gammas, numEV=300)

def Figure6():
    """ Figure 6 of Machine.pdf """
    rabi = Rabi(N=500, R=np.sqrt(20), gamma=2, adaptiveN=False)
    deltas = np.linspace(-1, 1, 201)
    rabi.LevelDynamics(deltas=deltas, numEV=300)

def Figure2E1():
    "Figure 2 of Extended1.pdf"
    rabi = Rabi(N=500, R=20, delta=0.5, mu=0.4, adaptiveN=False)
    gammas = np.linspace(0, 2, 501)

    rabi.LevelDynamics(gammas=gammas, numEV=300, limits=(-2, 0))

def Figure2E2():
    "Figure 2 of Extended1.pdf"
    rabi = Rabi(N=500, R=20, delta=0.5, mu=0.4, nu=0.4, adaptiveN=False)
    gammas = np.linspace(0, 2, 501)

    rabi.LevelDynamics(gammas=gammas, numEV=300, limits=(-2, 0))

def Wigner():
    rabi = Rabi(R=100, delta=0.5, gamma=1, adaptiveN=True)
    times = np.linspace(0, 50, 51)

    w = rabi.Wigner(times)

    x = np.linspace(-25, 25, 251)
    X,Y = np.meshgrid(x, x)

    fig = plt.figure(2, figsize=(9, 6))

    for state, t in zip(w.states, times):
        startTime = time.time()
        rho = ptrace(ket2dm(state), 0)
        W = wigner(rho, x, x)
        plt.clf()
        plt.contourf(X, Y, W, np.linspace(-0.1, 0.1, 101), cmap=plt.cm.seismic, extend = "both")        
        print(t, time.time() - startTime)
        #plt.colorbar(shrink=0.65, aspect=20)
        plt.colorbar()
        plt.title(f"T={t}")
        plt.show()
        plt.savefig(f"d:\\results\\Wigner\\{t}.png")
        plt.pause(1)

def WignerAnimation(gamma, times):
    rabi = Rabi(R=100, delta=0.5, gamma=gamma, adaptiveN=True)
    w = rabi.Wigner(times)

    limits= (-25 * gamma, 25 * gamma)
    x = np.linspace(*limits, 100)
    X,Y = np.meshgrid(x, x)

    fig, ax = plt.subplots()

    def init():
        ax.clear()
        cf = ax.contourf(X, Y, np.zeros(shape(X)), np.linspace(-0.4, 0.4, 201), cmap=plt.cm.seismic, extend = "both")
        fig.colorbar(cf, ax=ax)

        return ax,

    def animate(i):
        startTime = time.time()
        state = w.states[i]
        rho = ptrace(ket2dm(state), 0)
        W = wigner(rho, x, x)
        ax.clear()
        cf = ax.contourf(X, Y, W, np.linspace(-0.4, 0.4, 201), cmap=plt.cm.seismic, extend = "both")
        ax.set_title(f"t = {times[i]:.1f}")

        print(f"{time.time() - startTime}: {i}")
        return ax,

    interval = 40            # In milliseconds

    anim = FuncAnimation(fig, animate, init_func=init, interval=interval, blit=False, frames=len(w.states))
    anim.save(f"wigner_{gamma}.mp4", writer='ffmpeg')

# If we want to save the animation, we have to switch off blit!!!
#anim = FuncAnimation(fscrig, animate, interval=interval, frames=len(randomWalk3D))
#anim.save("randomwalk.gif", writer='imagemagick')
#anim.save("randomwalk.mp4", extra_args=['-vcodec', 'libx264'])
#plt.show()

"""
times = np.linspace(0, 200, 1001)

for gamma in np.linspace(0.2, 2, 10):
    WignerAnimation(gamma, times)
"""
#Wigner()
#LD()

