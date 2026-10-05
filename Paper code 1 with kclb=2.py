import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# ======================================
# Parameters from paper
# ======================================

omega = 1.0
Delta_c = -0.8*omega

Nc = 10      # photon cutoff
Nm = 10      # Landau cutoff

# Fig.1(a)
xi = 2.0

# For Fig.1(b) replace by:
# xi = np.sqrt(5)

x0_values = np.pi*np.array([
    -1,
    -3/4,
    -1/2,
    -1/4,
     0,
     1/4,
     1/2,
     3/4
])

eta_values = np.linspace(0,2,80)

# ======================================
# Operators
# ======================================

a = tensor(destroy(Nc), qeye(Nm))
adag = a.dag()

b = tensor(qeye(Nc), destroy(Nm))
bdag = b.dag()

# Free Hamiltonian
H0 = (
    -Delta_c*adag*a
    + omega*(bdag*b + 0.5)
)

# ======================================
# Figure 1(a,b)
# ======================================

plt.figure(figsize=(7,5))

for eta in eta_values:

    for x0 in x0_values:

        theta = (
            xi/np.sqrt(2)*(b+bdag)
            + x0
        )

        H = H0 + eta*(adag+a)*theta.cosm()

        E = H.eigenenergies()

        for n in range(32):
            plt.plot(eta, E[n],'.b',linewidth=0.5)

plt.xlabel(r'$\eta/\omega$')
plt.ylabel(r'$E/\omega$')
plt.title(r'$k_c l_B = 2$')
plt.show()