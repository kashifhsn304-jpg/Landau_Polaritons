import numpy as np
import matplotlib.pyplot as plt
from qutip import *

Na, Nb = 20, 20

a = tensor(destroy(Na), qeye(Nb))
b = tensor(qeye(Na), destroy(Nb))

Qb = (b + b.dag())/2
Pb = (b - b.dag())/(2j)

omega = 1.0
Delta_c = 0.8
xi = np.sqrt(0.4)

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

etas = np.linspace(0, 2, 50)

H0 = Delta_c*a.dag()*a + omega*b.dag()*b

plt.figure(figsize=(8,6))

for n, x0 in enumerate(x0_values):
    print(x0)
    Qb_var = np.zeros(len(etas))
    Pb_var = np.zeros(len(etas))

    theta = xi/np.sqrt(2)*(b + b.dag()) + x0
    C = theta.cosm()

    V = (a + a.dag())*C

    for i, eta in enumerate(etas):

        H = H0 + eta*V

        psi0 = H.groundstate()[1]

        Qb_var[i] = variance(Qb, psi0)
        Pb_var[i] = variance(Pb, psi0)

    if n == 0:

        plt.plot(etas, Qb_var, '-', label='Position')
        plt.plot(etas, Pb_var, '--', label='Momentum')

    else:

        plt.plot(etas, Qb_var, '-')
        plt.plot(etas, Pb_var, '--')

plt.xlabel(r'$\eta/\omega$')
plt.ylabel('Variance')
plt.legend()
plt.xlim(0, 2)
plt.ylim(0, 1.5)
plt.show()