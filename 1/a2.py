"Løsninger til oppgavene i a2, skrevet som et sammarbeid av både Bendik og Ole Kristian"

# a21

import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as cn

def gaussian(sigma, mu, vx):
	return 1/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*((vx - mu)/sigma)**2)

T = 3000     # gass temperature
N = 10**5    # number of H2 molecules
m = cn.m_H2  # mass of a hydrogen molecule

sigma = np.sqrt((T*cn.k_B)/m)
vx = np.linspace(-2.5*10**4, 2.5*10**4, 1000)

plt.plot(vx, gaussian(sigma, 0, vx))
plt.xlabel("vx [m/s]")
plt.ylabel("probability []")
plt.savefig("a21.png")
plt.clf
""()

# a22

vx = np.linspace(5*10**3, 30*10**3, 1000)
prob = np.trapz(gaussian(sigma, 0, vx), vx)
print(f"The probability of a H2-molecule having a velocity in the range is {prob*100:.3g} %")
print(f"Number of H2-molecules in the velocity range is {prob*N:.0f}")


# a23

def maxboltz(m, T, v):
	return (np.sqrt(m/(2*np.pi*cn.k_B*T))**3)*np.exp(-0.5*((m*v**2)/(cn.k_B*T)))*4*np.pi*v**2

v = np.linspace(0, 3*10**4, 1000)

plt.plot(v, maxboltz(m, T, v))
plt.xlabel("v [m/s]")
plt.ylabel("probability []")
plt.savefig("a23.png")

# Hvordor går ikke dette imot det vi viste i a21? tja, godt spørsmål egentlig.
