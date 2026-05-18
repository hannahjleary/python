import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.integrate import quad
import seaborn as sns
import h5py
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

plt.style.use('classic')

# Constants
pi = np.pi

# Safe exponential to avoid overflow
def safe_exp(x):
    x = np.clip(x, None, 700)  # Prevent overflow
    return np.exp(x)

# Bose energy density
def rho_B(T, m, g, P0, P1):
    def integrand(p):
        E = np.sqrt(p**2 + m**2)
        return (p**2 * E) / (safe_exp(E / T) - 1)
    result, _ = quad(integrand, P0, P1)
    return (g / (2 * pi**2)) * result

# Bose pressure
def P_B(T, m, g, P0, P1):
    def integrand(p):
        E = np.sqrt(p**2 + m**2)
        return (p**4 / (3 * E)) / (safe_exp(E / T) - 1)
    result, _ = quad(integrand, P0, P1)
    return (g / (2 * pi**2)) * result

# Fermi energy density
def rho_F(T, m, g, P0, P1):
    def integrand(p):
        E = np.sqrt(p**2 + m**2)
        return (p**2 * E) / (safe_exp(E / T) + 1)
    result, _ = quad(integrand, P0, P1)
    return (g / (2 * pi**2)) * result

# Fermi pressure
def P_F(T, m, g, P0, P1):
    def integrand(p):
        E = np.sqrt(p**2 + m**2)
        return (p**4 / (3 * E)) / (safe_exp(E / T) + 1)
    result, _ = quad(integrand, P0, P1)
    return (g / (2 * pi**2)) * result

# Composite g(T)
def g_T(T):
    term1 = P_F(T, 0.5, 4.0, 0, 200)
    print(P_F(T, 0.5, 4.0, 0, 200))
    term2 = P_B(T, 0, 2, 0, 200)
    print(P_B(T, 0, 2, 0, 200))
    term3 = rho_F(T, 0.5, 4.0, 0, 200)
    print(rho_B(T, 0, 2, 0, 200))
    term4 = rho_B(T, 0, 2, 0, 200)
    print(rho_B(T, 0, 2, 0, 200))
    return (45 / (T**4 * 2 * pi**2)) * (term1 + term2 + term3 + term4)

# Initial temperature in energy units (e.g., eV)
T_CMB = 2.725  # Kelvin
k_B = 8.617e-11  # MeV/K
T_G = T_CMB * k_B  # Convert to MeV

# Evaluate g(T_G)
g_value = g_T(T_G)
print(f"g(T_G) = {g_value:.4e}")

# Scale factor range
a_vals = np.linspace(10**-10, 10**-8, 100)
print(a_vals)

# Temperature evolution
T_gamma = T_G / a_vals
T_nu = (11 / 2)**(1/3) * g_value**(-1/3) * T_gamma

# Plotting
plt.figure(figsize=(6, 4))
plt.plot(a_vals, T_gamma, label='Photon Temperature $T_\\gamma$', color='orange')
plt.plot(a_vals, T_nu, label='Neutrino Temperature $T_\\nu$', color='blue', linestyle='--')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Scale Factor $a$')
plt.ylabel('Temperature (meV)')
plt.title('Photon and Neutrino Temperature vs Scale Factor')
plt.legend()
# plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('decoupling.png', dpi=300, bbox_inches='tight', pad_inches=0.2)