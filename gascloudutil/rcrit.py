import matplotlib.pyplot as plt
import numpy as np
import unyt
from cholla_cooling import ChollaEOS
from estimate_survival_radius import estimate_survival_radius

cholla_cie_eos = ChollaEOS()

# you can specify the cloud temperature
r_crit_a = estimate_survival_radius(
    eos = cholla_cie_eos, vwind = 1000*unyt.km/unyt.s,
    density_contrast = 100, alpha = 7,
    pressure = 1e4 * unyt.K/unyt.cm**3 * unyt.kboltz_cgs,
    temperature_cl = unyt.unyt_quantity(1e4, 'K'))

print('\n Critical Survival Radius:', r_crit_a)


# r_crit_b = estimate_survival_radius(
#     eos = cholla_cie_eos, vwind = 1000*unyt.km/unyt.s,
#     density_contrast = 1000, alpha = 7,
#     pressure = 1e3 * unyt.K/unyt.cm**3 * unyt.kboltz_cgs,
#     temperature_w = unyt.unyt_quantity(3e7, 'K'))

# print('Critical Survival Radius:', r_crit_b)