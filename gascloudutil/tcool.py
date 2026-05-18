import matplotlib.pyplot as plt
import numpy as np
import unyt
from cholla_cooling import ChollaEOS
from estimate_survival_radius import find_minmix

cholla_cie_eos = ChollaEOS()

# you can specify the cloud temperature
tcool = find_minmix(
    Tcl = unyt.unyt_quantity(1e4, 'K'),
    Tw = unyt.unyt_quantity(1e6, 'K'),
    assumed_p = 1e4 * unyt.K/unyt.cm**3 * unyt.kboltz_cgs,
    eos = cholla_cie_eos)

print('\n Cooling Time:', tcool)

#  Cooling Time: {'rho_minmix': unyt_quantity(5.76981054e-26, 'g/cm**3'), 'e_minmix': unyt_quantity(3.58932617e+13, 'cm**2/s**2'), 'T_minmix': unyt_quantity(173934.04947458, 'K'), 'tcool_minmix': unyt_array(-1.11447332e+12, 's')}

