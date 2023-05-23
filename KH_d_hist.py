import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) 

iend = 200
dnameout = '../../../plots/KH/d_hist/single/'

######################## v_05 ################################

dnamein = '../../../data/KH/KH_d/d_1_2/'

for i in range(iend):

    f = h5py.File(dnamein+str(i)+'.h5.0', 'r')

    head = f.attrs

    gamma = head['gamma'] #ratio of specific heats
    t = head['t'] #time of snapshot (kyr)
    nx = head['dims'][0] # number of cels in he x direction
    ny = head['dims'][1] # number of cells in the y direction
    nz = head['dims'][2] # number of cells in the z direction
    dx = head['dx'][0] #width of cells in the x direction
    dy = head['dx'][1] #width of cells in the y direction
    dz = head['dx'][2] #width of cells in the z direction
    l_c = head['length_unit']
    t_c = head['time_unit']
    m_c = head['mass_unit']
    d_c = head['density_unit']
    v_c = head['velocity_unit']
    e_c = head['energy_unit']
    p_c = e_c

    d = f['density'][:]

    f.close()

    d_cgs = d*d_c

    d_cgs = d_cgs.flatten()

    ####################### PLOTTING ################################

    with plt.style.context("ggplot"):
        
        fig, ax = plt.subplots(figsize=(6,4))

        ax.hist(np.log10(d_cgs), bins=50, alpha=0.8)
        ax.set_title("Density of Cells")

        ax.set_yticks([])
        ax.tick_params(left=False, bottom=False)
        ax.set_xlim(-31.25,-30.8)
        ax.set_ylim(0, 0.5*nx*nx)
        
        ax.set_ylabel('Frequency')
        ax.set_xlabel(r'$\mathrm{log}_{10}(\rho_A)$ [$\mathrm{g}\mathrm{cm}^{-2}$]')

        fig.tight_layout()
        # save the figure
        plt.savefig(dnameout + 'd_hist_' + str(i) + '.png', dpi=300, transparent=False)
        plt.close()