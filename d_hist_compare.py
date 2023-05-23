import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) 

iend = 20
dnameout = '../../../plots/KH/KH_d/d_hist/comparison/'

######################## v 1 ################################

for i in range(iend):

    dnamein = '../../../data/KH/KH_v/v_1/'

    f = h5py.File(dnamein+str(i)+'.h5.0', 'r')

    head = f.attrs

    #print(f.keys())

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

    d1_cgs = d*d_c
    d1_cgs = d1_cgs.flatten()

    ######################## v 05 ################################

    dnamein = '../../../data/KH/KH_v/v_05/' 

    f = h5py.File(dnamein+str(i)+'.h5.0', 'r')

    head = f.attrs

    #print(f.keys())

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

    d2_cgs = d*d_c
    d2_cgs = d2_cgs.flatten()


######################## v 075 ################################

    dnamein = '../../../data/KH/KH_v/v_075/'

    f = h5py.File(dnamein+str(i)+'.h5.0', 'r')

    head = f.attrs

    #print(f.keys())

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

    d3_cgs = d*d_c
    d3_cgs = d3_cgs.flatten()



    ####################### PLOTTING ################################

    with plt.style.context("Solarize_Light2"):
        
        fig, ax = plt.subplots(3, sharex=True, figsize=(6,4))

        ax[0].hist(np.log10(d1_cgs), alpha=0.7)
        ax[1].hist(np.log10(d2_cgs), alpha=0.7)
        ax[2].hist(np.log10(d3_cgs), alpha=0.7)

        fig.suptitle("Density of Cells")

        ax.set_yticks([])
        ax.tick_params(left=False, bottom=False)
        ax.set_xlim(-31.4,-30.7)
        
        ax.set_ylabel('Frequency')
        ax.set_xlabel(r'$\mathrm{log}_{10}(\rho_A)$ [$\mathrm{g}\mathrm{cm}^{-2}$]')

        fig.tight_layout()
        # save the figure
        plt.savefig(dnameout + 'd_hist_comp_' + str(i) + '.png', dpi=300, transparent=False)
        plt.close()