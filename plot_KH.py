import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable 

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K

iend = 200
dnamein = '../../../data/KH/KH_d/d_3_1/'
dnameout = '../../../plots/KH/KH_d/d_3_1/'

#plot = input("Enter 'd' 'P' or 'T': ")

for i in range(iend):

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

    d  = f['density'][:]
    #GE = f['GasEnergy'][:]

    f.close()

    mu = 0.6 # mean molecular weight (mu) (we should add this to the header, when relevant)

    d_cgs = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable

    n = d_cgs/(mu*mp) # number density, particles per cm^3


    #P = GE * (gamma - 1.0) * p_c #pressure

    #T = GE * (gamma - 1.0) * p_c / (n * kb) 

    #print(np.min(np.log10(d_cgs)))
    #print(np.max(np.log10(d_cgs)))

    fig, ax = plt.subplots()
    ax.set_xticks(ny*np.arange(0.25, 1, 0.25))
    ax.set_yticks(nz*np.arange(0.25, 1, 0.25))
    ax.tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1)
    ax.set_title('d1 = 3.0  d2 = 1.0')
    ax.text(0.03, 0.95, '512 res', transform=ax.transAxes, color='white')
    image = ax.imshow(np.log10(d_cgs).T, origin='lower', cmap='gist_heat', vmin=-31.3, vmax=-30.6) #originally -31.4 and -30.8

    # add a colorbar
    divider = make_axes_locatable(ax)
    cbax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(image, cax = cbax)
    cbax.tick_params(axis='y', direction='in')
    cb.solids.set_edgecolor('face')
    cbax.set_ylabel(r'$\mathrm{log}_{10}(\rho_A)$ [$\mathrm{g}\mathrm{cm}^{-2}$]')
    #plt.show()

    # save the figure
    plt.savefig(dnameout + 'd_3_1_' +str(i)+ '.png', dpi=300, transparent=False)
    plt.close()
        