import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) (we should add this to the header, when relevant)

iend = 200
dnamein = '../../../data/KH/KH_res/KH_256/'
dnameout = '../../../plots/KH/'

################################################## 256 res ######################################3

d_256 = np.empty(iend)

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

    d_cgs = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable
    dx_cgs = dx*l_c
    n = d_cgs/(mu*mp) # number density, particles per cm^3

    d_mean = np.mean(np.log10(d_cgs))

    d_256[i] = d_mean

##################################### 512 res #######################################

dnamein = '../../../data/KH/KH_res/KH_512/'

d_512 = np.empty(iend)

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

    d_cgs = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable
    dx_cgs = dx*l_c
    n = d_cgs/(mu*mp) # number density, particles per cm^3

    d_mean = np.mean(np.log10(d_cgs))

    d_512[i] = d_mean

################################ 1024 res #############################

dnamein = '../../../data/KH/KH_res/KH_1024/'

d_1024 = np.empty(iend)

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

    d_cgs = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable
    dx_cgs = dx*l_c
    n = d_cgs/(mu*mp) # number density, particles per cm^3

    d_mean = np.mean(np.log10(d_cgs))

    d_1024[i] = d_mean

times = np.arange(iend)

with plt.style.context("ggplot"):
    fig = plt.plot(times, d_256, label="256")
    plt.plot(times, d_512, label="512")
    plt.plot(times, d_1024, label="1024")
    plt.legend()
    plt.title('Density vs Time')
    plt.ylabel(r'$\mathrm{log}_{10}(\rho_A)$ [$\mathrm{g}\mathrm{cm}^{-2}$]')
    plt.tight_layout()
    # save the figure
    plt.savefig(dnameout + 'kh_d'+ '.png', dpi=300, transparent=False)
    plt.close()