## plots the x-direction velocity vs time for three different density ratios
## can use different resolutions instead of density ratios

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) (we should add this to the header, when relevant)

iend = 200 #timescale
dnameout = '../../../plots/KH/'

###################################### 1st Simulation ################################

dnamein = '../../../data/KH/KH_res/KH_256/'

vx_256 = np.empty(iend)

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

    d = f['density'][:]
    px  = f['momentum_x'][:]

    f.close()

    km = 1e-5

    Vx = (px*v_c*km)/d #velocity in the x direction

    vx_mean = np.mean(np.log10(np.abs(Vx)))

    vx_256[i] = vx_mean

################################# 2nd Simulation ##########################

dnamein = '../../../data/KH/KH_res/KH_512/'

vx_512 = np.empty(iend)

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

    d = f['density'][:]
    px  = f['momentum_x'][:]

    f.close()

    km = 1e-5

    Vx = (px*v_c*km)/d #velocity in the x direction

    vx_mean = np.mean(np.log10(np.abs(Vx)))

    vx_512[i] = vx_mean

########################### 3rd simulation ############################

dnamein = '../../../data/KH/KH_res/KH_1024/'

vx_1024 = np.empty(iend)

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

    d = f['density'][:]
    px  = f['momentum_x'][:]

    f.close()

    km = 1e-5

    Vx = (px*v_c*km)/d #velocity in the x direction

    vx_mean = np.mean(np.log10(np.abs(Vx)))

    vx_1024[i] = vx_mean


#print("n range    = %e %e" % (np.min(vx_avgs),np.max(vx_avgs)))
#vx_min = 5.57
#vx_max = 5.69
#print(vx_avgs)

times = np.arange(iend)
with plt.style.context("ggplot"):
    plt.plot(times, vx_256, label="256")
    plt.plot(times, vx_512, label="512")
    plt.plot(times, vx_1024, label="1024")

    plt.legend()
    plt.title('x-Direction Velocity vs Time')
    plt.ylabel(r'$\mathrm{log}_{10}(v_x)$ [$\mathrm{km}\mathrm{s}^{-1}$]')
    plt.tight_layout()
    #plt.ylim(vx_min, vx_max)

    # save the figure
    plt.savefig(dnameout + 'kh_vx_res'+ '.png', dpi=300, transparent=False)
    plt.close()