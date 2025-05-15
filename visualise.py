import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 0 # Dual Energy Flag

dnamein = '../../../../../ix/eschneider/hjl28/data/cloud_wind/4/8retry/hdf5/'
dnameout = '../../../../../ix/eschneider/hjl28/data/cloud_wind/4/8retry/png/'

CAT = 1

# t_cc = 4.89e4 # (vwind = 10 km/s)
# t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
istart = 0
iend = 1
time = 0

Tmin = 1e5
Tmax = 2e6

nmin = 14
nmax = 16

vmin = -200 #-200
vmax = 1200 #1200

for i in range(istart, iend):
    
    if CAT:
        f = h5py.File(dnamein + str(i) + '_slice.h5', 'r') # open the hdf5 file for reading
    else:
        f = h5py.File(dnamein + str(i) + '_slice.h5.0', 'r') # open the hdf5 file for reading
    head = f.attrs # read the header attributes into a structure, called head

    gamma = head['gamma'] # ratio of specific heats
    t  = head['t'] # time of this snapshot, in kyr
    nx = head['dims'][0] # number of cells in the x direction
    ny = head['dims'][1] # number of cells in the y direction
    nz = head['dims'][2] # number of cells in the z direction
    dx = head['dx'][0] # width of cell in x direction
    l_c = head['length_unit']
    t_c = head['time_unit']
    m_c = head['mass_unit']
    d_c = head['density_unit']                 
    v_c = head['velocity_unit']
    e_c = head['energy_unit']
    p_c = e_c # pressure units are the same as energy density units, density*velocity^2/length^3

    d  = f['d_xy'][:]
    px  = f['mx_xy'][:]
    py  = f['my_xy'][:]
    pz  = f['mz_xy'][:]
    E = f['E_xy'][:]

    if DE:
        GE = f['GE_xy'][:]

    f.close()

    n = d * d_c/ (mu*mp) # number density, particles per cm^3  

    vx = px/d
    vy = py/d
    vz = pz/d
        
    if not DE:    
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE

    T = GE*(gamma-1.0)*p_c / (n*kb) #temperature
    logT = np.log10(T)

    km = 1e-5

    Vx = vx*v_c*km #velocity in the x direction

    if CAT:
        f = h5py.File(dnamein + str(i) + '_proj.h5', 'r') # open the hdf5 file for reading
    else:
        f = h5py.File(dnamein + str(i) + '_proj.h5.0', 'r') # open the hdf5 file for reading
    head = f.attrs # read the header attributes into a structure, called head
    d  = f['d_xy'][:]
    # T = f['T_xy'][:]    


    d = d * m_c / (l_c**2)
    n = d / (mu*mp) # number density, particles per cm^3  
    logn = np.log10(n)

   
    P = np.log10(n*kb*T)

    f.close()

    # print(d_c)
    # print(p_c)
    # print('\t min \t\t\t max')
    # print('n: ', np.min(logn) , '\t' , np.max(logn))
    # print('T: ', np.min(logT) , '\t' , np.max(logT))
    # print('Vx: ', np.min(Vx) , '\t' , np.max(Vx))

    subplots = [logT.T, P.T, Vx.T]
    mins = [Tmin, nmin, vmin]
    maxs = [Tmax, nmax, vmax]
    cmaps = ['plasma', 'viridis', 'YlOrRd']
    labels = ['$log_{10}(K)$', '$log_{10}(N_{H})$ [$cm^{-2}$]', '$kms^{-1}$']

    fig, axs = plt.subplots(nrows=len(subplots), ncols=1)
    fig_color = 'white'
    bg_color = 'black'

    for j in range(len(subplots)):

        im = axs[j].imshow(subplots[j], cmap=cmaps[j]) #, vmin=mins[j], vmax = maxs[j]
        # axs[j].set_ylabel(labels[j], size=10, color=fig_color)
        axs[j].set_xticks(np.linspace(0,nx,9))
        axs[j].set_yticks(np.linspace(0,nz,5))
        axs[j].invert_yaxis()

        plt.setp(axs[j].spines.values(), color=fig_color)
        plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color=fig_color)

        if j == (len(subplots)-1):
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=6)
            axs[j].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.4),1))
            # print(nx*dx)
            [l.set_visible(False) for (i,l) in enumerate(axs[j].xaxis.get_ticklabels()) if i % 2 != 0]
            axs[j].set_xlabel('$kpc$', size=8, color=fig_color)
        else:
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0)
            
        divider = make_axes_locatable(axs[j])
        cax = divider.append_axes('right', size = 0.10, pad = 0.17)
        cb = plt.colorbar(im, cax=cax)
        cb.set_ticks(np.round(np.linspace(mins[j], maxs[j], 5), 2))
        cax.tick_params(axis='y', direction='out', color = fig_color, labelcolor=fig_color, labelsize=6)
        cax.set_ylabel(labels[j], size=8, color=fig_color)
        cb.outline.set_edgecolor(fig_color)

    fig.text(0.5, 0.9, str(int(t/t_cc))+r' $t_{cc}$', size=8, color=fig_color)
    # fig.text(0.5, 0.9, str(t)+r' $Myr$', size=8, color=fig_color)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)