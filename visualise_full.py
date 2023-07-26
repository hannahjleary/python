import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6

DE = 1 # Dual Energy Flag

dnamein='../../data/cloud_wind/4_largediff/16/hdf5/' # directory where the file is located
dnameout='../../data/cloud_wind/4_largediff/16/png/' # directory where the plot will be saved

iend = 500
t_cc = 4.89e2

for i in range(0, iend, 10):

    f = h5py.File(dnamein + str(i) + '.h5', 'r') # open the hdf5 file for reading
    head = f.attrs # read the header attributes into a structure, called head

    head.keys()

    f.keys()

    gamma = head['gamma'] # ratio of specific heats
    t  = head['t'] # time of this snapshot, in kyr
    nx = head['dims'][0] # number of cells in the x direction
    ny = head['dims'][1] # number of cells in the y direction
    nz = head['dims'][2] # number of cells in the z direction
    dx = head['dx'][0] # width of cell in x direction
    dy = head['dx'][1] # width of cell in y direction
    dz = head['dx'][2] # width of cell in z direction
    l_c = head['length_unit']
    t_c = head['time_unit']
    m_c = head['mass_unit']
    d_c = head['density_unit']                 
    v_c = head['velocity_unit']
    e_c = head['energy_unit']
    p_c = e_c # pressure units are the same as energy density units, density*velocity^2/length^3

    d  = f['density'][:]
    px  = f['momentum_x'][:]
    py  = f['momentum_y'][:]
    pz  = f['momentum_z'][:]
    E = f['Energy'][:]

    if DE:
        GE = f['GasEnergy'][:]

    f.close()

    n = d * d_c/(mu*mp) # number density, particles per cm^3  

    vx = px/d
    vy = py/d
    vz = pz/d
        
    if not DE:    
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE

    T = GE*(gamma-1.0)*p_c / (n*kb) #temperature

    km = 1e-5

    Vx = (px*v_c*km)/d #velocity in the x direction

    # print('\t min \t\t\t max')
    # print('n: ', np.min(n) , '\t' , np.max(n))
    # print('T: ', np.min(T) , '\t' , np.max(T))
    # print('Vx: ', np.min(Vx) , '\t' , np.max(Vx))

    #Temperature Projection
    # d_avg = np.average(d)
    # T_weighted = T*d/d_avg
    # T_y = np.sum(n, axis=1)
    # log_T_y = np.log10(T_y)

    #Temperature Slice
    T_slice_xz = T[:,int(ny/2),:]
    logT_slice_xz = np.log10(T_slice_xz)

    #Number Density Projection
    n_y = np.sum(n, axis=1)*dy*l_c
    log_n_y = np.log10(n_y)

    #Velocity in the x-direction Slice
    Vxslice_xz = Vx[:,int(ny/2),:]

    Tmin = 3.0
    Tmax = 7.0

    nmin = 18.0
    nmax = 30.0

    vmin = 0.0
    vmax = 1300.0


    subplots = [logT_slice_xz.T, log_n_y.T, Vxslice_xz.T]
    mins = [Tmin, nmin, vmin]
    maxs = [Tmax, nmax, vmax]
    cmaps = ['plasma', 'viridis', 'YlOrRd']
    labels = ['$log_{10}(K)$', '$log_{10}(N_{H})$ [$cm^{-2}$]', '$kms^{-1}$']

    fig, axs = plt.subplots(nrows=len(subplots), ncols=1)
    fig_color = 'white'
    bg_color = 'black'

    for j in range(len(subplots)):

        im = axs[j].imshow(subplots[j], cmap=cmaps[j], vmin=mins[j], vmax = maxs[j]) 
        # axs[j].set_ylabel(labels[j], size=10, color=fig_color)
        axs[j].set_xticks(np.linspace(0,nx,9))
        axs[j].set_yticks(np.linspace(0,nz,5))
        axs[j].invert_yaxis()

        plt.setp(axs[j].spines.values(), color=fig_color)
        plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color=fig_color)

        if j == (len(subplots)-1):
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=6)
            axs[j].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.2),1))
            [l.set_visible(False) for (i,l) in enumerate(axs[j].xaxis.get_ticklabels()) if i % 2 != 0]
            axs[j].set_xlabel('$kpc$', size=8, color=fig_color)
        else:
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0)
            
        divider = make_axes_locatable(axs[j])
        cax = divider.append_axes('right', size = 0.12, pad = 0.2)
        cb = plt.colorbar(im, cax=cax)
        cb.set_ticks(np.round(np.linspace(mins[j], maxs[j], 5), 2))
        cax.tick_params(axis='y', direction='out', color = fig_color, labelcolor=fig_color, labelsize=6)
        cax.set_ylabel(labels[j], size=8, color=fig_color)
        cb.outline.set_edgecolor(fig_color)


    fig.text(0.5, 0.9, str(int(t/t_cc))+r' $t_{cc}$', size=8, color=fig_color)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)