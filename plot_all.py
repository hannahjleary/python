import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 1.0 # mean molecular weight (mu) of 1

DE = 0 # Dual Energy Flag

dnamein='../../data/cloud_wind/3/' # directory where the file is located
dnameout='../../data/cloud_wind/3/plots3/' # directory where the plot will be saved

sims = ['4/', '8/', '16/', '32/', '48/']
res_labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$', '$R_{32}$', '$R_{48}$']
cat = [False, False, True, True, True]

###### Cloud Crushing Times ###################
# t_cc = 4.89e4 # (vwind = 10 km/s)
t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)

istart = 0
iend = 500
step = 1
time = 0

Tmin = 2.9
Tmax = 7.0

nmin = 18.9
nmax = 20.6


# Density, Temperature, Velocity
titles = ['Column Density', 'Temperature']
mins = [nmin, Tmin]
maxs = [nmax, Tmax]
cmaps = ['viridis', 'magma']
labels = ['$log_{10}(N_{H} \ [cm^{-2}])$', '$log_{10}(T \ [K])$']


for i in range(istart, iend):

    # print(i)
    fig, axs = plt.subplots(nrows=len(sims), ncols=2, figsize=(6, 4.6))
    fig_color = 'black'
    bg_color = '#DFE6F3'

    for j in range(len(sims)):

        # print(j)
        if cat[j]:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '_slice.h5', 'r') 
        else:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '_slice.h5.0', 'r') 
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

        if cat[j]:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' + str(i) + '_proj.h5', 'r') # open the hdf5 file for reading
        else:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' + str(i) + '_proj.h5.0', 'r') # open the hdf5 file for reading
        head = f.attrs # read the header attributes into a structure, called head
        d = f['d_xy'][:]

        d = d * m_c / (l_c**2)
        n = d / (mu*mp) # number density, particles per cm^3  
        logn = np.log10(n)

        subplots = [logn.T, logT.T]

        for k in range(len(subplots)):

            # print(k)

            im = axs[j][k].imshow(subplots[k], cmap=cmaps[k], vmin=mins[k], vmax = maxs[k]) 
            # axs[j].set_ylabel(labels[j], size=10, color=fig_color)
            axs[j][k].set_xticks(np.linspace(0,nx,9))
            axs[j][k].set_yticks(np.linspace(0,nz,5))
            axs[j][k].invert_yaxis()

            plt.setp(axs[j][k].spines.values(), color=fig_color)
            plt.setp([axs[j][k].get_xticklines(), axs[j][k].get_yticklines()], color=fig_color)

            if j == 0:
                axs[j][k].set_title(titles[k], fontsize=8, color=fig_color)

            if k == 0:
                axs[j][k].set_ylabel(res_labels[j], size=9, rotation='horizontal', ha='right', va='center', color=fig_color)

            if j == (len(sims)-1):
                axs[j][k].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                        labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=6)
                axs[j][k].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.3),1))
                # print(nx*dx)
                [l.set_visible(False) for (i,l) in enumerate(axs[j][k].xaxis.get_ticklabels()) if i % 2 != 0]
                axs[j][k].set_xlabel('$kpc$', size=8, color=fig_color)


                # fig.subplots_adjust(wspace=0.2, hspace=0.1)
                col = []
                for x in range(j+1):
                    col.append(axs[x][k])
                cb = fig.colorbar(im, ax=col, aspect=40, pad=.06)
                cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
                cb.ax.yaxis.set_tick_params(color=fig_color, labelsize=6)
                cb.outline.set_edgecolor(fig_color)
                plt.setp(cbar_yticks, color=fig_color)
                cb.ax.set_ylabel(labels[k], size=8, color=fig_color)
                    
            else:
                axs[j][k].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                        labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.5) #width=0.5 for cooling
    
    # fig.subplots_adjust(wspace=0., hspace=0.05)
    fig.text(0.5, 0.94, str(int(t/t_cc))+r' $t_{cc}$', size=9, color=fig_color)
    # fig.text(0.5, 0.9, str(t)+r' $Myr$', size=8, color=fig_color)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)
