import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 0.6

DE = 1 # Dual Energy Flag
CAT = 1 # Conactenated Flag

adiabatic_dir='../../data/cloud_wind/1/16/hdf5/' # directory where the file is located
cooling_dir='../../data/cloud_wind/3/16/hdf5/' # directory where the file is located
dnameout='../../plots/compare2/' # directory where the plot will be saved

# t_cc = 4.89e4 # (vwind = 10 km/s)
t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
iend = 10

TEMP = 0
DENS = 1
VEL = 0
adiabatic = 0
cooling = 0
plot = 0
vals = [adiabatic, cooling]
directories = [adiabatic_dir, cooling_dir]
labels = ["adiabatic", "radiative cooling"]
units = ['$log_{10}(K)$', '$log_{10}(N_{H})$ [$cm^{-2}$]', '$kms^{-1}$']

for i in range(iend):

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    fig_color = 'black'
    bg_color = 'white'

    for j in range(2):
        if TEMP or VEL:
            if CAT:
                f = h5py.File(directories[j] + str(i) + '_slice.h5', 'r') 
            else:
                f = h5py.File(directories[j] + str(i) + '_slice.h5.0', 'r') 
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
            vx = vx*v_c*km #velocity in the x direction

            if TEMP:
                plot = logT
                units = '$log_{10}(K)$'
            else:
                plot = vx
                units = '$km s^{-1}$'

            if j==1:
                plot = plot[0:2*nx//3,:,:]

        else:
            if CAT:
                f = h5py.File(directories[j] + str(i) + '_proj.h5', 'r') 
            else:
                f = h5py.File(directories[j] + str(i) + '_proj.h5.0', 'r') 
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
            d  = f['d_xy'][:]

            f.close()

            d = d * m_c / (l_c**2)
            n = d / (mu*mp) # number density, particles per cm^3  
            logn = np.log10(n)

            plot = logn
            units = units = '$log_{10}(N_{H} \ [cm^{-2}])$'

            if j==1:
                plot = plot[0:2*nx//3,:,:]

        im = axs[j].imshow(plot.T, cmap='viridis', vmin=19.6, vmax=21) 
        axs[j].text(0.03*nx,0.85*ny,labels[j], size=10, color='white')
        # axs[j].text(0.85*nx,0.85*ny, str(int(t/1000)) + ' Myr', size=10, color='white') 
        axs[j].text(0.85*nx,0.85*ny, str(int(t/t_cc))+r' $t_{cc}$', size=10, color='white')
        axs[j].set_xticks(np.linspace(0,nx,9))
        axs[j].set_yticks(np.linspace(0,nz,5))
        axs[j].invert_yaxis()

        plt.setp(axs[j].spines.values(), color=fig_color)
        plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color=fig_color)

        if j == (len(vals)-1):
            axs[j].tick_params(axis='both', which='both', direction='in', color='white', bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=8)
            axs[j].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.3),1))
            [l.set_visible(False) for (i,l) in enumerate(axs[j].xaxis.get_ticklabels()) if i % 2 != 0]
            axs[j].set_xlabel('$kpc$', size=10, color=fig_color)
        else:
            axs[j].tick_params(axis='both', which='both', direction='in', color='white', bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0)

        # cb = plt.colorbar(im, ax=axs[j], aspect=30, pad=0.02)
    fig.subplots_adjust(wspace=0.3, hspace=0.1)

    cb = fig.colorbar(im, ax=axs.ravel().tolist(), aspect=40, pad=.03)
    cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
    cb.ax.yaxis.set_tick_params(color=fig_color, labelsize=8)
    cb.outline.set_edgecolor(fig_color)
    plt.setp(cbar_yticks, color=fig_color)
    cb.ax.set_ylabel(units, size=10, color=fig_color) 

    # fig.subplots_adjust(wspace=0.3, hspace=0.1)

    # fig.text(0.5, 0.9, str(int(t))+r' $t_{cc}$', size=10, color=fig_color)
    plt.suptitle('Column Density', color=fig_color,fontsize = 10, y=.93)
    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color or transparent=True
    plt.close(fig)