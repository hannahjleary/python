import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6

DE = 0 #Dual Energy Flag

plt.style.use('classic')
plt.rcParams['mathtext.default']='regular'

dnamein='../../../../../ix/eschneider/hjl28/data/radiative/super/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/plots/radiative/super/bowshock/'  # directory where the plot will be saved


sims = ['48/']
labels = ['$R_{48}$']
cat = [True]

vmin = -200.0
vmax = 1200.0

# t_cc = 4.89e4 # (vwind = 10 km/s)
# t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
istart = 0
iend = 150
time = 0
for i in range(istart, iend):

    fig = plt.figure(figsize=(4,7))
    gs = gridspec.GridSpec(2,1,figure=fig, height_ratios=[2,1], hspace=0.1)
    gs_sub = gridspec.GridSpecFromSubplotSpec(3,1,subplot_spec = gs[0], hspace=0)
    ax0 = fig.add_subplot(gs_sub[0])
    ax1 = fig.add_subplot(gs_sub[1])
    ax2 = fig.add_subplot(gs_sub[2])
    ax3 = fig.add_subplot(gs[1])
    fig_color = 'black'
    bg_color = 'white'

    for j in range(len(sims)):
        if cat[j]:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '_slice.h5', 'r') 
        else:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '/' + str(i) + '_slice.h5.0', 'r') 
        head = f.attrs # read the header attributes into a structure, called head

        t  = head['t'] # time of this snapshot, in kyr

        gamma = head['gamma']
        nx = head['dims'][0] # number of cells in the x direction
        ny = head['dims'][1] # number of cells in the y direction
        nz = head['dims'][2] # number of cells in the z direction
        dx = head['dx'][0] # width of cell in x direction

        v_c = head['velocity_unit']
        d_c = head['density_unit']
        e_c = head['energy_unit']
        p_c = e_c
        d = f['d_xy'][:] #the line causing issues
        px  = f['mx_xy'][:]
        py  = f['my_xy'][:]
        pz  = f['mz_xy'][:]
        E = f['E_xy'][:]
        if DE:
            GE = f['GE_xy'][:]
        f.close()

        vx = px/d
        vy = py/d
        vz = pz/d
        if not DE:
            KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
            GE = E - KE
        km = 1e-5
        vx = vx*v_c*km #velocity in the x direction

        n = d * d_c/ (mu*mp) # number density, particles per cm^3  
        T = GE*(gamma-1.0)*p_c / (n*kb) #temperature
        logT = np.log10(T)

        vx_norm = vx/1000
        T_norm = T/1e6 
        n_norm = n/1e-2
        print(np.min(vx_norm[int(ny/2)]))
        plt.suptitle(str(int(t/t_cc))+r' $t_{cc}$')
        # axs[0].plot(np.arange(192), vx_norm[:int(0.375*nx),int(ny/2)], label='$v_{x}$')
        ax0.scatter(np.arange(0.375*nx), logT[:int(0.375*nx),int(ny/2)], color='black', s=0.5)
        ax0.set_ylim(3,9)
        ax0.set_ylabel("$log(T)$", rotation='horizontal', ha='right', va='center', size=8)
        ax0.tick_params(labelbottom=False, labelsize=8)
        ax1.scatter(np.arange(0.375*nx), np.log10(n[:int(0.375*nx),int(ny/2)]), color='black', s=0.5)
        ax1.set_ylim(-3, 2.1)
        ax1.set_ylabel("$log(n)$", rotation='horizontal', ha='right', va='center', size=8)
        ax1.tick_params(labelbottom=False, labelsize=8)
        ax2.scatter(np.arange(0.375*nx), vx_norm[:int(0.375*nx),int(ny/2)], color='black', s=0.5)
        ax2.set_ylim(-.4, 1.3)
        ax2.set_ylabel("$v_{x}/v_{x,w}$", rotation='horizontal', ha='right', va='center', size=8)
        ax2.tick_params(labelsize=8)

        ax0.set_xticks(np.linspace(0,0.375*nx,9))
 
        im = ax3.imshow(vx.T, cmap='magma', vmin=vmin, vmax=vmax) #, vmin=vmin, vmax = vmax
        ax3.set_ylabel(labels[j], size=8, rotation='horizontal', ha='right', va='center', color=fig_color)
        # axs[1].set_xticks(np.linspace(0,nx,9))
        ax3.set_xticks(np.linspace(0,nx,9))
        ax3.set_yticks(np.linspace(0,nz,5))
        ax3.tick_params(labelsize=8)
        # [l.set_visible(False) for (i,l) in enumerate(axs[1].xaxis.get_ticklabels()) if i % 2 != 0]
        # axs[1].invert_yaxis()
        ax3.axhline(ny/2, xmin=0, xmax=192/nx, color='white', zorder=1)

        # print(vx[:][int(5*ny/8)])
        # print("min: " + str(np.min(vx[:][int(5*ny/8)])))
        # print(np.where(np.min(vx[:][int(5*ny/8)])))
        # # axs[j].scatter(np.where(np.min(vx[:][int(5*ny/8)])), 5*ny/8, s=0.2, color='red')
        # print('\n')
        
        # plt.setp(axs[1].spines.values(), color=fig_color)
        # plt.setp([axs[1].get_xticklines(), axs[j].get_yticklines()], color=fig_color)

    #     if j == (len(sims)-1):
    #         axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
    #                 labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=6)
    #         axs[j].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.15),1)) #0.3
    #         [l.set_visible(False) for (i,l) in enumerate(axs[j].xaxis.get_ticklabels()) if i % 2 != 0]
    #         axs[j].set_xlabel('$kpc$', size=6, color=fig_color)
    #     else:
    #         axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
    #                 labelleft=0, labelbottom=0, labeltop=0, labelright=0)

    # cb = fig.colorbar(im, ax=axs.ravel().tolist(), aspect=40, pad=0.025)
    # cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
    # cb.ax.yaxis.set_tick_params(color=fig_color, labelsize=6)
    # cb.outline.set_edgecolor(fig_color)
    # plt.setp(cbar_yticks, color=fig_color)
    # cb.ax.set_ylabel('$kms^{-1}$', size=8, color=fig_color)

    # fig.text(0.65, 0.9, str(int(t/t_cc))+r' $t_{cc}$', size=8, color=fig_color)
    # axs[0].text(5, cells-18, "M ="+str(np.round(float(M), 2)), size=6, color='black')

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.1, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)