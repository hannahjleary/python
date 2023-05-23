import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K

DE = 0 # Dual Energy Flag

dnamein='../../data/cloud_wind/1/' # directory where the file is located
dnameout='../../data/cloud_wind/1/png/' # directory where the plot will be saved

sims = ['1_256/', '1_512/', '1_1024/']
labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$']
cat = [False, False, True]

t_cc = 4.89e3 # cloud crushing time in kyr
iend = 500
time = 0

for i in range(iend):

    fig, axs = plt.subplots(nrows=len(sims), ncols=1)
    fig_color = 'white'
    bg_color = 'black'

    for j in range(len(sims)):
        if cat[j]:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '_slice.h5', 'r') 
        else:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '_slice.h5.0', 'r') 
        head = f.attrs # read the header attributes into a structure, called head

        t  = head['t'] # time of this snapshot, in kyr

        nx = head['dims'][0] # number of cells in the x direction
        ny = head['dims'][1] # number of cells in the y direction
        nz = head['dims'][2] # number of cells in the z direction
        dx = head['dx'][0] # width of cell in x direction

        v_c = head['velocity_unit']

        r = ny/16

        d = f['d_xy'][:]
        px  = f['mx_xy'][:]

        f.close()

        vx = px/d

        km = 1e-5
        vx = vx*v_c*km #velocity in the x direction

        f.close()

        # print('\t min \t\t\t max')
        # print('Vx: ', np.min(Vx) , '\t' , np.max(Vx))

        vmin = 0.0
        vmax = 100.0

        im = axs[j].imshow(vx.T, cmap='YlOrRd', vmin=vmin, vmax = vmax) 
        axs[j].set_ylabel(labels[j], size=10, rotation='horizontal', ha='right', va='center', color=fig_color)
        axs[j].set_xticks(np.linspace(0,nx,17))
        axs[j].set_yticks(np.linspace(0,nz,5))
        axs[j].invert_yaxis()

        plt.setp(axs[j].spines.values(), color=fig_color)
        plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color=fig_color)

        if j == (len(sims)-1):
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=8)
            axs[j].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.2),1))
            [l.set_visible(False) for (i,l) in enumerate(axs[j].xaxis.get_ticklabels()) if i % 2 != 0]
            axs[j].set_xlabel('$kpc$', size=10, color=fig_color)
        else:
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0)

    cb = fig.colorbar(im, ax=axs.ravel().tolist(), aspect=30)
    cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
    cb.ax.yaxis.set_tick_params(color=fig_color, labelsize=8)
    cb.outline.set_edgecolor(fig_color)
    plt.setp(cbar_yticks, color=fig_color)
    cb.ax.set_ylabel('$kms^{-1}$', size=10, color=fig_color)

    fig.text(0.47, 0.9, str(int(t/t_cc))+r' $t_{cc}$', size=10, color=fig_color)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.3, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)