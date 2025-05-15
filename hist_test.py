import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 0 # Dual Energy Flag

dnamein='../../data/cloud_wind/3/' # directory where the file is located
dnameout='../../data/cloud_wind/3/hists/' # directory where the plot will be saved

istart = 110
iend = 120
n_step = 10
num = int(iend/n_step)

sims = ['4/', '8/', '16/', '32/'] #, '32/', '48/'
cutoffs = np.zeros(len(sims))
cat = [False, False, True, True]
labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$', '$R_{32}$'] #, '$R_{48}$'

DARK = 0

if DARK:
    fig_color = 'white'
    bg_color = 'grey'
else:
    fig_color = 'black'
    bg_color = 'white'

# fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize=(8,7))
# axs = axs.ravel()
    
fig, ax = plt.subplots()

for i in range(istart, iend, n_step):

    for j in range(len(sims)):
        if cat[j]:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '.h5', 'r') 
        else:
            f = h5py.File(dnamein + sims[j] + 'hdf5/' +str(i) + '.h5.0', 'r') 

        head = f.attrs # read the header attributes into a structure, called head

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
        n_init = 1.0 

        vx = px/d
        vy = py/d
        vz = pz/d

        km = 1e-5
        Vx = (px*v_c*km)/d #velocity in the x direction

        dx = 2.4 / nx
        mass = d * dx*dx*dx
        cloud_mass = mass[n >= (n_init / 3.0)]
        mass_tot = np.sum(cloud_mass)
        Vx_cloud = Vx[n >= (n_init / 3.0)]

        Vx_weighted = (Vx_cloud * cloud_mass) / mass_tot

        Vx_weighted = Vx_weighted.ravel()
        print(np.min(Vx_weighted))
        print(np.max(Vx_weighted))

        box_end = n[-1,:,:]

        if cutoffs[j] == 0:
            if np.any(box_end > (n_init/3)):
                cutoffs[j] = i/n_step

        print(str(i) + ": " + sims[j])

################################# Plot ##########################################


        sns.kdeplot(Vx_cloud, bw = 4, ax = ax, shade=True, alpha = 0.3, label=labels[j])
        ax.set_xlim(-25, 150)
        ax.set_xlabel('[km/s]')
        ax.set_ylabel('Density')
        # ax.legend()
        ax.set_title('Velocity Distribution of Cloud Cells')

        # plot hist midpoints and their values


    # for s in range(len(sims)):
    #     ax.plot(velocities_d[s], label=labels[s])
    #     if s == len(sims)-1:
    #         ax.plot(cutoffs[s], velocities_d[s][int(cutoffs[s])], c='black', linestyle=' ', marker='X', label='Mass Leaving Box')
    #     else:
    #         ax.plot(cutoffs[s], velocities_d[s][int(cutoffs[s])], c='black', linestyle=' ', marker='X')
    # ax.legend(loc='lower right', fontsize=10)
    # ax.set_xlabel("t $[Myr]$", color=fig_color, fontsize=12)
    # ax.set_ylabel("v $[kms^{-1}]$", color=fig_color, fontsize=12)
    # ax.set_title('Mean Cloud Velocity', fontsize=12)
    # ax.set_xticks(np.arange(0, 50+1, n_step))
    # ax.set_xticklabels(tick_labels)
    # ax.tick_params(labelsize=9)
    # # ax.set_xticklabels(n_step*np.arange(0, num))
    #     # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]

    # ax.set_facecolor(bg_color)
    # plt.setp(ax.spines.values(), color=fig_color)
    # plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)
