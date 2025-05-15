import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 1.0 # mean molecular weight (mu) of 1

dnamein='../../../../../ix/data/' # directory where the file is located
dnameout='../../../../../ix/plots/' # directory where the plot will be saved

iend = 500
n_step = 10
num = int(iend/n_step)

sims = ['4/', '8/', '16/', '32/', '48/']
cat = [False, False, True, True, True]
labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$', '$R_{32}$', '$R_{48}$']
tick_labels = np.arange(0, 50+1, 10)

masses = [np.empty(num), np.empty(num), np.empty(num), np.empty(num), np.empty(num)]
cutoffs = np.zeros(len(sims))

for i in range(0, iend, n_step):

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

        d  = f['density'][:]

        f.close()
        dx = 1.6 / nx
        mass = d * dx*dx*dx

        n = d * d_c/(mu*mp) # number density, particles per cm^3  
        n_init = 1.0
        cloud_mass = mass[n >= n_init/3]
        mass_tot = np.sum(cloud_mass)
        if i == 0 and j == 0:
            mass_init = mass_tot
        masses[j][int(i/n_step)] = mass_tot / mass_init

        box_end = n[-1,:,:]

        if cutoffs[j] == 0:
            if np.any(box_end > (n_init/3)):
                cutoffs[j] = i/n_step
                print(str(i/n_step))

        print(str(i) + ": " + sims[j] + " mass: " + str(mass_tot / mass_init))

################################# Plots ##########################################

DARK = 0

if DARK:
    fig_color = 'white'
    bg_color = 'grey'
else:
    fig_color = 'black'
    bg_color = 'white'

fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(14,3))

for s in range(len(sims)):
    ax.plot(masses[s], label=labels[s])
    if s == len(sims)-1:
      ax.plot(cutoffs[s], masses[s][int(cutoffs[s])], linestyle=' ', c='black', marker='X', label='Mass Leaving Box')
    else:
        ax.plot(cutoffs[s], masses[s][int(cutoffs[s])], linestyle=' ', c='black', marker='X')
    # print(str(cutoffs[s]) + " , " + str(masses[s][int(cutoffs[s])]))
    
ax.legend(loc = 'upper right', fontsize=10)
ax.set_xlabel("t $[Myr]$", color=fig_color, fontsize=12)
ax.set_ylabel("$[M/M_{i}]$", color=fig_color, fontsize=12)
ax.set_title('Change in Cloud Mass', color=fig_color, fontsize=12)
ax.set_xticks(np.arange(0, 50+1, n_step))
ax.set_xticklabels(tick_labels)
ax.tick_params(labelsize=9)
ax.set_facecolor(bg_color)
plt.setp(ax.spines.values(), color=fig_color)
plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)
# ax.axvline(x=20, ls='--', c='black', alpha=0.5)

plt.savefig(dnameout + 'M_cloud.png', dpi=300, 
        bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
plt.close(fig)


    