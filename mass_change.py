import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

dnamein='../../data/cloud_wind/4_lowest/' # directory where the file is located
dnameout='../../data/cloud_wind/4_lowest/png/' # directory where the plot will be saved

iend = 500
n_step = 10
num = int(iend/n_step)

sims = ['4/', '8/', '16/']
n_init = [np.empty(num), np.empty(num), np.empty(num)]
original_masses = np.empty(3)
masses = [np.empty(num), np.empty(num), np.empty(num)]
cat = [False, False, True]
labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$']

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

        n = d * d_c/(mu*mp) # number density, particles per cm^3  
        if i == 0:
            n_init[j] = n

        cloud_mass = np.nansum(np.where(n >= (n_init[j] / 3.0), n, np.NaN))
        if i == 0:
            original_masses[j] = cloud_mass
        masses[j][int(i/10)] = cloud_mass/ original_masses[j]

        print(str(i) + ": " + sims[j] + "mass: " + str(cloud_mass))

################################# Plots ##########################################

DARK = 0

if DARK:
    fig_color = 'white'
    bg_color = 'grey'
else:
    fig_color = 'black'
    bg_color = 'white'

fig, ax = plt.subplots()

for s in range(len(sims)):
    ax.plot(masses[s], label=labels[s])
ax.legend()
ax.set_xlabel("t $[Myr]$", color=fig_color)
ax.set_ylabel("Change in Cloud Mass $[M/M_{i}]$", color=fig_color)
ax.set_xticks(np.arange(0, 16, 3))
# ax.set_xticklabels([0,2,4,6,8,10])
ax.tick_params(labelsize=8)
# ax.set_xticklabels(n_step*np.arange(0, num))
    # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]

ax.set_facecolor(bg_color)
plt.setp(ax.spines.values(), color=fig_color)
plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

plt.savefig(dnameout + 'M_cloud.png', dpi=300, 
        bbox_inches='tight', pad_inches = 0.2, facecolor = bg_color) #facecolor=bg_color
plt.close(fig)


    