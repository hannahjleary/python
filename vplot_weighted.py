import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K

DE = 0 # Dual Energy Flag

dnamein='../../data/cloud_wind/2/' # directory where the file is located
dnameout='../../data/cloud_wind/2/png/' # directory where the plot will be saved

t_cc = 4.89e3 # cloud crushing time in kyr
iend = 500
n_step = 10
n = int(iend/n_step)

sims = ['128/', '256/', '512/']
v_128 = np.empty(n)
v_256 = np.empty(n)
v_512 = np.empty(n)
velocities = [v_128, v_256, v_512]
# labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$']
cat = [False, False, True]

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
        v_c = head['velocity_unit']
        e_c = head['energy_unit']
        p_c = e_c # pressure units are the same as energy density units, density*velocity^2/length^3

        d  = f['density'][:]
        px  = f['momentum_x'][:]
        py  = f['momentum_y'][:]
        pz  = f['momentum_z'][:]
        E = f['Energy'][:]

        f.close()

        if DE:
            GE = f['GasEnergy'][:]

        mu = 1.0 # mean molecular weight (mu) of 1
        n = d * d_c/(mu*mp) # number density, particles per cm^3  

        vx = px/d
        vy = py/d
        vz = pz/d
            
        if not DE:    
            KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
            GE = E - KE

        T = GE*(gamma-1.0)*p_c / (n*kb) #temperature
        # print(T)

        f.close()

        km = 1e-5
        Vx = (px*v_c*km)/d #velocity in the x direction
        Vx_clouds = np.where(T<4e4, Vx, np.NaN)
        Vx_avg = np.nanmean(Vx_clouds)
        velocities[j][int(i/10)] = Vx_avg

        ###########################################################################
dark = 0

if dark:
    fig_color = 'white'
    bg_color = 'grey'
else:
    fig_color = 'black'
    bg_color = 'white'

fig, ax = plt.subplots()

times = np.arange(0,iend,n_step)
labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$']

for s in range(len(sims)):
    ax.plot(times, velocities[s], label=labels[s])
ax.legend()

ax.set_xlabel("$kyr$", color=fig_color)
ax.set_ylabel("$kms^{-1}$", color=fig_color)

ax.set_facecolor(bg_color)
plt.setp(ax.spines.values(), color=fig_color)
plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

# ax.set_xticks(times)
# ax.set_yticks(np.linspace(0,nz,5))
[l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]


plt.savefig(dnameout + 'v_vs_t.png', dpi=300, 
        bbox_inches='tight', pad_inches = 0.2, facecolor = bg_color) #facecolor=bg_color
plt.close(fig)
