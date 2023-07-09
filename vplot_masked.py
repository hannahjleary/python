import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K

DE = 0 # Dual Energy Flag
T_CLOUD_PLOT = 0
D_CLOUD_PLOT = 1

dnamein='../../data/cloud_wind/4.1/' # directory where the file is located
dnameout='../../data/cloud_wind/4.1/png/' # directory where the plot will be saved

# t_cc = 4.89e3 # cloud crushing time in kyr
iend = 500
n_step = 10
num = int(iend/n_step)

sims = ['R4/', 'R8/', 'R16/']
if T_CLOUD_PLOT:
    velocities_T = [np.empty(num), np.empty(num), np.empty(num)]
if D_CLOUD_PLOT:
    velocities_d = [np.empty(num), np.empty(num), np.empty(num)]
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

        # print(j)
        # print(np.min(T))
        # print(np.max(T))
        # print()
        print(str(i) + ": " + sims[j])

        if T_CLOUD_PLOT:
            clouds = np.where(T<4e4, Vx, np.NaN)
            Vx_avg = np.nanmean(clouds)
            velocities_T[j][int(i/10)] = Vx_avg

        if D_CLOUD_PLOT:
            n_init = 1.0
            clouds = np.where(n >= (n_init / 3.0), Vx, np.NaN)
            Vx_avg = np.nanmean(clouds)
            velocities_d[j][int(i/10)] = Vx_avg

################################# Plots ##########################################

DARK = 0

if DARK:
    fig_color = 'white'
    bg_color = 'grey'
else:
    fig_color = 'black'
    bg_color = 'white'

times = np.arange(0,iend,n_step)

############ Cloud = < 2e4 K ########################

if T_CLOUD_PLOT:
    fig, ax = plt.subplots()

    for s in range(len(sims)):
        ax.plot(velocities_T[s], label=labels[s])
    ax.legend()
    ax.set_xlabel("t $[Myr]$", color=fig_color)
    ax.set_ylabel("Mean Cloud Velocity $[kms^{-1}]$", color=fig_color)
    ax.set_xticks(np.arange(0, num+1, n_step))
    ax.tick_params(labelsize=8)
    # ax.set_xticklabels(n_step*np.arange(0, num))
        # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]

    ax.set_facecolor(bg_color)
    plt.setp(ax.spines.values(), color=fig_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

    plt.savefig(dnameout + 'Vx_Tcloud_4.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, facecolor = bg_color) #facecolor=bg_color
    plt.close(fig)

############ Cloud = > 1/3 init. density ########################

if D_CLOUD_PLOT:
    fig, ax = plt.subplots()

    for s in range(len(sims)):
        ax.plot(velocities_d[s], label=labels[s])
    ax.legend()
    ax.set_xlabel("t $[Myr]$", color=fig_color)
    ax.set_ylabel("Mean Cloud Velocity $[kms^{-1}]$", color=fig_color)
    ax.set_xticks(np.linspace(0,num,11))
    # ax.set_xticks(np.arange(0, num+1, n_step))
    ax.set_xticklabels(np.arange(0,10.1,1.0))
    ax.tick_params(labelsize=8)
    # ax.set_xticklabels(n_step*np.arange(0, num))
        # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]

    ax.set_facecolor(bg_color)
    plt.setp(ax.spines.values(), color=fig_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

    plt.savefig(dnameout + 'Vx_dcloud.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, facecolor = bg_color) #facecolor=bg_color
    plt.close(fig)

    