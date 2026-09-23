import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1
SOLAR_METAL_MASS_FRAC = 0.01295

DE = 1 # Dual Energy Flag
CT = 0 #Cloud Tracking Flag
CAT = 0

dnamein='../../../../../ix/eschneider/hjl28/data/tests/dev-test/hdf5/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/data/tests/dev-test/png/'

t_cc = 4.89e4 # (vwind = 10 km/s)
#t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
istart = 0
iend = 10
time = 0

Tmin = 3.5
Tmax = 6.5

nmin = 19.0
nmax = 20.2

vmin = -20 #-200
vmax = 120 #1200

Zmin = -0.0
Zmax = 1.4

if CT:
    velocity_shifts = []

    with open(dnamein + "cloud_velocities.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            velocity_shifts.append(float(line))

for i in range(istart, iend):

    print(str(i))
    
    if CAT:
        f = h5py.File(dnamein + str(i) + '/' + str(i) + '_slice.h5', 'r') # open the hdf5 file for reading
    else:
        f = h5py.File(dnamein + str(i) + '/' + str(i) + '_slice.h5.0', 'r') # open the hdf5 file for reading
    head = f.attrs # read the header attributes into a structure, called head
    # print(list(f.keys()))

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
    d_metals = f['metal_density_xy'][:]
    # scalar = f['scalar0_xy'][:]

    if DE:
        GE = f['GE_xy'][:]

    f.close()

    Z = d_metals / d
    Z_sol = Z / SOLAR_METAL_MASS_FRAC

    # print(gamma)

    # scalar = scalar / d

    n = d * d_c/ (mu*mp) # number density, particles per cm^3  

    vx = px/d
    vy = py/d
    vz = pz/d
        
    if not DE:    
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE

    T = GE*(gamma-1.0)*p_c / (n*kb) #temperature
    logT = np.log10(T)

    print(np.min(T))

    km = 1e-5

    Px = px * v_c * km * d_c

    Vx = vx*v_c*km #velocity in the x direction
    
    if CT:
        Vx = Vx + velocity_shifts[i]
        print(str(Vx[10, int(ny/2)]))

    if CAT:
        f = h5py.File(dnamein + str(i) + '/' + str(i) + '_proj.h5', 'r') # open the hdf5 file for reading
    else:
        f = h5py.File(dnamein + str(i) + '/' + str(i) + '_proj.h5.0', 'r') # open the hdf5 file for reading
    head = f.attrs # read the header attributes into a structure, called head
    d  = f['d_xy'][:]
    # T = f['T_xy'][:]    


    d = d * m_c / (l_c**2)
    n = d / (mu*mp) # number density, particles per cm^3  
    logn = np.log10(n)

    logv = np.log10(Vx)
    P = np.log10(n*kb*T)

    f.close()

    # print(d_c)
    # print(p_c)
    # print('\t min \t\t\t max')
    # print('n: ', np.min(logn) , '\t' , np.max(logn))
    # print('T: ', np.min(logT) , '\t' , np.max(logT))
    # print('Vx: ', np.min(Vx) , '\t' , np.max(Vx))

    subplots = [logT.T, logn.T, Vx.T, Z_sol.T] #subplots = [logT.T, P.T, Vx.T]
    mins = [Tmin, nmin, vmin, Zmin]
    maxs = [Tmax, nmax, vmax, Zmax]
    cmaps = ['plasma', 'viridis', 'magma_r', 'cividis']
    labels = ['$log_{10}(K)$', '$log_{10}(N_{H})$ [$cm^{-2}$]', '$kms^{-1}$', '$Z/Z_{\odot}$']

    fig, axs = plt.subplots(nrows=len(subplots), ncols=1)
    fig_color = 'white'
    bg_color = 'black'

    for j in range(len(subplots)):

        im = axs[j].imshow(subplots[j], cmap=cmaps[j], vmin=mins[j], vmax = maxs[j]) #, vmin=mins[j], vmax = maxs[j]
        # axs[j].set_ylabel(labels[j], size=10, color=fig_color)
        axs[j].set_xticks(np.linspace(0,nx,9))
        axs[j].set_yticks(np.linspace(0,nz,5))
        axs[j].invert_yaxis()

        plt.setp(axs[j].spines.values(), color=fig_color)
        plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color=fig_color)

        if j == (len(subplots)-1):
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=6)
            axs[j].set_xticklabels(np.round(np.linspace(0,nx*dx+.01, 9),1))
            # print(nx*dx)
            [l.set_visible(False) for (i,l) in enumerate(axs[j].xaxis.get_ticklabels()) if i % 2 != 0]
            axs[j].set_xlabel('$kpc$', size=8, color=fig_color)
        else:
            axs[j].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0)
            
        divider = make_axes_locatable(axs[j])
        cax = divider.append_axes('right', size = 0.10, pad = 0.17)
        cb = plt.colorbar(im, cax=cax)
        cb.set_ticks(np.round(np.linspace(mins[j], maxs[j], 5), 2))
        cax.tick_params(axis='y', direction='out', color = fig_color, labelcolor=fig_color, labelsize=6)
        cax.set_ylabel(labels[j], size=8, color=fig_color)
        cb.outline.set_edgecolor(fig_color)

    # fig.text(0.5, 0.9, str(int(t/t_cc))+r' $t_{cc}$', size=8, color=fig_color)
    fig.text(0.5, 0.9, str(t)+r' $kyr$', size=8, color=fig_color)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)
