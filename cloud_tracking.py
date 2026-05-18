import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 1 # Dual Energy Flag
CREATE_VEL_FILE = 0
CAT = 0

dnamein='../../../../../ix/eschneider/hjl28/data/tests/cloud_tracking/' # directory where the file is located
dname_reg = dnamein + "hdf5_sub_shock/raw/"
dname_ct = dnamein + "hdf5_large_ct/raw/"
dnameout='../../../../../ix/eschneider/hjl28/plots/tests/cloud_tracking/png_sub_shock/'

velocity_shifts = []
if CREATE_VEL_FILE:
    with open(dnamein + "/ct_output.out", "r") as f:
        lines = f.readlines()

    for line in lines:
        if "Average cloud velocity" in line:
            match = re.search(r"Average cloud velocity = ([0-9eE+\-.]+) km/s", line)
            if match:
                velocity_shifts.append(float(match.group(1)))

        with open(dnamein + "hdf5_large_ct/cloud_velocities.txt", "w") as f:
            for v in velocity_shifts:
                f.write(str(v) + "\n")
else:
    with open(dnamein + "hdf5_large_ct/cloud_velocities.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            velocity_shifts.append(float(line))

print(f"Found {len(velocity_shifts)} velocity values")

# t_cc = 4.89e4 # (vwind = 10 km/s)
t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
istart = 0
iend = 300
time = 0

Tmin = 3.5
Tmax = 6.5

nmin = 19.1
nmax = 20.8

vmin = -25 #-200
vmax = 120 #1200

for i in range(istart, iend):

    print(str(i))
    
    if CAT:
        slice = h5py.File(dname_reg + str(i) + '/' + str(i) + '_slice.h5', 'r') # open the hdf5 file for reading
        proj = h5py.File(dname_reg + str(i) + '/' + str(i) + '_proj.h5', 'r') # open the hdf5 file for reading
    else:
        slice = h5py.File(dname_reg + str(i) + '/' + str(i) + '_slice.h5.0', 'r') # open the hdf5 file for reading
        proj = h5py.File(dname_reg + str(i) + '/' + str(i) + '_proj.h5.0', 'r') # open the hdf5 file for reading
    slice_head = slice.attrs # read the header attributes into a structure, called head
    proj_head = proj.attrs

    gamma = slice_head['gamma'] # ratio of specific heats
    t  = slice_head['t'] # time of this snapshot, in kyr
    nx = slice_head['dims'][0] # number of cells in the x direction
    ny = slice_head['dims'][1] # number of cells in the y direction
    nz = slice_head['dims'][2] # number of cells in the z direction
    dx = slice_head['dx'][0] # width of cell in x direction
    l_c = slice_head['length_unit']
    t_c = slice_head['time_unit']
    m_c = slice_head['mass_unit']
    d_c = slice_head['density_unit']                 
    v_c = slice_head['velocity_unit']
    e_c = slice_head['energy_unit']
    p_c = e_c # pressure units are the same as energy density units, density*velocity^2/length^3

    d  = proj['d_xy'][:]
    px  = slice['mx_xy'][:]
    py  = slice['my_xy'][:]
    pz  = slice['mz_xy'][:]
    E = slice['E_xy'][:]

    if DE:
        GE = slice['GE_xy'][:]

    f.close()

    # print(gamma)

    n_reg = d * m_c/ (l_c**2 * mu*mp) # number density, particles per cm^3  
    logn_reg = np.log10(n_reg)

    vx = px/d
    vy = py/d
    vz = pz/d
        
    if not DE:    
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE

    T = GE*(gamma-1.0)*p_c / (n_reg*kb) #temperature
    logT_reg = np.log10(T)

    km = 1e-5

    Px = px * v_c * km * d_c

    Vx_reg = vx*v_c*km #velocity in the x direction

    if CAT:
        slice = h5py.File(dname_ct + str(i) + '/' + str(i) + '_slice.h5', 'r') # open the hdf5 file for reading
        proj = h5py.File(dname_ct + str(i) + '/' + str(i) + '_proj.h5', 'r') # open the hdf5 file for reading
    else:
        slice = h5py.File(dname_ct + str(i) + '/' + str(i) + '_slice.h5.0', 'r') # open the hdf5 file for reading
        proj = h5py.File(dname_ct + str(i) + '/' + str(i) + '_proj.h5.0', 'r') # open the hdf5 file for reading
    slice_head = slice.attrs # read the header attributes into a structure, called head
    proj_head = proj.attrs

    d  = proj['d_xy'][:]
    px  = slice['mx_xy'][:]
    py  = slice['my_xy'][:]
    pz  = slice['mz_xy'][:]
    E = slice['E_xy'][:]

    if DE:
        GE = slice['GE_xy'][:]

    f.close()

    n_ct = d * m_c/ (l_c**2 * mu*mp) # number density, particles per cm^3  
    logn_ct = np.log10(n_ct)

    vx = px/d
    vy = py/d
    vz = pz/d
        
    if not DE:    
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE

    T = GE*(gamma-1.0)*p_c / (n_ct*kb) #temperature
    logT_ct = np.log10(T)

    km = 1e-5

    Px = px * v_c * km * d_c

    Vx_ct = vx*v_c*km #velocity in the x direction

    Vx_ct = vx*v_c*km #velocity in the x direction
    Vx_ct = Vx_ct + velocity_shifts[i]

    # subplots = [logT.T, logn.T, Vx.T] #subplots = [logT.T, P.T, Vx.T]
    # mins = [Tmin, nmin, vmin]
    # maxs = [Tmax, nmax, vmax]
    cmaps = ['plasma', 'viridis', 'magma_r']
    sns.set_palette("mako")
    labels = ['$log_{10}(K)$', '$log_{10}(N_{H})$ [$cm^{-2}$]', '$kms^{-1}$']

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8.3,3))
    fig_color = 'white'
    bg_color = 'black'

    im1 = axs[0][0].imshow(logn_reg.T, cmap='mako', vmin=nmin, vmax=nmax) #, vmin=mins[j], vmax = maxs[j]
    axs[0][0].set_xticks(np.linspace(0,nx,9))
    axs[0][0].set_yticks(np.linspace(0,nz,5))
    axs[0][0].invert_yaxis()
    axs[0][0].set_title('No Cloud Tracking', fontsize=8, color=fig_color)

    plt.setp(axs[0][0].spines.values(), color=fig_color)
    plt.setp([axs[0][0].get_xticklines(), axs[0][0].get_yticklines()], color=fig_color)
    axs[0][0].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0)

    im2 = axs[1][0].imshow(Vx_reg.T, cmap='rocket_r', vmin=vmin, vmax = vmax) #, vmin=mins[j], vmax = maxs[j]
    axs[1][0].set_xticks(np.linspace(0,nx,9))
    axs[1][0].set_yticks(np.linspace(0,nz,5))
    axs[1][0].invert_yaxis()

    plt.setp(axs[1][0].spines.values(), color=fig_color)
    plt.setp([axs[1][0].get_xticklines(), axs[1][0].get_yticklines()], color=fig_color)
    axs[1][0].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=6)
    axs[1][0].set_xticklabels(np.round(np.linspace(0,nx*dx+.01, 8),1))
    print(nx*dx)
    [l.set_visible(False) for (i,l) in enumerate(axs[1][0].xaxis.get_ticklabels()) if i % 2 != 0]
    axs[1][0].set_xlabel('$kpc$', size=8, color=fig_color)

    im3 = axs[0][1].imshow(logn_ct.T, cmap='mako', vmin=nmin, vmax = nmax) #, vmin=mins[j], vmax = maxs[j]
    axs[0][1].set_xticks(np.linspace(0,nx,9))
    axs[0][1].set_yticks(np.linspace(0,nz,5))
    axs[0][1].invert_yaxis()
    axs[0][1].set_title('Cloud Tracking', fontsize=8, color=fig_color)

    plt.setp(axs[0][1].spines.values(), color=fig_color)
    plt.setp([axs[0][1].get_xticklines(), axs[0][1].get_yticklines()], color=fig_color)
    axs[0][1].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0)

    im4 = axs[1][1].imshow(Vx_ct.T, cmap='rocket_r', vmin=vmin, vmax = vmax) #, vmin=mins[j], vmax = maxs[j]
    axs[1][1].set_xticks(np.linspace(0,nx,9))
    axs[1][1].set_yticks(np.linspace(0,nz,5))
    axs[1][1].invert_yaxis()

    plt.setp(axs[1][1].spines.values(), color=fig_color)
    plt.setp([axs[1][1].get_xticklines(), axs[1][1].get_yticklines()], color=fig_color)
    axs[1][1].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=6)
    axs[1][1].set_xticklabels(np.round(np.linspace(0,nx*dx+.01, 8),1))
            # print(nx*dx)
    [l.set_visible(False) for (i,l) in enumerate(axs[1][1].xaxis.get_ticklabels()) if i % 2 != 0]
    axs[1][1].set_xlabel('$kpc$', size=8, color=fig_color)

            
    divider = make_axes_locatable(axs[0][1])
    cax = divider.append_axes('right', size = 0.10, pad = 0.17)
    cb1 = plt.colorbar(im3, cax=cax)
    cb1.set_ticks(np.round(np.linspace(nmin, nmax, 5), 2))
    cax.tick_params(axis='y', direction='out', color = fig_color, labelcolor=fig_color, labelsize=6)
    cax.set_ylabel('$log_{10}(N_{H})$ [$cm^{-2}$]', size=8, color=fig_color)
    cb1.outline.set_edgecolor(fig_color)

    divider = make_axes_locatable(axs[1][1])
    cax = divider.append_axes('right', size = 0.10, pad = 0.17)
    cb2 = plt.colorbar(im4, cax=cax)
    cb2.set_ticks(np.round(np.linspace(vmin, vmax, 5), 2))
    cax.tick_params(axis='y', direction='out', color = fig_color, labelcolor=fig_color, labelsize=6)
    cax.set_ylabel('$kms^{-1}$', size=8, color=fig_color)
    cb2.outline.set_edgecolor(fig_color)

    # fig.text(0.5, 0.9, str(int(t/t_cc))+r' $t_{cc}$', size=8, color=fig_color)
    # fig.text(0.5, 0.9, str(int(t/1000))+r' $Myr$', size=8, color=fig_color)
    fig.suptitle(str(int(t/1000))+r' $Myr$', color=fig_color, fontsize=8)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)