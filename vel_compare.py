import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 0 # Dual Energy Flag
CAT = 1 # Conactenated Flag

slower_dir='../../data/cloud_wind/3/16/hdf5/' # directory where the file is located
faster_dir='../../data/cloud_wind/4_high/16/hdf5/' # directory where the file is located
dnameout='../../plots/velocity_compare/' # directory where the plot will be saved
time = 250
# t_cc = 4.89e4 # (vwind = 10 km/s)
t_cc_slower = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
t_cc_faster = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
iend = 500

#####################################################################
for i in range(iend):
    f = h5py.File(slower_dir + '0_slice.h5', 'r') 
    head = f.attrs # read the header attributes into a structure, called head

    gamma = head['gamma']
    nx = head['dims'][0] # number of cells in the x direction
    ny = head['dims'][1] # number of cells in the y direction
    nz = head['dims'][2] # number of cells in the z direction
    dx = head['dx'][0] # width of cell in x direction

    v_c = head['velocity_unit']
    d_c = head['density_unit']
    e_c = head['energy_unit']
    m_c = head['mass_unit']
    l_c = head['length_unit']
    p_c = e_c

    d = f['d_xy'][:]
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

    windT = T[0,0]
    windn = n[0,0]
    windv = vx[0,0]
    cells = ny
    P = windT * windn * kb / p_c
    rho = windn * mu * mp / d_c
    c = np.sqrt(gamma*P/ rho)
    c = c * v_c * 1e-5
    M_slower = windv / c

    f = h5py.File(slower_dir + str(i) + '_proj.h5', 'r') 
    head = f.attrs # read the header attributes into a structure, called head
    t_slower = head['t']
    d = f['d_xy'][:]
    d = d * m_c / (l_c**2)
    n = d / (mu*mp) # number density, particles per cm^3  
    logn_slower= np.log10(n)


    #############################################################
    f = h5py.File(faster_dir + '0_slice.h5', 'r') 
    head = f.attrs # read the header attributes into a structure, called head

    gamma = head['gamma']
    nx = head['dims'][0] # number of cells in the x direction
    ny = head['dims'][1] # number of cells in the y direction
    nz = head['dims'][2] # number of cells in the z direction
    dx = head['dx'][0] # width of cell in x direction

    v_c = head['velocity_unit']
    d_c = head['density_unit']
    e_c = head['energy_unit']
    m_c = head['mass_unit']
    l_c = head['length_unit']
    p_c = e_c

    d = f['d_xy'][:]
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

    windT = T[0,0]
    windn = n[0,0]
    windv = vx[0,0]
    cells = ny
    P = windT * windn * kb / p_c
    rho = windn * mu * mp / d_c
    c = np.sqrt(gamma*P/ rho)
    c = c * v_c * 1e-5
    M_faster = windv / c

    f = h5py.File(faster_dir + str(i) + '_proj.h5', 'r') 
    head = f.attrs # read the header attributes into a structure, called head
    t_faster = head['t']
    d = f['d_xy'][:]
    d = d * m_c / (l_c**2)
    n = d / (mu*mp) # number density, particles per cm^3  
    logn_faster= np.log10(n)

    ###############################################################################

    units = '$log_{10}(N_{H} [cm^{-2}])$'

    fig, axs = plt.subplots(nrows=2, ncols=1)
    fig_color = 'white'
    bg_color = 'black'

    im = axs[0].imshow(logn_slower.T, cmap='viridis', vmin=19.6, vmax=21.0) 
    # axs[0].text(0.03*nx,0.85*ny,'adiabatic', size=10, color='white')
    axs[0].text(0.03*nx, 0.85*ny, "M = "+str(np.round(float(M_slower), 2)), size=10, color='white')
    axs[0].text(0.85*nx,0.85*ny, str(int(t_slower/t_cc_slower)) + r' $t_{cc}$', size=10, color='white')
    axs[0].set_xticks(np.linspace(0,nx,9))
    axs[0].set_yticks(np.linspace(0,nz,5))
    axs[0].invert_yaxis()
    axs[0].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
            labelleft=0, labelbottom=0, labeltop=0, labelright=0)


    axs[1].imshow(logn_faster.T, cmap='viridis', vmin=19.6, vmax=21.0) 
    # axs[1].text(0.03*nx,0.85*ny,'cooling', size=10, color='white')
    axs[1].text(0.03*nx, 0.85*ny, "M = "+str(np.round(float(M_faster), 2)), size=10, color='white')
    axs[1].text(0.85*nx,0.85*ny, str(int(t_faster/t_cc_faster)) + r' $t_{cc}$', size=10, color='white')
    axs[1].set_xticks(np.linspace(0,nx,9))
    axs[1].set_yticks(np.linspace(0,nz,5))
    axs[1].invert_yaxis()
    axs[1].tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
            labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=8)
    axs[1].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.3),1))
    [l.set_visible(False) for (i,l) in enumerate(axs[1].xaxis.get_ticklabels()) if i % 2 != 0]
    axs[1].set_xlabel('$kpc$', size=10, color=fig_color)

    for j in range(2):
        plt.setp(axs[j].spines.values(), color=fig_color)
        plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color=fig_color)


    fig.subplots_adjust(wspace=0.3, hspace=0.1)

    cb = fig.colorbar(im, ax=axs.ravel().tolist(), aspect=40, pad=.03)
    cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
    cb.ax.yaxis.set_tick_params(color=fig_color, labelsize=8)
    cb.outline.set_edgecolor(fig_color)
    plt.setp(cbar_yticks, color=fig_color)
    cb.ax.set_ylabel(units, size=10, color=fig_color) 

    # fig.subplots_adjust(wspace=0.3, hspace=0.1)

    # fig.text(0.5, 0.9, str(int(t))+r' $t_{cc}$', size=10, color=fig_color)
    plt.suptitle('Density Projections', color=fig_color,fontsize = 12, y=.93)
    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, transparent=True) #facecolor=bg_color
    plt.close(fig)