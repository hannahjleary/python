import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1
km = 1e-5

DE = 1 # Dual Energy Flag

dnamein='../../../../../ix/eschneider/hjl28/data/radiative/super/48_2/hdf5/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/plots/radiative/super/vfields/' # directory where the plot will be saved

CAT = 1

# t_cc = 4.89e4 # (vwind = 10 km/s)
# t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)

istart = 0
iend = 500
step = 1

vmin = 3.5
vmax = 6.7


for i in range(istart, iend, step):

    fig, ax = plt.subplots(figsize=(5.5,2))
    fig_color = 'white'
    bg_color = 'black'

    if CAT:
        f = h5py.File(dnamein + str(i) + '_slice.h5', 'r') 
    else:
        f = h5py.File(dnamein + str(i) + '_slice.h5.0', 'r') 
    head = f.attrs # read the header attributes into a structure, called head

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

    if DE:
        GE = f['GE_xy'][:]

    f.close()

    mass = np.sum(d)

    vx = (px/d)*v_c*km
    vy = (py/d)*v_c*km
    vz = (pz/d )*v_c*km
        
    if not DE:   
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE

    n = d * d_c/ (mu*mp) # number density, particles per cm^3  
    T = GE*(gamma-1.0)*p_c / (n*kb) #temperature
    logT = np.log10(T)

    # Sound speed in the wind
    # windT = T[0,0]
    # windn = n[0,0]
    # P = windT * windn * kb / p_c
    # rho = windn * mu * mp / d_c
    # c = np.sqrt((gamma-1.0)*P/ rho)
    # c = c * v_c * 1e-5

    fs = 10
    ls = 8

    print(np.min(vx))
    x = np.arange(0, nx, 1)
    y = np.arange(0,ny,1)
    X, Y = np.meshgrid(x,y, indexing='ij')
    skip=60
    im = ax.imshow(logT.T, cmap='magma', vmin=vmin, vmax = vmax) #, vmin=vmin, vmax = vmax
    # axs[j].set_ylabel(labels[j], size=fs, rotation='horizontal', ha='right', va='center', color=fig_color)
    ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], vx[::skip, ::skip], vy[::skip, ::skip], scale=95000, width=0.0015, headwidth=3, headaxislength=4,  color='white')
    ax.set_xticks(np.linspace(0,nx,9))
    ax.set_yticks(np.linspace(0,nz,5))
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='both', direction='in', color=fig_color, bottom=1, left=1, top=1, right=1, 
            labelleft=0, labelbottom=0, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=ls)

    plt.setp(ax.spines.values(), color=fig_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

    divider = make_axes_locatable(ax)
    cbax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(im, cax = cbax)
    cbax.tick_params(axis='y', direction='in', color=fig_color, labelcolor=fig_color)
    cb.solids.set_edgecolor('face')
    cb.outline.set_edgecolor(fig_color)
    cbax.set_ylabel(r'$\mathrm{log}_{10}(T)$ [K]', color=fig_color)

    fig.text(0.48, 0.9, str(int(t/t_cc))+r' $t_{cc}$', size=fs, color=fig_color)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
                bbox_inches='tight', facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)