import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable 

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K

istart=0
iend = 10
dnamein='../../../../../ix/eschneider/hjl28/data/KH/mixinglayer/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/plots/KH/'

#plot = input("Enter 'd' 'P' or 'T': ")

for i in range(istart, iend):

    f = h5py.File(dnamein+str(i)+'/' + str(i) + '.h5.0', 'r')

    head = f.attrs

    #print(f.keys())

    gamma = head['gamma'] #ratio of specific heats
    t = head['t'] #time of snapshot (kyr)
    nx = head['dims'][0] # number of cels in he x direction
    ny = head['dims'][1] # number of cells in the y direction
    nz = head['dims'][2] # number of cells in the z direction
    dx = head['dx'][0] #width of cells in the x direction
    dy = head['dx'][1] #width of cells in the y direction
    dz = head['dx'][2] #width of cells in the z direction
    l_c = head['length_unit']
    t_c = head['time_unit']
    m_c = head['mass_unit']
    d_c = head['density_unit']
    v_c = head['velocity_unit']
    e_c = head['energy_unit']
    p_c = e_c

    d  = f['density'][:]
    px  = f['momentum_x'][:]
    py  = f['momentum_y'][:]
    # pz  = f['momemtum_z'][:]
    E = f['Energy'][:]
    GE = f['GasEnergy'][:]

    f.close()

    vx = px/d
    vy = py/d
    # vz = pz/d
  
    KE = 0.5 * d * (vx*vx + vy*vy)
    # GE = E - KE

    mu = 0.6 # mean molecular weight (mu) (we should add this to the header, when relevant)

    d_cgs = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable

    n = d_cgs/(mu*mp) # number density, particles per cm^3


    #P = GE * (gamma - 1.0) * p_c #pressure

    T = GE * (gamma - 1.0) * p_c / (n * kb) 
    logT = np.log10(T)

    print(np.min(logT))
    print(np.max(logT))

    # print('p_c =', p_c)
    # print('d_c =', d_c)
    # print('n min/max =', np.min(n), np.max(n))
    # print('GE_cgs min/max =', np.min(GE*p_c), np.max(GE*p_c))

    fig, ax = plt.subplots(figsize=(3,6))
    ax.set_yticks(np.linspace(0,nx,13))
    ax.set_xticks(np.linspace(0,ny,9))
    plt.setp(ax.spines.values(), color='white')
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='white')
    ax.tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1)
    # ax.set_title('d1 = 3.0  d2 = 1.0')
    # ax.text(0.03, 0.95, '512 res', transform=ax.transAxes, color='white')
    image = ax.imshow(np.rot90(logT[:, :ny//2].T), cmap='magma', vmin=18.3, vmax=21.5) #originally -31.4 and -30.8
    plt.hlines(y=11*nx//12, xmin=ny//8, xmax=2*ny//8, linewidth=1, color='white')
    plt.text(2.3*ny//8, 11.1*nx//12, '1 pc', fontsize=8, color='white')

    # add a colorbar
    # divider = make_axes_locatable(ax)
    # cbax = divider.append_axes('right', size='5%', pad=0.05)
    # cb = plt.colorbar(image, cax = cbax)
    # cbax.tick_params(axis='y', direction='in', labelcolor='white')
    # cb.solids.set_edgecolor('face')
    # cbax.set_ylabel(r'$\mathrm{log}_{10}(\rho_A)$ [$\mathrm{g}\mathrm{cm}^{-2}$]')
    # plt.show()

    # save the figure
    plt.savefig(dnameout +str(i)+ '.png', dpi=300, bbox_inches='tight', pad_inches = 0, facecolor='black')
    plt.close()
        