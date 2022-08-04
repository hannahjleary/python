import numpy as np
import matplotlib.pyplot as plt
import h5py 
from mpl_toolkits.axes_grid1 import make_axes_locatable

#dnamein = '../Data/'
dnameout = '../Plots/'

mp = 1.672622e-24   # mass of hydrogren atom (grams)
kb = 1.38065e-16    # boltzmann constant (ergs/K)

f = h5py.File('200_512_f32.h5', 'r')
head = f.attrs

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

d = f['density'][:]
GE = f['GasEnergy'][:]
px = f['momentum_x'][:]
py = f['momentum_y'][:]
pz = f['momentum_z'][:]

f.close()

mu = 1.0

d = d*d_c
n = d / (mp * mu)

P = GE * (gamma - 1.0) * p_c #pressure

T = GE * (gamma - 1.0) * p_c / (n * kb) 
T = np.log10(T)

km = 1e-5

Vx = (px*v_c*km)/d #velocity in the x direction
Vy = (py*v_c*km)/d #velocity in the y direction
Vz = (pz*v_c*km)/d #velocity in the z direction

for i in range(ny):

    nslice_xz = np.log10(n[:,int(ny/2),:])
    Tslice_xz = np.log10(T[:,int(i),:])
    Pslice_xz = P[:,i,:]
    Vxslice_xz = Vx[:,i,:]
    Vyslice_xz = Vy[:,i,:]
    Vzslice_xz = Vz[:,i,:]

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(64,16))

    im = ax1.imshow(Vxslice_xz.T, label='X Velocity', cmap='coolwarm')
    plt.colorbar(im, ax=ax1)

    im = ax2.imshow(nslice_xz.T, label='Number Density')
    plt.colorbar(im, ax=ax2)

    im = ax3.imshow(Vyslice_xz.T, label='Y Velocity', cmap='coolwarm')
    plt.colorbar(im, ax=ax3)

    im = ax4.imshow(Tslice_xz.T, label='Temperature', cmap='plasma')
    plt.colorbar(im, ax=ax4)

    im = ax5.imshow(Vzslice_xz.T, label='Z Velocity', cmap='coolwarm')
    plt.colorbar(im, ax=ax5)

    im = ax6.imshow(Pslice_xz.T, label='Pressure')
    plt.colorbar(im, ax=ax6)

    plt.tight_layout()
    plt.savefig(dnameout + 'all_slices_' + str(i) + '.png', bbox_inches='tight', dpi=300)
    plt.close(fig)