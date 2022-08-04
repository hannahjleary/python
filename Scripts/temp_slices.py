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

f.close()

mu = 1.0

d = d*d_c
n = d / (mp * mu)

T = GE * (gamma - 1.0) * p_c / (n * kb) 
T = np.log10(T)

for i in range(ny):
    Tslice_xz = T[:,int(i),:]

    fig = plt.figure(figsize=(16,4))
    image = plt.imshow(Tslice_xz.T, origin ='lower', cmap='plasma')
    cb = plt.colorbar(image, label='$log_{10}(K)$',ticks=np.arange(2.5, 7.5, 0.5))
    
    plt.savefig('../Plots/t_xz_' + str(i) + '.png', dpi=300)
    plt.close(fig)