import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 1 # Dual Energy Flag
CAT = 1 # Conactenated Flag

adiabatic_dir='../../data/cloud_wind/1.2/16/hdf5/250_proj.h5' # directory where the file is located
cooling_dir='../../data/cloud_wind/3/16/hdf5/250_proj.h5' # directory where the file is located
dnameout='../../plots/compare100/' # directory where the plot will be saved

# t_cc = 4.89e4 # (vwind = 10 km/s)
t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
# iend = 102

#####################################################################
f = h5py.File(adiabatic_dir, 'r') 
head = f.attrs # read the header attributes into a structure, called head

gamma = head['gamma'] # ratio of specific heats
t_adiabatic  = head['t'] # time of this snapshot, in kyr
nx = head['dims'][0] # number of cells in the x direction
ny = head['dims'][1] # number of cells in the y direction
nz = head['dims'][2] # number of cells in the z direction
dx = head['dx'][0] # width of cell in x direction
l_c = head['length_unit']
t_c = head['time_unit']
m_c = head['mass_unit']
d_c = head['density_unit']          
d  = f['d_xy'][:]

f.close()

d = d * m_c / (l_c**2)
n = d / (mu*mp) # number density, particles per cm^3  
logn_adiabatic = np.log10(n)[:,int(ny/6):int(5*ny/6)]

#############################################################
f = h5py.File(cooling_dir, 'r') 
head = f.attrs # read the header attributes into a structure, called head

gamma = head['gamma'] # ratio of specific heats
t_cooling  = head['t'] # time of this snapshot, in kyr
nx = head['dims'][0] # number of cells in the x direction
ny = head['dims'][1] # number of cells in the y direction
nz = head['dims'][2] # number of cells in the z direction
dx = head['dx'][0] # width of cell in x direction
l_c = head['length_unit']
t_c = head['time_unit']
m_c = head['mass_unit']
d_c = head['density_unit']          
d  = f['d_xy'][:]

f.close()

d = d * m_c / (l_c**2)
n = d / (mu*mp) # number density, particles per cm^3  
logn_cooling = np.log10(n)

###############################################################################

units = '$log_{10}(N_{H} [cm^{-2}])$'

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7,4))
fig_color = 'black'
bg_color = 'white'

im = axs[0].imshow(logn_adiabatic.T, cmap='viridis', vmin=19.6, vmax=20) 
axs[0].text(0.03*nx,0.85*ny,'adiabatic', size=10, color='white')
axs[0].text(0.85*nx,0.85*ny, str(int(t_adiabatic/1000)) + r' $Myr$', size=10, color='white')
axs[0].set_xticks(np.linspace(0,nx,9))
axs[0].set_yticks(np.linspace(0,nz,5))
axs[0].invert_yaxis()
axs[0].tick_params(axis='both', which='both', direction='in', color='white', bottom=1, left=1, top=1, right=1, 
        labelleft=0, labelbottom=0, labeltop=0, labelright=0)


axs[1].imshow(logn_cooling.T, cmap='viridis', vmin=19.6, vmax=20.6) 
axs[1].text(0.03*nx,0.85*ny,'cooling', size=10, color='white')
axs[1].text(0.85*nx,0.85*ny, str(int(t_cooling/1000)) + r' $Myr$', size=10, color='white')
axs[1].set_xticks(np.linspace(0,nx,9))
axs[1].set_yticks(np.linspace(0,nz,5))
axs[1].invert_yaxis()
axs[1].tick_params(axis='both', which='both', direction='in', color='white', bottom=1, left=1, top=1, right=1, 
        labelleft=0, labelbottom=1, labeltop=0, labelright=0, labelcolor=fig_color, labelsize=8)
axs[1].set_xticklabels(np.round(np.arange(0,nx*dx+.01,0.3),1))
[l.set_visible(False) for (i,l) in enumerate(axs[1].xaxis.get_ticklabels()) if i % 2 != 0]
axs[1].set_xlabel('$kpc$', size=10, color=fig_color)

for j in range(2):
    plt.setp(axs[j].spines.values(), color='white')
    plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color='white')


fig.subplots_adjust(wspace=0.3, hspace=0.1)

cb = fig.colorbar(im, ax=axs.ravel().tolist(), aspect=40, pad=.03)
cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
cb.ax.yaxis.set_tick_params(color=fig_color, labelsize=8)
cb.outline.set_edgecolor(fig_color)
plt.setp(cbar_yticks, color=fig_color)
cb.ax.set_ylabel(units, size=10, color=fig_color) 

# fig.subplots_adjust(wspace=0.3, hspace=0.1)

# fig.text(0.5, 0.9, str(int(t))+r' $t_{cc}$', size=10, color=fig_color)
# plt.suptitle('Density Projections', color=fig_color,fontsize = 12, y=.93)
plt.savefig(dnameout + 'plot.png', dpi=300, 
        bbox_inches='tight', pad_inches = 0.2, transparent=True) #facecolor=bg_color
plt.close(fig)