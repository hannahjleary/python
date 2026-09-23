import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import font_manager
import os
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

plt.style.use('classic')
font_path = os.path.expanduser("~/.fonts/Helvetica.ttf")
font_manager.fontManager.addfont(font_path)
helvetica = font_manager.FontProperties(fname=font_path)
plt.rcParams['font.family'] = helvetica.get_name()
plt.rcParams.update({'font.family': 'Helvetica'})
plt.rcParams['mathtext.default']='regular'

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 0 # Dual Energy Flag
DARKMODE=0
LIGHTMODE=1

dnamein='../../../../../ix/eschneider/hjl28/data/archive/paper2026-PPMP/radiative/sub/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/plots/poster_plots/BUGSposter/' # directory where the plot will be saved

res = ['4/', '8/', '16/', '32/', '48/']
res_labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$', '$R_{32}$', '$R_{48}$']
cat = [False, False, True, True, True]

###### Cloud Crushing Times ###################
# t_cc = 4.89e4 # (vwind = 10 km/s)
t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)

Rcl = .05 # in kpc

istart = 0
iend = 1

Tmin = 0
Tmax = 200

nmin = 19.2 #19.55
nmax = 20.6

vmin = -200 #-200
vmax = 1200 #1200

# Density, Temperature, Velocity
titles = ['Column Density', 'Temperature', '$\hat{x}$ Velocity']
mins = [nmin, Tmin, vmin]
maxs = [nmax, Tmax, vmax]
cmaps = ['viridis', 'plasma', 'magma']
labels = ['$log_{10}(N_{H} \ [cm^{-2}])$', '$log_{10}(T \ [K])$', '$kms^{-1}$']


fig, axs = plt.subplots(nrows=len(res), ncols=1, figsize=(6,9)) #4.8
# print(str(axs.shape))

if LIGHTMODE:
    spine_color='white'
    text_color='black'
    bg_color='white'
    cb_spine_color='black'

if DARKMODE:
    spine_color='white'
    cb_spine_color=spine_color
    text_color = 'white'
    bg_color = 'black'

for i in range(istart, iend):    

    for j in range(len(res)):
        if cat[j]:
            f = h5py.File(dnamein + res[j] + 'hdf5/' + str(i) + '_proj.h5', 'r')
        else:
            f = h5py.File(dnamein + res[j] + 'hdf5/' + str(i) + '_proj.h5.0', 'r') 
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
        T = f['T_xy'][:]

        if DE:
            GE = f['GE_xy'][:]

        f.close()

        dimensions = np.arange(0,nx*dx+.01,(nx*dx+0.01)/8.01)/Rcl

        d = d * m_c / (l_c)**2
        n = d / (mu*mp) # number density, particles per cm^3  
        n_tot = np.sum(n)

        logT = np.log10(T)
        # print(np.min(np.log10((n))), np.max(np.log10(n)))

        fsize = 12
        
        im = axs[j].imshow(np.log10(n).T, sns.color_palette("rocket", as_cmap=True), vmin = nmin, vmax = nmax) 
        # axs[j].set_ylabel(labels[j], size=10, color=fig_color)
        axs[j].set_xticks(np.linspace(0,nx,9))
        axs[j].set_yticks(np.linspace(0,nz,5))
        axs[j].invert_yaxis()
        axs[j].tick_params(axis='both', which='both', direction='in', color=bg_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.6, length=3)

        axs[len(res)-1].tick_params(axis='both', which='both', direction='in', color=bg_color, bottom=1, left=1, top=1, right=1, 
                labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.6, length=3)
        # axs[j][i].set_xticklabels(dimensions.astype(int))
        # print(nx*dx)
        [l.set_visible(False) for (i,l) in enumerate(axs[j].xaxis.get_ticklabels()) if i % 2 != 0]
  
    

        plt.setp(axs[j].spines.values(), color=spine_color)
        plt.setp([axs[j].get_xticklines(), axs[j].get_yticklines()], color=bg_color)

        axs[0].set_title(str(int((t/t_cc, 0)[0]))+r' $t_{cc}$', fontsize=fsize, color=text_color)

    cbar_ax = fig.add_axes([0.91, 0.104, 0.017, 0.794])
    cb = fig.colorbar(im, cax=cbar_ax)
    cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
    cb.ax.yaxis.set_tick_params(color=text_color, labelsize=fsize)
    cb.outline.set_edgecolor(cb_spine_color)
    plt.setp(cbar_yticks, color=text_color)
    cb.ax.set_ylabel(labels[0], size=fsize, color=text_color)
    def truncate(value, _):
        return f"{math.floor(value * 10) / 10:.1f}"
    cb.ax.yaxis.set_major_formatter(mticker.FuncFormatter(truncate))

    fig.text(0.5, 0.94, str(int(t/t_cc))+r' $t_{cc}$', size=9, color=text_color)
    # fig.text(0.5, 0.9, str(t)+r' $Myr$', size=8, color=fig_color)

    plt.savefig(dnameout + 'projections.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.15, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)
