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
ADIABATIC=0
RADIATIVE=1

dnamein='../../../../../ix/eschneider/hjl28/data/radiative/super/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/plots/radiative/super/png/' # directory where the plot will be saved

res = ['4/', '8/', '16/', '32/', '48/']
res_labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$', '$R_{32}$', '$R_{48}$']
cat = [False, False, True, True, True]

###### Cloud Crushing Times ###################
# t_cc = 4.89e4 # (vwind = 10 km/s)
# t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)

Rcl = .05 # in kpc

# Density, Temperature, Velocity
titles = ['Column Density', 'Temperature', '$\hat{x}$ Velocity']
cols = ['50', '150', '300'] #195, 
labels = ['$log_{10}(N_{H} \ [cm^{-2}])$', '$log_{10}(T \ [K])$', '$kms^{-1}$']

if ADIABATIC:
    figsize = (4.1,4.3)
    col1_slice = 2
    vlims = [19.3, 20.0]
    cb_left = 0.913
    cb_width = 0.019
    model = 'a1000'
    #add + '/' + cols[i] 
if RADIATIVE:
    figsize=(5.9,4.4)
    col1_slice = 3
    vlims = [19.2, 20.6]
    cb_left = 0.91
    cb_width = 0.017

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

                                                            #r: 5.9, 4.4 
fig, axs = plt.subplots(nrows=len(res), ncols=3, figsize=figsize, gridspec_kw={'width_ratios':[1, col1_slice, col1_slice], 'wspace':0.0, 'hspace':0.1}) #4.8
# print(str(axs.shape))

for i in range(len(cols)):

    for j in range(len(res)):
        if cat[j]:
            f = (h5py.File(dnamein + res[j] + 'hdf5/' + cols[i] + '_proj.h5', 'r'))
        else:
            f = h5py.File(dnamein + res[j] + 'hdf5/' + cols[i] + '_proj.h5.0', 'r') 
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

        if i == 0:
            im = axs[j][i].imshow(np.log10(n)[0:int(nx/col1_slice),:].T, cmap=sns.color_palette("rocket", as_cmap=True), vmin = vlims[0], vmax = vlims[1]) # 19.27 20.5
            axs[j][i].set_ylabel(res_labels[j], size=fsize, labelpad=8, rotation='horizontal', ha='right', va='center', color=text_color)
            axs[j][i].set_xticks(np.linspace(0,nz,5))
            axs[j][i].set_yticks(np.linspace(0,nz,5))
            if j == len(res)-1:
                axs[j][i].tick_params(axis='both', which='both', direction='in', color=bg_color, bottom=1, left=1, top=1, right=1, 
                            labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.6, length=3)
                # axs[j][i].set_xlabel('$R_{cl}$', size=fsize, color=fig_color)
                [l.set_visible(False) for (i,l) in enumerate(axs[j][i].xaxis.get_ticklabels()) if i % 2 != 0]
                # axs[j][i].set_xticklabels(dimensions.astype(int))
            else:
                axs[j][i].tick_params(axis='both', which='both', direction='in', color=bg_color, bottom=1, left=1, top=1, right=1, 
                        labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.6, length=3)
        else:
            im = axs[j][i].imshow(np.log10(n).T, sns.color_palette("rocket", as_cmap=True), vmin = vlims[0], vmax = vlims[1]) 
            # axs[j].set_ylabel(labels[j], size=10, color=fig_color)
            axs[j][i].set_xticks(np.linspace(0,nx,9))
            axs[j][i].set_yticks(np.linspace(0,nz,5))
            axs[j][i].invert_yaxis()
            # if i==2 & j == 0:
            #     axs[0][2].text(0.09*nx, 0.25*ny, model, fontsize=8, color='white')

            if j == (len(res)-1):
                axs[j][i].tick_params(axis='both', which='both', direction='in', color=bg_color, bottom=1, left=1, top=1, right=1, 
                        labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.6, length=3)
                axs[j][i].set_xticklabels(dimensions.astype(int))
                # print(nx*dx)
                [l.set_visible(False) for (i,l) in enumerate(axs[j][i].xaxis.get_ticklabels()) if i % 2 != 0]
                # axs[j][i].set_xlabel('$R_{cl}$', size=fsize, color=fig_color)


                # fig.subplots_adjust(wspace=-0.03, hspace=0.08)
                # if i == 2:
                #     col = []
                #     for x in range(j+1):
                #         col.append(axs[x][i])
                # cb = fig.colorbar(im, ax=axs[:,2], aspect=40, pad=.06)
                # cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
                # cb.ax.yaxis.set_tick_params(color=text_color, labelsize=fsize)
                # cb.outline.set_edgecolor(cb_spine_color)
                # plt.setp(cbar_yticks, color=text_color)
                # def truncate(value, _):
                #     return f"{math.floor(value * 10) / 10:.1f}"
                # cb.ax.yaxis.set_major_formatter(mticker.FuncFormatter(truncate))
                # # [l.set_visible(False) for (i,l) in enumerate(cb.ax.get_yticklabels()) if i % 2 != 0]
                # cb.ax.set_ylabel(labels[0], size=fsize, color=text_color)
                
            else:
                axs[j][i].tick_params(axis='both', which='both', direction='in', color=bg_color, bottom=1, left=1, top=1, right=1, 
                        labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.6, length=3)


        plt.setp(axs[j][i].spines.values(), color=spine_color)
        plt.setp([axs[j][i].get_xticklines(), axs[j][i].get_yticklines()], color=bg_color)

        if j == 0:
            # axs[j][i].set_title(str(int(t/1000))+r' $Myr$', fontsize=8, color=fig_color) 
            axs[j][i].set_title(str(int((t/t_cc, 0)[0]))+r' $t_{cc}$', fontsize=fsize, color=text_color)
        
# divider = make_axes_locatable(axs[:,2])
# fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([cb_left, 0.104, cb_width, 0.794])
cb = fig.colorbar(im, cax=cbar_ax)
cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
cb.ax.yaxis.set_tick_params(color=text_color, labelsize=fsize)
cb.outline.set_edgecolor(cb_spine_color)
plt.setp(cbar_yticks, color=text_color)
cb.ax.set_ylabel(labels[0], size=fsize, color=text_color)
def truncate(value, _):
    return f"{math.floor(value * 10) / 10:.1f}"
cb.ax.yaxis.set_major_formatter(mticker.FuncFormatter(truncate))

# fig.subplots_adjust(wspace=-0.02, hspace=0.05)
# fig.text(0.5, 0.94, str(int(t/t_cc))+r' $t_{cc}$', size=9, color=fig_color)
# fig.text(0.5, 0.9, str(t)+r' $Myr$', size=8, color=fig_color)

plt.savefig(dnameout + 'projections.png', dpi=300, 
        bbox_inches='tight', pad_inches = 0.01, facecolor=bg_color) #facecolor=bg_color
plt.close(fig)
