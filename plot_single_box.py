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

# plt.style.use('classic')
# font_path = os.path.expanduser("~/.fonts/Helvetica.ttf")
# font_manager.fontManager.addfont(font_path)
# helvetica = font_manager.FontProperties(fname=font_path)
# plt.rcParams['font.family'] = helvetica.get_name()
# plt.rcParams.update({'font.family': 'Helvetica'})
# plt.rcParams['mathtext.default']='regular'

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 0 # Dual Energy Flag
DARKMODE=1
LIGHTMODE=0

dnamein='../../../../../ix/eschneider/hjl28/data/radiative/super/48/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/plots/radiative/super/boxes/' # directory where the plot will be saved

CAT = 1

###### Cloud Crushing Times ###################
# t_cc = 4.89e4 # (vwind = 10 km/s)
t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)

Rcl = .05 # in kpc

istart = 0
iend = 500

vlims_T = [3.2,6.8]
vlims_n = [19.35,20.6]
vlims_v = [-200,1200]

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

    if CAT:
        f = h5py.File(dnamein + 'hdf5/' + str(i) + '_slice.h5', 'r')
    else:
        f = h5py.File(dnamein + 'hdf5/' + str(i) + '_slice.h5.0', 'r') 
        # f = h5py.File(dnamein + 'hdf5/' + str(i) + '/' + str(i) + '_proj.h5.0', 'r') 
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
    # T = f['T_xy'][:]
    px  = f['mx_xy'][:]
    py  = f['my_xy'][:]
    pz  = f['mz_xy'][:]
    E = f['E_xy'][:]

    if DE:
        GE = f['GE_xy'][:]

    if not DE:   
        vx = px/d
        vy = py/d
        vz = pz/d 
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE

    f.close()

    dimensions = np.arange(0,nx*dx+.01,(nx*dx+0.01)/8.01)/Rcl

    # d = d * m_c / (l_c)**2
    n = d * d_c/ (mu*mp) # number density, particles per cm^3  
    n_tot = np.sum(n)

    T = GE*(gamma-1.0)*p_c / (n*kb) #temperature
    logT = np.log10(T)

    print(np.min(logT),np.max(logT))

    fig, ax = plt.subplots(figsize=(3,6)) #4.8
    fsize = 12

    #  cmap=sns.color_palette("magma", as_cmap=True)
    im = ax.imshow(np.rot90(logT.T), cmap='magma', vmin = vlims_T[0], vmax = vlims_T[1]) # 19.27 20.5
    # axs[j][i].set_ylabel(res_labels[j], size=fsize, labelpad=8, rotation='horizontal', ha='right', va='center', color=text_color)
    ax.set_yticks(np.linspace(0,nx,25))
    ax.set_xticks(np.linspace(0,nz,9))
    ax.tick_params(axis='both', which='both', direction='in', color=spine_color, bottom=1, left=1, top=1, right=1, 
                    labelleft=0, labelbottom=0, labeltop=0, labelright=0, width=0.8, length=4) 
        # axs[j][i].set_xlabel('$R_{cl}$', size=fsize, color=fig_color)
    # [l.set_visible(False) for (i,l) in enumerate(axs[j][i].xaxis.get_ticklabels()) if i % 2 != 0]
        # axs[j][i].set_xticklabels(dimensions.astype(int))
    # plt.arrow(2*nx//16, ny//8, nx//16, 0, linewidth=1, head_width=5, color='white')
    plt.hlines(y=23*nx//24, xmin=ny//8, xmax=2*ny//8, linewidth=1, color='white')
    plt.text(2.5*ny//8, 23.2*nx//24, '100 pc', fontsize=8, color='white')
    plt.text(5.2*ny//8, 1.2*nx//24, str(int(t/1000 + 20))+' Myr', fontsize=8, color='white')
    # plt.text(nx//16, 1.2*ny//8, '$v_w$', fontsize=10, color='white')
    plt.setp(ax.spines.values(), color=spine_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=spine_color)

    # ax.set_title(str(int((t/t_cc, 0)[0]))+r' $t_{cc}$', fontsize=fsize, color=text_color)

    # cbar_ax = fig.add_axes([0.91, 0.104, 0.017, 0.794])
    # cb = fig.colorbar(im, ax=ax, aspect=40, pad=.06)
    # cbar_yticks = plt.getp(cb.ax.axes, 'yticklabels')
    # cb.ax.yaxis.set_tick_params(color=text_color, labelsize=fsize)
    # cb.outline.set_edgecolor(cb_spine_color)
    # plt.setp(cbar_yticks, color=text_color)
    # # cb.ax.set_ylabel('label', size=fsize, color=text_color)
    # def truncate(value, _):
    #     return f"{math.floor(value * 10) / 10:.1f}"
    # cb.ax.yaxis.set_major_formatter(mticker.FuncFormatter(truncate))

    # fig.text(0.5, 0.94, str(int(t/t_cc))+r' $t_{cc}$', size=9, color=text_color)
    # fig.text(0.5, 0.9, str(t)+r' $Myr$', size=8, color=fig_color)

    plt.savefig(dnameout + str(i) + '.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)