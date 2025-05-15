import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
import seaborn as sns
import h5py
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

# font_path = "../../fonts/Helvetica.ttf"  # Your font path goes here
# # font_manager.fontManager.addfont(font_path)
# prop = font_manager.FontProperties(fname=font_path)
# plt.rcParams['font.family'] = prop.get_name()

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

DE = 0 # Dual Energy Flag
DARKMODE = 0
MASS = 1
VELOCITY = 0

dnamein='../../data/cloud_wind/4/' # directory where the file is located
dnameout='../../data/cloud_wind/4/png/' # directory where the plot will be saved

num = 50
nstep = 10

masses = [np.full(num, np.nan), np.full(num, np.nan), np.full(num, np.nan), np.full(num, np.nan), np.full(num, np.nan)]
velocities = [np.full(num, np.nan), np.full(num, np.nan), np.full(num, np.nan), np.full(num, np.nan), np.full(num, np.nan)]
t = []
cutoffs = [0, 0, 0, 0, 0]

cloud_thresh = 3
# t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
# vwind = 100
t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
vwind = 1000

if MASS:
    j = 0
    i = 0
    with open(os.path.join(dnamein, "dm_" + str(cloud_thresh) + ".csv"), "r") as f:
        for line in f:
            if line == '\n':
                j+= 1
                i = 0
                continue
            line = line.split(",")
            # print(j, i, line[2])
            if cutoffs[j] == 0:
                if int(line[0]) == 1:
                    cutoffs[j] = int(i/nstep)
            # t.append(float(line[1]))
            masses[j][int(i/nstep)] = line[2]
            i += 10

if VELOCITY:
    j = 0
    i = 0
    k = 0
    with open(os.path.join(dnamein, "dv_" + str(cloud_thresh) + ".csv"), "r") as f:
        for line in f:
            if line == '\n':
                j+= 1
                i = 0
                continue
            line = line.split(",")
            print(j, i, line[2])
            if cutoffs[j] == 0:
                if int(line[0]) == 1:
                    cutoffs[j] = int(i/nstep)
            # t.append(float(line[1]))
            velocities[j][int(i/nstep)] = float(line[2])/vwind
            i += 10
            k += 1

################### Plotting #################################

if DARKMODE:
    fig_color = 'white'
    bg_color = 'grey'
else:
    fig_color = 'black'
    bg_color = 'white'

colors = sns.color_palette('rocket_r', len(velocities))
# plt.rcParams.update({"font.family" : "Helvetica"})

labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$', '$R_{32}$', '$R_{48}$']
tick_labels = np.arange(0, 20+1, 4) # 0, upper limit + 1, upper limit / 5
fsize = 22


if MASS and VELOCITY:
    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(8,15), sharex = True)
    ax1.axvline(x=20, color='gray', linestyle='dashed')
    ax2.text(21, .03, 't = 8.0 $t_{cc}$', color='gray', rotation=90, fontsize=18)
    ax1.axvline(x=35, color='gray', linestyle='dashed')
    ax2.text(36, .03, 't = 14.0 $t_{cc}$', color='gray', rotation=90, fontsize=18)
    ax2.axvline(x=20, color='gray', linestyle='dashed')
    ax2.axvline(x=35, color='gray', linestyle='dashed')
    for s in range(len(velocities)):
        ax1.plot(masses[s], label=labels[s], color=colors[s], linewidth=2.5)
    for s in range(len(masses)):
        ax2.plot(velocities[s], label=labels[s], color=colors[s], linewidth=2.5)   

    ax2.set_ylabel("$[v_{c}/v_{w}]$", color=fig_color, fontsize=fsize, labelpad=15)
    ax1.tick_params(labelsize=fsize, direction='in', length=10, width=2)
    ax1.set_ylabel(r"$[M(\rho > \rho_{\rm cl}/3) / M_i]$", color=fig_color, fontsize=fsize, labelpad=15)
    ax2.legend(loc='upper left', fontsize=20)
    ax2.set_xlabel("t $[t_{cc}]$", color=fig_color, fontsize=fsize, labelpad=10)
    ax2.set_xticks(np.arange(0, 50+1, nstep))
    ax2.set_xticklabels(tick_labels)
    ax2.tick_params(labelsize=fsize, direction='in', length=10, width=2)

    fig.set_facecolor(bg_color)
    plt.setp(ax1.spines.values(), color=fig_color)
    plt.setp(ax2.spines.values(), color=fig_color)
    plt.setp([ax1.get_xticklines(), ax1.get_yticklines()], color=fig_color)
    plt.setp([ax2.get_xticklines(), ax2.get_yticklines()], color=fig_color)

    fig.subplots_adjust(wspace=0., hspace=0)
    plt.savefig(dnameout + 'summary_' + str(cloud_thresh) + '.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.3, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)


if MASS and not VELOCITY:
    for s in range(len(masses)):
        fig, ax = plt.subplots(figsize=(8,7))
        ax.plot(masses[s], label=masses[s])
        # if s == len(velocities)-1:
        #     ax.plot(cutoffs[s], velocities[s][int(cutoffs[s])], c='black', linestyle=' ', marker='X', label='Mass Leaving Box')
        # else:
        #     ax.plot(cutoffs[s], velocities[s][int(cutoffs[s])], c='black', linestyle=' ', marker='X')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlabel("t $[Myr]$", color=fig_color, fontsize=12)
    ax.set_ylabel("$[M/M_{i}]$", color=fig_color, fontsize=12)
    ax.set_title('Change in Cloud Mass', fontsize=12)
    ax.set_xticks(np.arange(0, 50+1, nstep))
    ax.set_xticklabels(tick_labels)
    ax.tick_params(labelsize=9)
    # ax.set_xticklabels(n_step*np.arange(0, num))
        # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]

    ax.set_facecolor(bg_color)
    plt.setp(ax.spines.values(), color=fig_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

    plt.savefig(dnameout + 'vm.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)

if VELOCITY and not MASS:
    fig, ax = plt.subplots(figsize=(8,7))
    for s in range(len(velocities)):
        ax.plot(velocities[s], label=labels[s], color=colors[s])
        if s == len(velocities)-1:
            ax.scatter(cutoffs[s], velocities[s][int(cutoffs[s])], c='black', marker='X', label='Mass Leaving Box')
        else:
            ax.scatter(cutoffs[s], velocities[s][int(cutoffs[s])], c='black', marker='X')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlabel("t $[Myr]$", color=fig_color, fontsize=12)
    ax.set_ylabel("v $[kms^{-1}]$", color=fig_color, fontsize=12)
    ax.set_title('Mean Cloud Velocity', fontsize=12)
    ax.set_xticks(np.arange(0, 50+1, nstep))
    ax.set_xticklabels(tick_labels)
    ax.tick_params(labelsize=9)
    # ax.set_xticklabels(n_step*np.arange(0, num))
        # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]

    ax.set_facecolor(bg_color)
    plt.setp(ax.spines.values(), color=fig_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=fig_color)

    plt.savefig(dnameout + 'vx_dcloud_new.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)
