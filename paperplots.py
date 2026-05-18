import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
import seaborn as sns
import h5py
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

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

DE = 1 # Dual Energy Flag
DARKMODE = 1
MASS = 1
VELOCITY = 1

# model = input("Enter a100, a1000, r100, or r1000: ")
model = sys.argv[1]

if model == "a100":
    dnamein='../../../../../ix/eschneider/hjl28/data/adiabatic_PPMP/sub/plotting_data/' # directory where the file is located
    dnameout='../../../../../ix/eschneider/hjl28/plots/adiabatic_PPMP/sub/png/' # directory where the plot will be saved
    vwind = 100
    t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
    tick_labels = np.arange(0, 10+1, 2)
    linelabel1 = "t = 4 $t_{cc}$"
    linelabel2 = "t = 7 $t_{cc}$"
    textheight1 = 0.01 #Mass: 0.5, Vel: 0.01
    textheight2 = 0.01
    ram_scale = 1.0
    min = -0.02 #-0.04 for mass, -0.02 for vel
    vmax = 0.46
    mmax = 1.03
    LEGEND = 1
    CUTOFFS = 0
if model == "a1000":
    dnamein='../../../../../ix/eschneider/hjl28/data/adiabatic_PPMP/super/plotting_data/' # directory where the file is located
    dnameout='../../../../../ix/eschneider/hjl28/plots/adiabatic_PPMP/super/png/' # directory where the plot will be saved
    vwind = 1000
    t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
    tick_labels = np.arange(0, 10+1, 2)
    linelabel1 = "t = 4 $t_{cc}$"
    linelabel2 = "t = 7 $t_{cc}$"
    textheight1 = 0.001
    textheight2 = 0.001
    ram_scale = 4.0
    min = -0.04
    vmax = 0.46
    mmax = 1.03
    LEGEND = 0
    CUTOFFS = 0
if model == "r100":
    # dnamein='../../../../../ix/eschneider/hjl28/data/radiative/sub/plotting_data/' # directory where the file is located
    # dnameout='../../../../../ix/eschneider/hjl28/plots/radiative/sub/png/' # directory where the plot will be saved
    dnamein='../../../../../ix/eschneider/hjl28/data/tests/cloud_tracking/hdf5_large/' # directory where the file is located
    dnameout='../../../../../ix/eschneider/hjl28/plots/tests/cloud_tracking/' # directory where the plot will be saved
    vwind = 100
    t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
    tick_labels = np.arange(0, 10+1, 2)
    # linelabel1 = "t = 4 $t_{cc}$"
    # linelabel2 = "t = 7 $t_{cc}$"
    # textheight1 = 0.1
    # textheight2 = 0.1
    # ram_scale = 1.0
    min = -0.07
    vmax = 1.19
    mmax = 4.6
    LEGEND = 1
    CUTOFFS = 1
if model == "r1000":
    dnamein='../../../../../ix/eschneider/hjl28/data/radiative/super/plotting_data/' # directory where the file is located
    dnameout='../../../../../ix/eschneider/hjl28/plots/radiative/super/png/' # directory where the plot will be saved
    vwind = 1000
    t_cc = 4.89e2 # cloud crushing time in kyr (vwind = 1000 km/s)
    tick_labels = np.arange(0, 20+1, 4)
    linelabel1 = "t = 8 $t_{cc}$"
    linelabel2 = "t = 14 $t_{cc}$"
    textheight1 = 0.1
    textheight2 = 0.1
    ram_scale = 4.0
    min = 0.0
    vmax = 0.54
    mmax = 1.19
    LEGEND = 1
    CUTOFFS = 1

num = 50
nstep = 10

# masses3 


if MASS:
    j = 0
    i = 0
    with open(os.path.join(dnamein, "dm_" + str(3) + ".csv"), "r") as f:
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
            t.append(float(line[1]))
            masses3[j][int(i/nstep)] = float(line[2])
            i += 10
    # j = 0
    # i = 0
    # with open(os.path.join(dnamein, "dm_" + str(5) + ".csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         # print(j, i, line[2])
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         t.append(float(line[1]))
    #         masses5[j][int(i/nstep)] = line[2]
    #         i += 10
    # j = 0
    # i = 0
    # with open(os.path.join(dnamein, "dm_" + str(7) + ".csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         # print(j, i, line[2])
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         t.append(float(line[1]))
    #         masses7[j][int(i/nstep)] = line[2]
    #         i += 10
    # j = 0
    # i = 0
    # with open(os.path.join(dnamein, "dm_" + str(10) + ".csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         # print(j, i, line[2])
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         t.append(float(line[1]))
    #         masses10[j][int(i/nstep)] = line[2]
    #         i += 10
    # j = 0
    # i = 0
    # with open(os.path.join(dnamein, "dm_3wind.csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         # print(j, i, line[2])
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         t.append(float(line[1]))
    #         masseswind[j][int(i/nstep)] = line[2]
    #         i += 10

# if MASS:
#     j = 0
#     i = 0
#     for thresh in threshs:
#         with open(os.path.join(dnamein, "dm_" + str(thresh) + ".csv"), "r") as f:
#             for line in f:
#                 if line == '\n':
#                     j+= 1
#                     i = 0
#                     continue
#                 line = line.split(",")
#                 # print(j, i, line[2])
#                 if cutoffs[j] == 0:
#                     if int(line[0]) == 1:
#                         cutoffs[j] = int(i/nstep)
#                 t.append(float(line[1]))
#                 choice[threshs.index(thresh)][j][int(i/nstep)] = line[2]
#                 i += 10

if VELOCITY:
    j = 0
    i = 0
    k = 0
    with open(os.path.join(dnamein, "dv_" + str(3) + ".csv"), "r") as f:
        for line in f:
            if line == '\n':
                j+= 1
                i = 0
                continue
            line = line.split(",")
            if cutoffs[j] == 0:
                if int(line[0]) == 1:
                    cutoffs[j] = int(i/nstep)
            if j==0:
                t.append(float(line[1]))
            # print(line[0], line[1], line[2])
            velocities3[j][int(i/nstep)] = float(line[2])/vwind
            i += 10
            k += 1
    # j = 0
    # i = 0
    # k = 0
    # with open(os.path.join(dnamein, "dv_" + str(5) + ".csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         if j==0:
    #             t.append(float(line[1]))
    #         # print(line[0], line[1], line[2])
    #         velocities5[j][int(i/nstep)] = float(line[2])/vwind
    #         i += 10
    #         k += 1
    # j = 0
    # i = 0
    # k = 0
    # with open(os.path.join(dnamein, "dv_" + str(7) + ".csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         if j==0:
    #             t.append(float(line[1]))
    #         # print(line[0], line[1], line[2])
    #         velocities7[j][int(i/nstep)] = float(line[2])/vwind
    #         i += 10
    #         k += 1
    # j = 0
    # i = 0
    # k = 0
    # with open(os.path.join(dnamein, "dv_" + str(10) + ".csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         if j==0:
    #             t.append(float(line[1]))
    #         # print(line[0], line[1], line[2])
    #         velocities10[j][int(i/nstep)] = float(line[2])/vwind
    #         i += 10
    #         k += 1
    # j = 0
    # i = 0
    # k = 0
    # with open(os.path.join(dnamein, "dv_3wind.csv"), "r") as f:
    #     for line in f:
    #         if line == '\n':
    #             j+= 1
    #             i = 0
    #             continue
    #         line = line.split(",")
    #         if cutoffs[j] == 0:
    #             if int(line[0]) == 1:
    #                 cutoffs[j] = int(i/nstep)
    #         if j==0:
    #             t.append(float(line[1]))
    #         # print(line[0], line[1], line[2])
    #         velocitieswind[j][int(i/nstep)] = float(line[2])/vwind
    #         i += 10
    #         k += 1

################### Plotting #################################

if DARKMODE:
    text_color = 'white'
    spine_color = 'white'
    bg_color = 'black'
else:
    text_color = 'black'
    spine_color = 'black'
    bg_color = 'white'

colors = sns.color_palette('rocket_r', len(velocities3)+2)
# plt.rcParams.update({"font.family" : "Helvetica"})

labels = ['$R_{4}$', '$R_{8}$', '$R_{16}$', '$R_{32}$', '$R_{48}$']
# tick_labels = np.arange(0, 10+1, 2) # 0, upper limit + 1, upper limit / 5
#0, 10+1, 2 for all sims except radiatice supersonic is 0, 20+1, 4
fsize1 = 10
fsize2 = 12

time = np.array(t)
ram_v = (float(vwind)/(ram_scale * float(10**2) * 3.086e16 * 0.05)) * (3.154e10 * time)
# print((float(vwind)* 3.154e10 * time[5])/(float(10**2) * 3.086e16 * 0.05))

if MASS and VELOCITY:
    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(8,15), sharex = True)
    # ax1.axvline(x=20, color='gray', linestyle='dashed')
    # ax2.text(21, .03, 't = 8.0 $t_{cc}$', color='gray', rotation=90, fontsize=18)
    # ax1.axvline(x=35, color='gray', linestyle='dashed')
    # ax2.text(36, .03, 't = 14.0 $t_{cc}$', color='gray', rotation=90, fontsize=18)
    # ax2.axvline(x=20, color='gray', linestyle='dashed')
    # ax2.axvline(x=35, color='gray', linestyle='dashed')
    for s in range(len(velocities)):
        ax1.plot(masses[s], label=labels[s], color=colors[s], linewidth=2.5)
    for s in range(len(masses)):
        ax2.plot(velocities[s], label=labels[s], color=colors[s], linewidth=2.5)   

    ax2.set_ylabel("$[v_{c}/v_{w}]$", color=fig_color, fontsize=fsize2, labelpad=15)
    ax1.tick_params(labelsize=fsize2, direction='in', length=10, width=2)
    ax1.set_ylabel(r"$[M(\rho > \rho_{\rm cl}/3) / M_i]$", color=fig_color, fontsize=fsize2, labelpad=15)
    ax2.legend(loc='upper left', fontsize=20)
    ax2.set_xlabel("t $[t_{cc}]$", color=fig_color, fontsize=fsize2, labelpad=10)
    ax2.set_xticks(np.arange(0, 50+1, nstep))
    ax2.set_xticklabels(tick_labels)
    ax2.tick_params(labelsize=fsize2, direction='in', length=10, width=2)

    fig.set_facecolor(bg_color)
    plt.setp(ax1.spines.values(), color=fig_color)
    plt.setp(ax2.spines.values(), color=fig_color)
    plt.setp([ax1.get_xticklines(), ax1.get_yticklines()], color=fig_color)
    plt.setp([ax2.get_xticklines(), ax2.get_yticklines()], color=fig_color)

    fig.subplots_adjust(wspace=0., hspace=0)
    plt.savefig(dnameout + 'both_' + model + str(cloud_thresh) + '.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.3, facecolor=bg_color) #facecolor=bg_color
    plt.close(fig)


if MASS and not VELOCITY:
    # sns.set_palette('Set1')
    fig, ax = plt.subplots(figsize=(4,2.3))
    ax.set_facecolor(bg_color)
    # plt.plot(masses3[2], label=r'$\rho > \rho_{cl,i}/3$', linewidth=1.3,zorder=2)
    # plt.plot(masses5[2], label=r'$\rho > \rho_{cl,i}/5$', linewidth=1.3,zorder=2)
    # plt.plot(masses7[2], label=r'$\rho > \rho_{cl,i}/7$', linewidth=1.3,zorder=2)
    # plt.plot(masses10[2], label=r'$\rho > \rho_{cl,i}/10$', linewidth=1.3,zorder=2)
    # plt.plot(masseswind[2], label=r'$\rho > 3\rho_{w}$', linewidth=1.3,zorder=2)
    # plt.axvline(x=5, color='gray', linestyle='dashed', linewidth=1)
    # plt.axvline(x=15, color='gray', linestyle='dashed', linewidth=1)
    # # plt.text(21, textheight1, linelabel1, color='gray', rotation=90, fontsize=fsize1, zorder=1) #17.25
    # plt.axvline(x=30, color='gray', linestyle='dashed', linewidth=1)
    # plt.text(36, textheight2, linelabel2, color='gray', rotation=90, fontsize=fsize1, zorder=1) #32.35
    for s in range(len(masses3)):
        ax.plot(masses3[s], label=labels[s], color=colors[s], linewidth=1.3, zorder=2)
        if CUTOFFS:
            ax.plot(cutoffs[s], masses3[s][int(cutoffs[s])], c=text_color, linestyle=' ', marker='X', markersize=4, zorder=2)
    if LEGEND:
        ax.legend(loc='upper center', ncol=3, fontsize=fsize1, bbox_to_anchor=(0.5,1.36), facecolor=bg_color, labelcolor=text_color, edgecolor=spine_color)
    plt.xlabel("t $[t_{cc}]$", color=text_color, fontsize=fsize2) #$[t_{cc}]$
    # ax.set_xlabel("time", color=text_color, fontsize=fsize2)
    # plt.ylabel(r"$[M(\rho > \rho_{cl}/3)/M_{i}]$", color=text_color, fontsize=fsize2)
    plt.ylabel(r"$[M_{cl}/M_{i}]$", color=text_color, fontsize=fsize2)
    # ax.set_ylabel(r"Cloud Mass", color=text_color, fontsize=fsize2)
    ax.set_ylim(min, mmax)
    plt.text(0.76, 0.79, model, transform=fig.transFigure, fontsize=fsize1)
    ax.set_xticks(np.arange(0, 50+1, nstep))
    ax.set_xticklabels(tick_labels)
    ax.tick_params(labelsize=fsize1, labelbottom=1, color=spine_color, labelcolor=text_color)
    plt.setp(ax.spines.values(), color=spine_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=spine_color)

    plt.savefig(dnameout + 'dm-' + model + '-slide.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color

if VELOCITY and not MASS:
    # sns.set_palette('Set1')
    fig, ax = plt.subplots(figsize=(4,2.3))
    ax.set_facecolor(bg_color)
    # plt.plot(velocities3[2], label=r'$\rho > \rho_{cl,i}/3$', linewidth=1.3,zorder=2)
    # plt.plot(velocities5[2], label=r'$\rho > \rho_{cl,i}/5$', linewidth=1.3,zorder=2)
    # plt.plot(velocities7[2], label=r'$\rho > \rho_{cl,i}/7$', linewidth=1.3,zorder=2)
    # plt.plot(velocities10[2], label=r'$\rho > \rho_{cl,i}/10$', linewidth=1.3,zorder=2)
    # plt.plot(velocitieswind[2], label=r'$\rho > 3\rho_{w}$', linewidth=1.3,zorder=2)
    # plt.axvline(x=5, color='gray', linestyle='dashed', linewidth=1)
    # plt.axvline(x=15, color='gray', linestyle='dashed', linewidth=1)
    # # plt.text(21, textheight1, linelabel1, color='gray', rotation=90, fontsize=fsize1, zorder=1) #21, 0.5
    # plt.axvline(x=30, color='gray', linestyle='dashed', linewidth=1)
    # plt.text(36, textheight2, linelabel2, color='gray', rotation=90, fontsize=fsize1, zorder=1) #36, 0.5
    for s in range(len(velocities3)):
        ax.plot(velocities3[s], label=labels[s], color=colors[s], linewidth=1.3, zorder=2)
        if CUTOFFS:
            ax.plot(cutoffs[s], velocities3[s][int(cutoffs[s])], c=text_color, linestyle=' ', marker='X', markersize=4, zorder=2)
    ax.plot(ram_v, color='darkgray', linewidth=1, label='$v_{ram}$', zorder=1)
    if LEGEND:
        ax.legend(loc='upper center', ncol=3, fontsize=fsize1, bbox_to_anchor=(0.5,1.36), facecolor=bg_color, labelcolor=text_color, edgecolor=spine_color) #columnspacing=0.7
    ax.set_xlabel("t $[t_{cc}]$", color=text_color, fontsize=fsize2)
    # ax.set_xlabel("time", color=fig_color, fontsize=fsize2)
    # ax.set_ylabel(r"$[\bar{v}_x(\rho > \rho_{cl}/3)/v_{w}]$", color=text_color, fontsize=fsize2)
    ax.set_ylabel(r"$[\bar{v}_x/v_{w}]$", color=text_color, fontsize=fsize2)
    ax.set_ylim(min, vmax)
    plt.text(0.23, 0.79, model, transform=fig.transFigure, fontsize=fsize1)
    ax.set_xticks(np.arange(0, 50+1, nstep))
    ax.set_xticklabels(tick_labels)
    ax.tick_params(labelsize=fsize1, labelbottom=1, color=spine_color, labelcolor=text_color)
    plt.setp(ax.spines.values(), color=spine_color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=spine_color)
    # ax.set_xticklabels(n_step*np.arange(0, num))
        # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 10 != 0]
    plt.savefig(dnameout + 'dv-' + model + '-slide.png', dpi=300, 
            bbox_inches='tight', pad_inches = 0.2, facecolor=bg_color) #facecolor=bg_color
