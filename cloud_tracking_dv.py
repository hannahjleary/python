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

dnamein='../../../../../ix/eschneider/hjl28/data/tests/cloud_tracking/' # directory where the file is located
dnameout='../../../../../ix/eschneider/hjl28/plots/tests/cloud_tracking/' # directory where the plot will be saved

iend = 300
nstep = 10

vwind = 100
t_cc = 4.89e3 # cloud crushing time in kyr (vwind = 100 km/s)
tick_labels = np.arange(0, 60+1, 10)
# linelabel1 = "t = 4 $t_{cc}$"
# linelabel2 = "t = 7 $t_{cc}$"
# textheight1 = 0.1
# textheight2 = 0.1
# ram_scale = 1.0
cloud_thresh = 3
min = -0.07
vmax = 1.19
mmax = 4.6
LEGEND = 1
CUTOFFS = 1
DARKMODE = 1

masses_reg = np.zeros(int(iend/nstep))
velocities_reg = np.zeros(int(iend/nstep))

masses_ct = np.zeros(int(iend/nstep))
velocities_ct = np.zeros(int(iend/nstep))

cutoffs = []
t = []

i = 0
with open(os.path.join(dnamein+  "hdf5_large/dm_" + str(3) + ".csv"), "r") as f:
    for line in f:
        if line == '\n':
            break
        line = line.split(",")
        # print(len(line))
        # print(j, i, line[2])
        if not cutoffs:
            if int(line[0]) == 1:
                cutoffs.append(int(i/nstep))
        t.append(float(line[1]))
        masses_reg[i] = float(line[2])
        i += 1

i = 0
with open(os.path.join(dnamein+  "hdf5_large/dv_" + str(3) + ".csv"), "r") as f:
    for line in f:
        if line == '\n':
            break
        line = line.split(",")
        val = line[2].strip()
        if val.lower() == 'nan' or val == '':
            i += 1
            continue
        velocities_reg[i] = float(val)/vwind
        i += 1

i = 0
with open(os.path.join(dnamein+  "hdf5_large_ct/dm_" + str(3) + ".csv"), "r") as f:
    for line in f:
        if line == '\n':
            break
        line = line.split(",")
        # print(len(line))
        # print(j, i, line[2])
        if not cutoffs:
            if int(line[0]) == 1:
                cutoffs.append(int(i/nstep))
        t.append(float(line[1]))
        masses_ct[i] = float(line[2])
        i += 1

i=0
with open(dnamein + "hdf5_large_ct/cloud_velocities.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            if i % 10 == 0 and i//10 < 30:
                velocities_ct[int(i/10)] = float(line)/vwind
            i += 1

# print(masses_reg)
# print(velocities_ct)
# print(t)

if DARKMODE:
    text_color = 'white'
    spine_color = 'white'
    bg_color = 'black'
else:
    text_color = 'black'
    spine_color = 'black'
    bg_color = 'white'

fsize1 = 10
fsize2 = 12

fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(5,11), sharex = True)

ax1.plot(masses_reg, lw=2, label='No cloud tracking')
ax1.plot(masses_ct, lw=2, label='Cloud tracking')

ax2.plot(velocities_reg, lw=2)
ax2.plot(velocities_ct, lw=2)

ax2.set_ylabel("$[v_{c}/v_{w}]$", color=text_color, fontsize=fsize2, labelpad=15)
ax1.tick_params(labelsize=fsize2, labelcolor=text_color, direction='in', length=10, width=2)
ax1.set_ylabel(r"$[M(\rho > \rho_{\rm cl}/3) / M_i]$", color=text_color, fontsize=fsize2, labelpad=15)
leg = ax1.legend(loc='upper left', fontsize=fsize2, facecolor=bg_color, edgecolor=spine_color)
for text in leg.get_texts():
    text.set_color(text_color)
ax2.set_xlabel("t $[t_{cc}]$", color=text_color, fontsize=fsize2, labelpad=10)
ax2.set_xticks(np.linspace(0, 29, 7))  # 7 tick marks from index 0 to 29
ax2.set_xticklabels([0, 1, 2, 3, 4, 5, 6])
ax2.set_xlabel("t [Myr]", color=text_color, fontsize=fsize2, labelpad=10)
ax2.tick_params(labelsize=fsize2, labelcolor=text_color, direction='in', length=10, width=2)

fig.set_facecolor(bg_color)
ax1.set_facecolor(bg_color)
ax2.set_facecolor(bg_color)
plt.setp(ax1.spines.values(), color=spine_color)
plt.setp(ax2.spines.values(), color=spine_color)
plt.setp([ax1.get_xticklines(), ax1.get_yticklines()], color=spine_color)
plt.setp([ax2.get_xticklines(), ax2.get_yticklines()], color=spine_color)

fig.subplots_adjust(wspace=0., hspace=0)
plt.savefig(dnameout + 'ct_dv_dm_' + str(cloud_thresh) + '.png', dpi=300, 
        bbox_inches='tight', pad_inches = 0.3, facecolor=bg_color) #facecolor=bg_color
plt.close(fig)