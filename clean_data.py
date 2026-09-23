import os
from csv import writer
import h5py
import numpy as np

datadir = '../../../../../projects/jcombar1/betty-testing/HJL/data/cw/test/hdf5_super_ct/'
outdir = '../../../../../projects/jcombar1/betty-testing/HJL/data/cw/test/hdf5_super_ct/'
CAT = 1

box_length   = 2.4
cloud_thresh = 3
n_init       = 1.0

MASS     = 1
VELOCITY = 0

ns    = 0
ne    = 300
nstep = 10

# Constants
mp = 1.672622e-24  # mass of hydrogen atom (g)
kb = 1.380658e-16  # Boltzmann constant (erg/K)
mu = 0.6           # mean molecular weight
km = 1e-5          # cm/s to km/s conversion

# --- Output files ---
mass_file = os.path.join(outdir, "cloud_masses.txt")
vel_file  = os.path.join(outdir, "cloud_velocities.txt")

# --- Count already written lines ---
def count_lines(filepath):
    if not os.path.exists(filepath):
        return 0
    with open(filepath, "r") as f:
        return sum(1 for line in f if line.strip())

mass_done = count_lines(mass_file) if MASS     else 0
vel_done  = count_lines(vel_file)  if VELOCITY else 0

if mass_done >= ne:
    MASS = 0
    print("Mass file exists, creating velocity file.")
if vel_done >= ne: 
    VELOCITY = 0
    print("Velocity file exists, creating mass file.")

if MASS and VELOCITY:
    lines_done = min(mass_done, vel_done)
elif MASS:
    lines_done = mass_done
elif VELOCITY:
    lines_done = vel_done
else:
    print("Both files already complete, exiting.")
    exit()

if lines_done > 0:
    print(f"Resuming from line {lines_done}.")

# --- Get initial cloud mass from snapshot 0 ---
if CAT:
    f = h5py.File(datadir + str(0) + '/' + str(0) + '.h5', 'r')
else:
    f = h5py.File(datadir + 'raw/' + str(0) + '/' + str(0) + '.h5.0', 'r')

head    = f.attrs
nx      = head['dims'][0]
d_c     = head['density_unit']
d       = f['density'][:]
f.close()

dx        = box_length / nx
n         = d * d_c / (mu * mp)
mass      = d * dx**3
mass_init = np.sum(mass[n > n_init / cloud_thresh])

# --- Main loop ---
f_mass = open(mass_file, "a") if MASS     else None
f_vel  = open(vel_file,  "a") if VELOCITY else None

for idx, i in enumerate(range(ns, ne, nstep)):
    if idx < lines_done:
        print(f"Skipping step {i} (already written).")
        continue

    print(str(i))

    if CAT:
        f = h5py.File(datadir + str(i) + '/' + str(i) + '.h5', 'r')
    else:
        f = h5py.File(datadir + 'raw/' + str(i) + '/' + str(i) + '.h5.0', 'r')

    head = f.attrs
    t    = head['t'][0]
    nx   = head['dims'][0]
    d_c  = head['density_unit']
    v_c  = head['velocity_unit']

    d  = f['density'][:]
    px = f['momentum_x'][:]
    f.close()

    dx   = box_length / nx
    n    = d * d_c / (mu * mp)
    mass = d * dx**3
    vx   = (px * v_c * km) / d

    cloud_mask = n > n_init / cloud_thresh
    cloud_mass = mass[cloud_mask]
    mass_tot   = np.sum(cloud_mass)
    mass_cur   = mass_tot / mass_init
    v_avg      = np.sum(vx[cloud_mask] * cloud_mass) / mass_tot
    escaped    = int(np.any(n[-1, :, :] > n_init / cloud_thresh))

    if f_mass:
        f_mass.write(f"{escaped} {t} {mass_cur}\n")
    if f_vel:
        f_vel.write(f"{escaped} {t} {v_avg}\n")

if f_mass:
    f_mass.close()
if f_vel:
    f_vel.close()