import os
from csv import writer
import h5py
import numpy as np

datadir = '../../data/cloud_wind/4/'
csvdir = '../../data/cloud_wind/4/'

sims = ['48/'] #, '4/', '8/', '16/', '32/', '48/'
cat = [True] # False, False, True, True, True

box_length = 2.4 #2.4 for sims 3 and 4, 1.6 for sims 1 and 2 
cloud_thresh = 3

MASS = 1
VELOCITY = 0

if MASS:
     filename = "dm_" + str(cloud_thresh) + "copy.csv"
if VELOCITY:
     filename = "dv_" + str(cloud_thresh) + "copy.csv"

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6 # mean molecular weight (mu) of 1

f_csv = open(os.path.join(csvdir, filename), "a")
f_csv.close()

ns = 410
ne = 500
nstep = 10

if cat[0]:
    f = h5py.File(datadir + sims[0] + 'hdf5/' +str(0) + '.h5', 'r') 
else:
    f = h5py.File(datadir + sims[0] + 'hdf5/' +str(0) + '.h5.0', 'r') 

head = f.attrs
nx = head['dims'][0] 
d  = f['density'][:]
d_c = head['density_unit'] 
n = d * d_c/(mu*mp) # number density, particles per cm^3 
n_init = 1.0

dx = box_length / nx
mass = d * dx*dx*dx
cloud_mass = mass[n >= (n_init / cloud_thresh)]
mass_tot = np.sum(cloud_mass)
mass_init = mass_tot

for i in range(len(sims)):
    for j in range(ns, ne, nstep):
        print (str(i) + ": " + str(j))
        if cat[i]:
            f = h5py.File(datadir + sims[i] + 'hdf5/' +str(j) + '.h5', 'r') 
        else:
            f = h5py.File(datadir + sims[i] + 'hdf5/' +str(j) + '.h5.0', 'r') 
        head = f.attrs
        t = head['t'][0]
        nx = head['dims'][0] # number of cells in the x direction
        ny = head['dims'][1] # number of cells in the y direction
        nz = head['dims'][2] # number of cells in the z direction
        dx = head['dx'][0] # width of cell in x direction
        dy = head['dx'][1] # width of cell in y direction
        dz = head['dx'][2] # width of cell in z direction
        l_c = head['length_unit']
        t_c = head['time_unit']
        m_c = head['mass_unit']
        d_c = head['density_unit']                 
        v_c = head['velocity_unit']
        e_c = head['energy_unit']
        p_c = e_c # pressure units are the same as energy density units, density*velocity^2/length^3

        d  = f['density'][:]
        px  = f['momentum_x'][:]
        py  = f['momentum_y'][:]
        pz  = f['momentum_z'][:]

        f.close()

        n = d * d_c/(mu*mp) # number density, particles per cm^3 
        n_init = 1.0 

        km = 1e-5
        vx = (px*v_c*km)/d #velocity in the x direction

        dx = box_length / nx
        mass = d * dx*dx*dx
        cloud_mass = mass[n >= (n_init / cloud_thresh)]
        mass_tot = np.sum(cloud_mass)
        if i == 0 and j == 0:
             mass_init = mass_tot

        box_end = n[-1,:,:]

        v_avg = np.sum(vx[n >= (n_init / cloud_thresh)] * cloud_mass) / mass_tot
        mass_cur = mass_tot / mass_init

        with open(os.path.join(csvdir, filename), "a") as f_csv:
            writer_obj = writer(f_csv)
            if MASS:
                if np.any(box_end > (n_init/cloud_thresh)):
                    writer_obj.writerow([1, t, mass_cur])
                else:
                    writer_obj.writerow([0, t, mass_cur]) 
            if VELOCITY:
                if np.any(box_end > (n_init/cloud_thresh)):
                    writer_obj.writerow([1, t, v_avg])
                else:
                    writer_obj.writerow([0, t, v_avg])
            f_csv.close()

        f.close()    

    with open(os.path.join(csvdir, filename), "a") as f_csv:
            writer_obj = writer(f_csv)
            writer_obj.writerow([])
            f_csv.close()
