#!/usr/bin/env python3
# Example file for concatenating on-axis projection data
# created when the -DPROJECTION flag is turned on

import h5py
import numpy as np

ns = 0
ne = 1
# step = 10 # n_hydro
n_procs = 8 # number of processors that did the cholla calculation
dnamein = '../../../../../ix/eschneider/hjl28/data/cloud_wind/4/48retry/hdf5/raw/'
dnameout = '../../../../../ix/eschneider/hjl28/data/cloud_wind/4/48retry/hdf5/'

# loop over the output times
for n in range(ns, ne):

  # open the output file for writing
  fileout = h5py.File(dnameout+str(n)+'_proj.h5', 'w')

  # loop over files for a given output time
  for i in range(0, n_procs):

    # open the input file for reading
    filein = h5py.File(dnamein+str(n)+'/'+str(n)+'_proj.h5.'+str(i), 'r')
    # read in the header data from the input file
    head = filein.attrs
    # if it's the first input file, write the header attributes
    # and create the datasets in the output file
    if (i == 0):
      gamma = head['gamma']
      t = head['t']
      dt = head['dt']
      dx = head['dx']
      n_step = head['n_step']
      nx = head['dims'][0]
      ny = head['dims'][1]
      nz = head['dims'][2]
      length_unit = head["length_unit"]
      mass_unit = head["mass_unit"]
      time_unit = head["time_unit"]
      density_unit = head["density_unit"]
      velocity_unit = head["velocity_unit"]
      energy_unit = head["energy_unit"]

      fileout.attrs['gamma'] = gamma
      fileout.attrs['t'] = t
      fileout.attrs['dt'] = dt
      fileout.attrs['dx'] = dx
      fileout.attrs['n_step'] = n_step
      fileout.attrs['dims'] = [nx, ny, nz]
      fileout.attrs['length_unit'] = length_unit
      fileout.attrs['time_unit'] = time_unit
      fileout.attrs['mass_unit'] = mass_unit
      fileout.attrs['density_unit'] = density_unit
      fileout.attrs['velocity_unit'] = velocity_unit
      fileout.attrs['energy_unit'] = energy_unit

      dxy = np.zeros((nx,ny))
      dxz = np.zeros((nx,nz))
      Txy = np.zeros((nx,ny))
      Txz = np.zeros((nx,nz))

    # write data from individual processor file to
    # correct location in concatenated file
    nxl = head['dims_local'][0]
    nyl = head['dims_local'][1]
    nzl = head['dims_local'][2]
    xs = head['offset'][0]
    ys = head['offset'][1]
    zs = head['offset'][2]

    dxy[xs:xs+nxl,ys:ys+nyl] += filein['d_xy']
    dxz[xs:xs+nxl,zs:zs+nzl] += filein['d_xz']
    Txy[xs:xs+nxl,ys:ys+nyl] += filein['T_xy']
    Txz[xs:xs+nxl,zs:zs+nzl] += filein['T_xz']

    filein.close()

  # write out the new datasets
  fileout.create_dataset('d_xy', data=dxy)
  fileout.create_dataset('d_xz', data=dxz)
  fileout.create_dataset('T_xy', data=Txy)
  fileout.create_dataset('T_xz', data=Txz)

  fileout.close()