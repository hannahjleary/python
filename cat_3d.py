#!/usr/bin/env python3
# Example file for concatenating 3D hdf5 datasets

import h5py
import numpy as np

ns = 0
ne = 1
step = 10 # n_hydro
n_procs = 28 # number of processors that did the cholla calculation
dnamein = '../../../../../ix/eschneider/hjl28/data/cloud_wind/4/48retry/hdf5/raw/'
dnameout = '../../../../../ix/eschneider/hjl28/data/cloud_wind/4/48retry/hdf5/'

DE = 0

# loop over outputs
for n in range(ns, ne, step):

  # loop over files for a given output
  for i in range(0, n_procs):

    # open the output file for writing (don't overwrite if exists)
    fileout = h5py.File(dnameout+str(n)+'.h5', 'a')
    # open the input file for reading
    filein = h5py.File(dnamein+str(n)+'/'+str(n)+'_slice.h5.'+str(i), 'r')
    # read in the header data from the input file
    head = filein.attrs

    # if it's the first input file, write the header attributes 
    # and create the datasets in the output file
    if (i == 0):
      nx = head['dims'][0]
      ny = head['dims'][1]
      nz = head['dims'][2]
      dx = head['dx'][0]
      dy = head['dx'][1]
      dz = head['dx'][2]
      fileout.attrs['dx'] = [dx,dy,dz]
      fileout.attrs['dims'] = [nx, ny, nz]
      fileout.attrs['gamma'] = [head['gamma'][0]]
      fileout.attrs['t'] = [head['t'][0]]
      fileout.attrs['dt'] = [head['dt'][0]]
      fileout.attrs['n_step'] = [head['n_step'][0]]

      units = ['time_unit', 'mass_unit', 'length_unit', 'energy_unit', 'velocity_unit', 'density_unit']
      for unit in units:
        fileout.attrs[unit] = [head[unit][0]]

      d  = fileout.create_dataset("density", (nx, ny, nz), chunks=True)
      mx = fileout.create_dataset("momentum_x", (nx, ny, nz), chunks=True)
      my = fileout.create_dataset("momentum_y", (nx, ny, nz), chunks=True)
      mz = fileout.create_dataset("momentum_z", (nx, ny, nz), chunks=True)
      E  = fileout.create_dataset("Energy", (nx, ny, nz), chunks=True)
      if (DE):
        GE = fileout.create_dataset("GasEnergy", (nx, ny, nz), chunks=True)

    # write data from individual processor file to
    # correct location in concatenated file
    nxl = head['dims_local'][0]
    nyl = head['dims_local'][1]
    nzl = head['dims_local'][2]
    xs = head['offset'][0]
    ys = head['offset'][1]
    zs = head['offset'][2]
    fileout['density'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl]  = filein['density']
    fileout['momentum_x'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['momentum_x']
    fileout['momentum_y'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['momentum_y']
    fileout['momentum_z'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['momentum_z']
    fileout['Energy'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl]  = filein['Energy']
    if (DE):
      fileout['GasEnergy'][xs:xs+nxl,ys:ys+nyl,zs:zs+nzl] = filein['GasEnergy']
      
    filein.close()

  fileout.close()