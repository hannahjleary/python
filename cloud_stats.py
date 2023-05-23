import h5py
import numpy as np
import matplotlib
#matplotlib.use('Agg')
matplotlib.rcParams['mathtext.default']='regular'
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage


# define some constants
mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K


istart = 200 # starting file number
iend   = 200 # ending file number
t_cc = 39.8 # cloud crushing time in kyr, n = 0.5 cloud
#t_cc = 56.4 # cloud crushing time in kyr, n = 1 cloud
dnamein='../data/' # directory where the file is located
dnameout='../plots/' # directory where the plot will be saved

# begin loop over input files
for i in range(1):

  print(i)
  f = h5py.File(dnamein + '200_512_f32.h5', 'r') # open the hdf5 file for reading
  head = f.attrs # read the header attributes into structure, head
  gamma = head['gamma'] # ratio of specific heats
  t  = head['t'] # time of this snapshot, in kyr
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
  p_c = e_c # pressure units are the same as energy, density*velocity^2
  
  
  d  = f['density'][:]
  GE = f['GasEnergy'][:]
  px  = f['momentum_x'][:]
  py = f['momentum_y'][:]
  pz = f['momentum_z'][:]

  f.close()

  km = 1e-5

  Vx = (px*v_c*km)/d #velocity in the x direction
  Vy = (py*v_c*km)/d #velocity in the y direction
  Vz = (pz*v_c*km)/d #velocity in the z direction

  mu = 1.0 # mean molecular weight (mu) of 1

  # for these simulations I used mu = 1 and a code mass unit of 1 mp,
  # so the following conversion from density to number density isn't really necessary,
  # but I'm including it for completeness
  d = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable
  n = d/(mu*mp) # number density, particles per cm^3
  T = GE*(gamma - 1.0)*p_c / (n*kb)

  d_avg = np.average(d)

  n_clouds = np.where(T<2e4, np.log10(n), 0)
  com = ndimage.center_of_mass(n_clouds)

  #Make array of density weighted distances from center of mass
  cloud_indices = np.where(T<2e4)
  data = np.array(zip(cloud_indices[0], cloud_indices[1], cloud_indices[2]))
  cloud_coords = list(data[()]) 

  #print(cloud_coords[1])
  print(com)
  #print(cloud_coords)

  print(np.shape(cloud_coords))
  print(np.shape(com))
  dist = np.abs(np.subtract(cloud_coords[:], com))
  
  print(dist)

  dist_weighted = dist*d[cloud_coords]/d_avg

  #mean_dist = np.nanmean(dist_weighted)
  #std_dist = np.nanstd(dist_weighted)

  #print(i)
  #print("c.o.m.: " + str(com))
  #print("Mean distance from c.o.m.: " + str(mean_dist))
  #print("Standard dev. distance from c.o.m.: " + str(std_dist))



  ##Before zipping, you have three arrays, the x coords, y coords, and z coords
  ##Keep that and then convert them into distances from center of mass
  ##Density weight... you need coords then in order to access the density at
  ## that specific coord. Maybe unzip after this.

  #is this even useful
