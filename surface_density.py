import h5py
import numpy as np
import matplotlib
#matplotlib.use('Agg')
matplotlib.rcParams['mathtext.default']='regular'
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
  f.close()
  mu = 1.0 # mean molecular weight (mu) of 1
  r = ny/8 # these simulations had an initial cloud radius of 1/8 of the size of the domain in the y-z directions

  # for these simulations I used mu = 1 and a code mass unit of 1 mp,
  # so the following conversion from density to number density isn't really necessary,
  # but I'm including it for completeness
  d = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable
  n = d/(mu*mp) # number density, particles per cm^3
  T = GE*(gamma - 1.0)*p_c / (n*kb)

  print(np.min(T), np.max(T))

  # to create a density projection, integrate the 3D density array along one axis
  pn_x = np.sum(n, axis=0)*dx*l_c
  pn_y = np.sum(n, axis=1)*dx*l_c

  # set the column density scale
  log_pn_x = np.log10(pn_x)
  log_pn_y = np.log10(pn_y)
  print("n range    = %e %e" % (np.min(n),np.max(n)))
  print("N range    = %5.2f %5.2f" % (np.min(log_pn_x), np.max(log_pn_x)))
  pn_min=17.5
  pn_max=21.0


  # plot the surface density, x-axis projection
  fig, ax = plt.subplots(figsize=(4,3))
  ax.set_xticks(ny*np.arange(0.25, 1, 0.25))
  ax.set_yticks(nz*np.arange(0.25, 1, 0.25))
  ax.tick_params(axis='both', which='both', direction='in', color='white', labelleft=0, labelbottom=0, top=1, right=1)
  image = ax.imshow(log_pn_x.T, origin='lower', cmap='viridis', vmin=pn_min, vmax=pn_max)

  # add a circle showing the original extent of the cloud
  circle=plt.Circle((ny/2,nz/2),r,fill=False,edgecolor='white',linestyle='dashed',linewidth=0.5)
  fig.gca().add_artist(circle)

  # add a colorbar
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(image, cax = cax, ticks=np.arange(17.5, 21.5, 0.5))
  cax.tick_params(axis='y', direction='in')
  cb.solids.set_edgecolor('face')
  cax.set_ylabel('$log_{10}(N_{H})$ [$cm^{-2}$]')
  #plt.show()

  # save the figure
  plt.savefig(dnameout + 'sd_x_'+str(i)+'_1024.png', dpi=300)
  plt.close(fig)


  # plot the surface density, y-axis projection
  fig = plt.figure(figsize=(8,2), frameon=False)
  a0 = fig.add_axes([0,0,1,1])
  a0.set_xticks(nx*np.arange(0.0625, 1, 0.0625))
  a0.set_yticks(ny*np.arange(0.25, 1, 0.25))
  a0.tick_params(axis='both', which='both', color='white', direction='in', labelleft=0, labelbottom=0, top=1, right=1)
  image = a0.imshow(log_pn_y.T, origin='lower', cmap='viridis', vmin=pn_min, vmax=pn_max)
  a0.hlines(0.15*ny, 2*r, 4*r, color='white')
  a0.text(4.1*r, r, '10 pc', color='white') 
  #plt.xlim(0, 0.75*nx)
  a0.text(2*r, 0.8*ny, r'$\tilde{n} = 0.5$ cloud ', color='white')
  a0.text(0.85*nx, 0.8*ny, str(int(t/t_cc))+r'$t_{cc}$', color='white')
  #plt.show()

  # save the figure
  plt.savefig(dnameout + 'sd_y_'+str(i)+'_512.png', dpi=300)
  plt.close(fig)
