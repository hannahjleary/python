import h5py
import numpy as np
import matplotlib
#matplotlib.use('Agg')
matplotlib.rcParams['mathtext.default']='regular'
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

iend = 100 #time
dnamein = '../data/'
dnameout = '../plots/'

DE = 0 #Dual Energy Flag

mp = 1.672622e-24   # mass of hydrogren atom (grams)
kb = 1.38065e-16    # boltzmann constant (ergs/K)

for i in range(iend):

      f = h5py.File(dnamein + str(i) + '.h5.0', 'r') #200_512_f32.h5
      head = f.attrs

      gamma = head['gamma'] #ratio of specific heats
      t = head['t'] #time of snapshot (kyr)
      nx = head['dims'][0] # number of cels in he x direction
      ny = head['dims'][1] # number of cells in the y direction
      nz = head['dims'][2] # number of cells in the z direction
      dx = head['dx'][0] #width of cells in the x direction
      dy = head['dx'][1] #width of cells in the y direction
      dz = head['dx'][2] #width of cells in the z direction
      l_c = head['length_unit']
      t_c = head['time_unit']
      m_c = head['mass_unit']
      d_c = head['density_unit']
      v_c = head['velocity_unit']
      e_c = head['energy_unit']
      p_c = e_c

      d = f['density'][:]
      E = f['Energy'][:]
      px = f['momentum_x'][:]
      py = f['momentum_y'][:]
      pz = f['momentum_z'][:]
      if DE:
            GE = f['GasEnergy'][:]
            f.close()
      else:
            f.close()
            vx = px/d
            vy = py/d
            vz = pz/d

            mu = 1.0

            n = d * d_c/ (mp * mu)

            KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
            GE = E - KE
            
      T = GE * (gamma - 1.0) * p_c / (n * kb) 
      T_log = np.log10(T)

      Tslice_xz = T_log[:,int(ny/2),:]

      fig, ax = plt.subplots()

      ax.set_yticks(nz*np.arange(0.25, 1, 0.25))
      ax.tick_params(axis='both', which='both', direction='in', color='black', bottom=1, left=1, top=0, right=0, labelleft=0, labelbottom=0, labeltop=0, labelright=0)
      ax.text(6, 6, 'y: ' + str(i+1), size=16, color='black')

      image = plt.imshow(Tslice_xz.T, origin ='lower', cmap='plasma')

      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.2)
      cb = plt.colorbar(image, cax=cax, orientation='vertical',ticks=np.arange(2.5, 7.5, 0.5))
      cax.tick_params(axis='y', direction='in')
      cb.solids.set_edgecolor('face')
      cax.set_ylabel('$log_{10}(K)$') 

      # plt.tight_layout()
      
      plt.savefig(dnameout + 't_xz_' + str(i) + '.png', dpi=300)
      plt.close(fig)
