import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K
mu = 0.6

DE = 1 #Dual Energy Flag

dnamein='../../data/cloud_wind/3/test/hdf5/' # directory where the file is located

iend = 60
istart = 0
step = 20

for i in range(istart, iend, step):

    print('timestep: ' + str(i))

    f = h5py.File(dnamein + str(i) + '_slice.h5.0', 'r') 
    head = f.attrs # read the header attributes into a structure, called head

    t  = head['t'] # time of this snapshot, in kyr

    gamma = head['gamma']
    nx = head['dims'][0] # number of cells in the x direction
    ny = head['dims'][1] # number of cells in the y direction
    nz = head['dims'][2] # number of cells in the z direction
    dx = head['dx'][0] # width of cell in x direction

    v_c = head['velocity_unit']
    d_c = head['density_unit']
    e_c = head['energy_unit']
    p_c = e_c

    d = f['d_xy'][:]
    px  = f['mx_xy'][:]
    py  = f['my_xy'][:]
    pz  = f['mz_xy'][:]
    E = f['E_xy'][:]
    if DE:
        GE = f['GE_xy'][:]

    f.close()

    vx = px/d
    vy = py/d
    vz = pz/d
    if not DE:
        KE = 0.5 * d * (vx*vx + vy*vy + vz*vz)
        GE = E - KE
    km = 1e-5
    vx = vx*v_c*km #velocity in the x direction

    n = d * d_c/ (mu*mp) # number density, particles per cm^3  
    T = GE*(gamma-1.0)*p_c / (n*kb) #temperature
    logT = np.log10(T)

    for k in range(len(vx)):
        if k == 0:
            print("Boundary Cell: " )
        else:
            print("x index " + str(k) + ": " )
        print("x Momentum: " + str(px[k]))
        print("Energy: " + str(E[k]))
        print("Density: " + str(d[k]))
        print("\n")

    # Sound speed in the wind
    if i == 0 and j==0:
        windT = T[0,0]
        windn = n[0,0]
        windv = vx[0,0]
        cells = ny
        P = windT * windn * kb / p_c
        rho = windn * mu * mp / d_c
        c = np.sqrt(gamma*P/ rho)
        c = c * v_c * 1e-5
        print("Sound Speed: " + str(int(c)))
        print("Wind speed: " + str(int(windv)))
        if int(windv) > int(c):
            speed = 'Supersonic Wind'
        elif int(windv) == int(c):
            speed = 'Transsonic Wind'
        else:
            speed = 'Subsonic Wind'
        print(speed)
    #Mach Number
        M = windv / c
        print("Mach Number: " + str(np.round(float(M), 2)))
