import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage

mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K

f = h5py.File('400_512_f32.h5', 'r')
head = f.attrs

head.keys()

f.keys()

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
p_c = e_c # pressure units are the same as energy density units, density*velocity^2/length^3

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

print(np.min(Vx), np.max(Vx))

mu = 1.0 # mean molecular weight (mu) of 1

d = d*d_c # to convert from code units to cgs, multiply by the code unit for that variable
n = d/(mu*mp) # number density, particles per cm^3

T = GE*(gamma - 1.0)*p_c / (n*kb) #Temperature

print(np.min(n), np.max(n))

## Masking
n_clouds = np.where(T<2e4, n, 10e-10)
n_int = np.where((T >= 2e4) & (T <= 5e5), n, 10e-10)
n_wind = np.where(T>5e5, n, 10e-10)

n_clouds_x = np.log10(np.sum(n_clouds, axis=0)*dx*l_c)
n_clouds_y = np.log10(np.sum(n_clouds, axis=1)*dx*l_c)
n_clouds_z = np.log10(np.sum(n_clouds, axis=2)*dx*l_c)

n_int_x = np.log10(np.sum(n_int, axis=0)*dx*l_c)
n_int_y = np.log10(np.sum(n_int, axis=1)*dx*l_c)
n_int_z = np.log10(np.sum(n_int, axis=2)*dx*l_c)

n_wind_x = np.log10(np.sum(n_wind, axis=0)*dx*l_c)
n_wind_y = np.log10(np.sum(n_wind, axis=1)*dx*l_c)
n_wind_z = np.log10(np.sum(n_wind, axis=2)*dx*l_c)

Vx_clouds = np.where(T<2e4, Vx, np.NaN)
Vx_int = np.where((T >= 2e4) & (T <= 5e5), Vx, np.NaN)
Vx_wind = np.where(T>5e5, Vx, np.NaN)

Vy_clouds = np.where(T<2e4, Vy, np.NaN)
Vy_int = np.where((T >= 2e4) & (T <= 5e5), Vy, np.NaN)
Vy_wind = np.where(T>5e5, Vy, np.NaN)

Vz_clouds = np.where(T<2e4, Vz, np.NaN)
Vz_int = np.where((T >= 2e4) & (T <= 5e5), Vz, np.NaN)
Vz_wind = np.where(T>5e5, Vz, np.NaN)

########################## Cloud ############################

fig = plt.figure(figsize=(16,4))
image = plt.imshow(n_clouds_x.T, origin='lower', cmap='viridis')
cb = plt.colorbar(image, label='$log_{10}(N_{H})$ [$cm^{-2}$]')

fig = plt.figure(figsize=(16,4))
image = plt.imshow(n_clouds_y.T, origin='lower', cmap='viridis')
cb = plt.colorbar(image, label='$log_{10}(N_{H})$ [$cm^{-2}$]')

#########################  ntermediate Gas ###########################33

fig = plt.figure(figsize=(16,4))
image = plt.imshow(n_int_x.T, origin='lower', cmap='viridis')
cb = plt.colorbar(image, label='$log_{10}(N_{H})$ [$cm^{-2}$]')

fig = plt.figure(figsize=(16,4))
image = plt.imshow(n_int_y.T, origin='lower', cmap='viridis')
cb = plt.colorbar(image, label='$log_{10}(N_{H})$ [$cm^{-2}$]')

#############################  Wind #################################

fig = plt.figure(figsize=(16,4))
image = plt.imshow(n_wind_x.T, origin='lower', cmap='viridis')
cb = plt.colorbar(image, label='$log_{10}(N_{H})$ [$cm^{-2}$]')

fig = plt.figure(figsize=(16,4))
image = plt.imshow(n_wind_y.T, origin='lower', cmap='viridis')
cb = plt.colorbar(image, label='$log_{10}(N_{H})$ [$cm^{-2}$]')


######################### Cloud Stats #########################################

print("X Direction Velocities")

print("Mean: ", np.nanmean(Vx_clouds), " km/s")
print("Median: ", np.nanmedian(Vx_clouds), " km/s")
print("75th Percentile: ", np.nanpercentile(Vx_clouds, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vx_clouds, 90), " km/s")

print("Y Direction Velocities")

print("Mean: ", np.nanmean(Vy_clouds), " km/s")
print("Median: ", np.nanmedian(Vy_clouds), " km/s")
print("75th Percentile: ", np.nanpercentile(Vy_clouds, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vy_clouds, 90), " km/s")

print("Z Direction Velocities")

print("Mean: ", np.nanmean(Vz_clouds), " km/s")
print("Median: ", np.nanmedian(Vz_clouds), " km/s")
print("75th Percentile: ", np.nanpercentile(Vz_clouds, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vz_clouds, 90), " km/s")

V_clouds = np.sqrt(np.power(Vx_clouds, 2) + np.power(Vy_clouds, 2) + np.power(Vz_clouds, 2))

print("Speeds")

print("Mean: ", np.nanmean(V_clouds), " km/s")
print("Median: ", np.nanmedian(V_clouds), " km/s")
print("75th Percentile: ", np.nanpercentile(V_clouds, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(V_clouds, 90), " km/s")

<h2>Intermediate Gas Stats</h2>

print("X Direction Velocities")

print("Mean: ", np.nanmean(Vx_int), " km/s")
print("Median: ", np.nanmedian(Vx_int), " km/s")
print("75th Percentile: ", np.nanpercentile(Vx_int, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vx_int, 90), " km/s")

print("Y Direction Velocities")

print("Mean: ", np.nanmean(Vy_int), " km/s")
print("Median: ", np.nanmedian(Vy_int), " km/s")
print("75th Percentile: ", np.nanpercentile(Vy_int, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vy_int, 90), " km/s")

print("Z Direction Velocities")

print("Mean: ", np.nanmean(Vz_int), " km/s")
print("Median: ", np.nanmedian(Vz_int), " km/s")
print("75th Percentile: ", np.nanpercentile(Vz_int, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vz_int, 90), " km/s")

V_int = np.sqrt(np.power(Vx_int, 2) + np.power(Vy_int, 2) + np.power(Vz_int, 2))

print("Speeds")

print("Mean: ", np.nanmean(V_int), " km/s")
print("Median: ", np.nanmedian(V_int), " km/s")
print("75th Percentile: ", np.nanpercentile(V_int, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(V_int, 90), " km/s")

<h2>Wind Stats</h2>

print("X Direction Velocities")

print("Mean: ", np.nanmean(Vx_wind), " km/s")
print("Median: ", np.nanmedian(Vx_wind), " km/s")
print("75th Percentile: ", np.nanpercentile(Vx_wind, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vx_wind, 90), " km/s")

print("Y Direction Velocities")

print("Mean: ", np.nanmean(Vy_wind), " km/s")
print("Median: ", np.nanmedian(Vy_wind), " km/s")
print("75th Percentile: ", np.nanpercentile(Vy_wind, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vy_wind, 90), " km/s")

print("Z Direction Velocities")

print("Mean: ", np.nanmean(Vz_wind), " km/s")
print("Median: ", np.nanmedian(Vz_wind), " km/s")
print("75th Percentile: ", np.nanpercentile(Vz_wind, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(Vz_wind, 90), " km/s")

V_wind = np.sqrt(np.power(Vx_wind, 2) + np.power(Vy_wind, 2) + np.power(Vz_wind, 2))

print("Speeds")

print("Mean: ", np.nanmean(V_wind), " km/s")
print("Median: ", np.nanmedian(V_wind), " km/s")
print("75th Percentile: ", np.nanpercentile(V_wind, 75), " km/s")
print("90th Percentile: ", np.nanpercentile(V_wind, 90), " km/s")