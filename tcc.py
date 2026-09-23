import matplotlib.pyplot as plt
import numpy as np
import h5py
import seaborn as sns

pc = 3.086e16 # meters per pc
year = 3.154e7 # seconds per year
gamma = 1.66667
mu = 0.6
mp = 1.672622e-24 # mass of hydrogren atom, in grams
kb = 1.380658e-16 # boltzmann constant in ergs/K

P = 501
T_w = 501000
Chi = 100
Rcl = 10 #pc
Rcl_m = Rcl * pc 
v_wind = 100 * 1000 #m/s

c = np.sqrt((gamma * kb * T_w)/(mu * mp)) * 1e-5

print('Wind sound speed: ' + str(c) + " km/s")

t_cc = Chi**(1/2) * Rcl_m / v_wind
t_cc = t_cc / year / 1000

print("Cloud Crushing Time: " + str(t_cc) + " kyr")

