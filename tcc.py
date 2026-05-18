import matplotlib.pyplot as plt
import numpy as np
import h5py
import seaborn as sns

pc = 3.086e16 # meters per pc
year = 3.154e7 # seconds per year

Chi = 100
Rcl = 20 #pc
Rcl_m = Rcl * pc 
v_wind = 100 * 1000 #m/s

t_cc = Chi**(1/2) * Rcl_m / v_wind
t_cc = t_cc / year / 1000

print("Cloud Crushing Time: " + str(t_cc) + " kyr")

