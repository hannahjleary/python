import matplotlib.pyplot as plt
import numpy as np
import h5py
import seaborn as sns

## Plots cooling curves

DARKMODE = 1

if DARKMODE:
  bg_color='black'
  text_color='white'
  spine_color='white'
else:
  bg_color='white'
  text_color='black'
  spine_color='black'

dnamein = "../../../../../ix/eschneider/hjl28/data/cooling_tables/"
dnameout = "../../../../../ix/eschneider/hjl28/plots/" 

###### Read in Grackle library rates #############################################################################

f = h5py.File(dnamein + 'CloudyData_UVB=HM2012.h5','r')
print(list(f.keys()))
print(list(f['CoolingRates'].keys()))
print(list(f['CoolingRates']['Metals'].keys()))
metals_cool = f['CoolingRates']['Metals']['Cooling']
primordial_cool = f['CoolingRates']['Primordial']['Cooling']

print(list(metals_cool.attrs.keys()))

print("Axis 0:", metals_cool.attrs['Parameter1_Name'], metals_cool.attrs['Parameter1'].shape)
print("Axis 1:", metals_cool.attrs['Parameter2_Name'], metals_cool.attrs['Parameter2'].shape)
print("Axis 2:", "Temperature", metals_cool.attrs['Temperature'].shape)
print("min: ", np.min(metals_cool.attrs['Parameter2']), "max: ", np.max(metals_cool.attrs['Parameter2']))
print("shape: ", metals_cool.shape)

metals_temps = metals_cool.attrs['Temperature']
# print(metals_cool.attrs['Parameter2'])
prim_temps = primordial_cool.attrs['Temperature'] # index 20 is log(hden)=0

metals_cool = metals_cool[20,0,:]
prim_cool = primordial_cool[20,0,:]

f.close()

###### Read in Cloudy test data  #############################################################################

T_new = []
Lam_new = []
with open("../cloudy/c25.00/tsuite/programs/hazy_coolingcurve/hazy_coolingcurve.txt") as f2:
  i=0
  for line in f2:
    if i==0:
      i+=1
      continue
    T_new.append(float(line.split()[0]))
    Lam_new.append(float(line.split()[1]))
    i+=1

T_new = np.array(T_new)
# print(T_new.shape)
Lam_new = np.array(Lam_new)

# print(T_new)
# print(Lam_new)

###### Read in Wiersma KI cooling table #############################################################################

f = h5py.File(dnamein + 'z_0.000.hdf5','r')
T = np.array(f['Solar/Temperature_bins/'])
L = np.array(f['Solar/Net_cooling/'])
f.close()

###### Parabolic Fit  #############################################################################

H = 10**np.arange(0, 4, 0.01)
T1 = 10**np.arange(0, 4, 0.01)
T2 = np.arange(4, 9, 0.01)
cool = np.zeros(np.size(T2))
cool1 = 2e-26 * (1e7 * np.exp(-1.148e5 / (T1+1000)) + 1.4e-2 * np.sqrt(T1) * np.exp(-92.0/T1))

i = 0
for temp in T2:
#  if (temp > 4.0 and temp < 4.43):
#    cool[i] = np.power(10.0, (-15.0 * (temp - 4.30) * (temp - 4.30) - 21.85))
  if (temp > 4.0 and temp < 5.9):
    cool[i] = np.power(10.0, (-1.3 * (temp - 5.25) * (temp - 5.25) - 21.25))
  if (temp > 5.9 and temp < 7.4):
    cool[i] = np.power(10.0, (0.7 * (temp - 7.1) * (temp - 7.1) - 22.8))
  if (temp > 7.4):
    cool[i] = np.power(10.0, (0.45*temp - 26.065))
  i += 1

LK1 = 1e-24
al1 = 2.3
Tk1 = 8e3
c1 = LK1*np.power(10**T2 / Tk1, al1)
LK2 = 6e-22
al2 = -0.65
Tk2 = 1e5
c2 = LK2*np.power(10**T2 / Tk2, al2)
LK3 = 1.75e-23
al3 = 0.45
Tk3 = 2e7
c3 = LK3*np.power(10**T2 / Tk3, al3)
Tref = 1e9
Lref = 1e-24
al = -0.75
c = Lref*np.power((10**T2/Tref), al)

####### Plot Cholla primordial cooling function (Katz 1996) #######################

def primordial_cholla(n, T):
    # Helium abundance by mass
    Y = 0.24
    y = Y / (4 - 4 * Y)

    # Hydrogen number density
    n_h = n

    # Recombination and collisional ionization rates (Table 2, Katz 1996)
    alpha_hp   = (8.4e-11) * (1.0 / np.sqrt(T)) * (T / 1e3)**(-0.2) * (1.0 / (1.0 + (T / 1e6)**0.7))
    alpha_hep  = (1.5e-10) * (T**(-0.6353))
    alpha_d    = (1.9e-3) * (T**(-1.5)) * np.exp(-470000.0 / T) * (1.0 + 0.3 * np.exp(-94000.0 / T))
    alpha_hepp = (3.36e-10) * (1.0 / np.sqrt(T)) * (T / 1e3)**(-0.2) * (1.0 / (1.0 + (T / 1e6)**0.7))
    gamma_eh0  = (5.85e-11) * np.sqrt(T) * np.exp(-157809.1 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5)))
    gamma_ehe0 = (2.38e-11) * np.sqrt(T) * np.exp(-285335.4 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5)))
    gamma_ehep = (5.68e-12) * np.sqrt(T) * np.exp(-631515.0 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5)))

    # Photoionization rates (assuming J(nu) = 10^-22 (nu_L/nu))
    gamma_lh0  = 3.19851e-13
    gamma_lhe0 = 3.13029e-13
    gamma_lhep = 2.00541e-14

    # Heating rates
    e_h0  = 2.4796e-24
    e_he0 = 6.86167e-24
    e_hep = 6.21868e-25

    # Solve for number densities (no photoionization)
    heat_flag = False

    n_e = n_h  # initial guess
    if heat_flag:
        n_iter = 20
        tol = 1.0e-6
        for _ in range(n_iter):
            n_e_old = n_e
            n_h0    = n_h * alpha_hp / (alpha_hp + gamma_eh0 + gamma_lh0 / n_e)
            n_hp    = n_h - n_h0
            n_hep   = y * n_h / (1.0 + (alpha_hep + alpha_d) / (gamma_ehe0 + gamma_lhe0 / n_e)
                                      + (gamma_ehep + gamma_lhep / n_e) / alpha_hepp)
            n_he0   = n_hep * (alpha_hep + alpha_d) / (gamma_ehe0 + gamma_lhe0 / n_e)
            n_hepp  = n_hep * (gamma_ehep + gamma_lhep / n_e) / alpha_hepp
            n_e     = n_hp + n_hep + 2 * n_hepp
            if abs(n_e_old - n_e) < tol:
                break
    else:
        n_h0   = n_h * alpha_hp / (alpha_hp + gamma_eh0)
        n_hp   = n_h - n_h0
        n_hep  = y * n_h / (1.0 + (alpha_hep + alpha_d) / gamma_ehe0 + gamma_ehep / alpha_hepp)
        n_he0  = n_hep * (alpha_hep + alpha_d) / gamma_ehe0
        n_hepp = n_hep * gamma_ehep / alpha_hepp
        n_e    = n_hp + n_hep + 2 * n_hepp

    # Cooling rates for various processes (Table 1, Katz 1996)
    le_h0   = (7.50e-19) * np.exp(-118348.0 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5))) * n_e * n_h0
    le_hep  = (5.54e-17) * T**(-0.397) * np.exp(-473638.0 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5))) * n_e * n_hep
    li_h0   = (1.27e-21) * np.sqrt(T) * np.exp(-157809.1 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5))) * n_e * n_h0
    li_he0  = (9.38e-22) * np.sqrt(T) * np.exp(-285335.4 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5))) * n_e * n_he0
    li_hep  = (4.95e-22) * np.sqrt(T) * np.exp(-631515.0 / T) * (1.0 / (1.0 + np.sqrt(T / 1e5))) * n_e * n_hep
    lr_hp   = (8.70e-27) * np.sqrt(T) * (T / 1e3)**(-0.2) * (1.0 / (1.0 + (T / 1e6)**0.7)) * n_e * n_hp
    lr_hep  = (1.55e-26) * T**0.3647 * n_e * n_hep
    lr_hepp = (3.48e-26) * np.sqrt(T) * (T / 1e3)**(-0.2) * (1.0 / (1.0 + (T / 1e6)**0.7)) * n_e * n_hepp
    ld_hep  = (1.24e-13) * T**(-1.5) * np.exp(-470000.0 / T) * (1.0 + 0.3 * np.exp(-94000.0 / T)) * n_e * n_hep
    g_ff    = 1.1 + 0.34 * np.exp(-(5.5 - np.log(T))**2 / 3.0)  # Gaunt factor
    l_ff    = (1.42e-27) * g_ff * np.sqrt(T) * (n_hp + n_hep + 4 * n_hepp) * n_e

    # Total cooling rate (erg s^-1 cm^-3)
    cool = le_h0 + le_hep + li_h0 + li_he0 + li_hep + lr_hp + lr_hep + lr_hepp + ld_hep + l_ff

    return cool/(n**2)

cholla_cool = []
for tp in prim_temps:
  cholla_cool.append(primordial_cholla(1, tp))
cholla_cool = np.array(cholla_cool)

####### Calculate cooling rate and cooling time at a given temp ################

T_mix = np.sqrt(10**4*10**6)
T_mix_log = np.log10(T_mix)
n_mix = np.sqrt(0.1*1e-3)
T_wind = 10**6
n_wind = 1e-2
# print(T_mix_log)
# print(n_mix)
# print(prim_temps)
# print(T2)
# coolrate=cool[np.argmin(np.abs(T2-T_mix_log))]
coolrate=cool[np.argmin(np.abs(T2-np.log10(T_mix)))]
lam = coolrate * n_mix**2
print(np.min(np.abs(T2-6.0)))
print(lam)
coolratecheck = 1.114295e-26
lamcheck = 1.114295e-26
kb = 1.380658e-16 # ergs/K
year = 3.154e7 # sec in a year
energy = n_wind* kb * T_wind
tcool = (n_mix * kb * T_mix) / (lam)
# tcoolcheck = (n_wind * kb * T_wind) / (lamcheck)

print(str(tcool/year/1000) + " kyr")
# print(str(tcoolcheck/year/1000) + " kyr")

############ Plotting #######################################

palette = sns.color_palette('Blues')
palette2 = sns.color_palette('Reds')

fig = plt.figure(figsize=(4,3), dpi=300)
ax = fig.add_axes([0.18,0.17,0.75,0.75])
ax.set_facecolor(bg_color)
# ax.set_xticks(np.linspace(0,nx,9))
# ax.set_yticks(np.linspace(0,nz,5))
ax.tick_params(axis='both', which='both', direction='in', color=spine_color, 
              left=1, right=1, top=1, bottom=1, labelcolor=text_color, labelbottom=1, labelleft=1)
plt.setp(ax.spines.values(), color=spine_color)
plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=spine_color)
line1, = ax.plot(T[101:], L[101:,80], color=palette[4], label='Solar Abundance CIE', zorder=2)
line2, = ax.plot(T1, cool1, color=palette[2], linestyle='dashed', label='KI 02 Fit', zorder=1)
line3, = ax.plot(10**T2, cool, color=palette[1], label='Parabolic Fit', zorder=1)
# line4, = ax.plot(10**T_new, 10**Lam_new, color=palette[3], label='Cloudy Test', zorder=1) #green
# line5, = ax.plot(prim_temps, cholla_cool, color='yellow', linewidth=1, label='Cholla Primordial', zorder=1)
# line6, = ax.plot(prim_temps, prim_cool, color=palette[1], label='Primordial', zorder=1)
# line7, = ax.plot(metals_temps, metals_cool, color=palette[2], label='Metals', zorder=1)
# line8, = ax.plot(prim_temps, prim_cool+metals_cool, color=palette[3], label='Total', zorder=1)
ax.scatter(T_mix, cool[np.argmin(np.abs(T2-np.log10(T_mix)))])
#ax.plot(10**T2, c1, color="Red") 
#ax.plot(10**T2, c2, color="Red") 
#ax.plot(10**T2, c3, color="Red") 
#ax.plot(10**T2, c, color="Cyan")
# leg = plt.legend(fontsize=6, facecolor='black', loc='lower right') 
# for text in leg.get_texts(): text.set_color('white')
ax.set_xlabel('T [K]', color=text_color)
ax.set_ylabel('$\Lambda/n^2$ [erg $s^{-1}$ $cm^3$]', color=text_color)
plt.axis([1e1, 1e9, 1e-28, 1e-20])
plt.xscale('log')
plt.yscale('log')
plt.show()
fig.savefig(dnameout + 'cooling.png', dpi=300, facecolor=bg_color)




