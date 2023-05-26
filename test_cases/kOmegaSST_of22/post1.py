# %%
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# %%
import os
endTime = max(os.listdir('postProcessing/samples_u'))

U = {}
k = {}
omega = {}
nut = {}
Z = {}

for i in np.arange(0,2200,200):
    U.update({"x"+str(i):pd.read_csv('postProcessing/samples_u/'+endTime+'/x_'+str(i)+'_U.csv')['U_0']})
    k.update({"x"+str(i):pd.read_csv('postProcessing/samples_k/'+endTime+'/x_'+str(i)+'_k.csv')['k']})
    omega.update({"x"+str(i):pd.read_csv('postProcessing/samples_omega/'+endTime+'/x_'+str(i)+'_omega.csv')['omega']})
    nut.update({"x"+str(i):pd.read_csv('postProcessing/samples_nut/'+endTime+'/x_'+str(i)+'_nut.csv')['nut']})
    Z.update({"x"+str(i):pd.read_csv('postProcessing/samples_u/'+endTime+'/x_'+str(i)+'_U.csv')['z']})

#%%
Z_inlet = Z["x0"]
U_inlet = U["x0"]
k_inlet = k["x0"]
omega_inlet = omega["x0"]
nut_inlet = nut["x0"]

# Set plot limits
ylim_min = 0
ylim_max = Z["x0"].max()*1.05

# Plot y in log scale?
ylog = False

# Write error?
err = True


# Plot profiles.
fig, axes = plt.subplots(1, 4, figsize=(15, 4), dpi=100)

# Mean wind speed
axes[0].set(xlabel='$U$ $[m/s]$',ylabel='$z$ (m)', xlim=(0,8), ylim=(ylim_min, ylim_max), title='Mean wind speed')
for x in U:
    axes[0].plot(U[x], Z[x], label='x='+x[1:], lw=0.75)


# Turbulence kinetic energy
axes[1].set(xlabel='$k$ $[m^2/s^2]$', xlim=(0.5*k_inlet.min(), 1.5*k_inlet.max()), ylim=(ylim_min, ylim_max),title='Turbulence kinetic energy')
for x in k:
    axes[1].plot(k[x], Z[x], label='x='+x[1:], lw=0.75)

# Turbulence dissipation rate
axes[2].set(xlabel='$\omega$ $[1/s]$', xlim=(0.8*omega_inlet.min(), 1.2*omega_inlet.max()), ylim=(ylim_min, ylim_max),title='Specific dissipation rate')
for x in omega:
    axes[2].semilogx(omega[x], Z[x], label='x='+x[1:], lw=0.75)


# Turbulence viscosity
axes[3].set(xlabel='$\\nu_T$ $[m^2/s$]', xlim=(0.8*nut_inlet.min(), 1.2*nut_inlet.max()), ylim=(ylim_min, ylim_max),title='Turbulence viscosity')
for x in omega:
    axes[3].plot(nut[x], Z[x], label='x='+x[1:], lw=0.75)

# axes[0].legend(['Inlet','Incident'],loc='upper left', framealpha=0)
axes[0].legend()

for axis in axes:
    axis.tick_params(axis="y",direction="in")
    axis.tick_params(axis="x",direction="in")
    if ylog:
        axis.set_yscale('log')
        axis.set(ylim=(Z_inlet.min()*0.75, Z_inlet.max()*1.2))


plt.tight_layout()
plt.savefig('plot1_ABL.jpg', dpi=300)


# %%
