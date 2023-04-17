# %%
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# %%
import os
endTime = max(os.listdir('postProcessing/incidentGraph'))

# Load inlet data
data1 = pd.read_csv('postProcessing/inletGraph/'+endTime+'/line_U.csv')
data2 = pd.read_csv('postProcessing/inletGraph/'+endTime+'/line_p_k_epsilon_nut.csv')

Z_inlet = data1['z']
U_inlet = data1['U_0']
k_inlet = data2['k']
epsilon_inlet = data2['epsilon']
nut_inlet = data2['nut']
p_inlet = data2['p']

# Load incident data
data1 = pd.read_csv('postProcessing/incidentGraph/'+endTime+'/line_U.csv')
data2 = pd.read_csv('postProcessing/incidentGraph/'+endTime+'/line_p_k_epsilon_nut.csv')

Z_incident = data1['z']
U_incident = data1['U_0']
k_incident = data2['k']
epsilon_incident = data2['epsilon']
nut_incident = data2['nut']
p_incident = data2['p']

# %%
# Function for calculating average inhomogeneity error
def calcError(gen, trg, weighted=True):

    avgError = (abs(gen-trg)*100/trg).mean()

    if weighted:

        # Calculate error weights by height
        Zw = (Z_inlet.diff(1))
        Zw[0] = Z_inlet[0]-0
        Zw /= Z_inlet.max()
        
        avgError = (Zw*(abs(gen-trg)*100/trg)).sum()

    return round(avgError,2)


# # %%
# # Extract yPlus value
# with open('postProcessing/yPlus/' + endTime + '/yPlus.dat') as f:
#     for line in f:
#         if line[0]!='#':
#             entry = line.split('\t')[1:]
#             if entry[0] == 'ground':
#                 yPlus = entry[1:]
#                 yPlus[-1] = yPlus[-1].split('\n')[0]
#                 yPlus = [round(float(x),0) for x in yPlus]

#%%

# Set plot limits
ylim_min = 0
ylim_max = Z_inlet.max()*1.05

# Plot y in log scale?
ylog = False

# Write error?
err = True


# Plot profiles.
fig, axes = plt.subplots(1, 4, figsize=(15, 4), dpi=100)

# Mean wind speed
axes[0].set(xlabel='$U$ $[m/s]$',ylabel='$z$ (m)', xlim=(0*U_inlet.min(), 1.1*U_inlet.max()), ylim=(ylim_min, ylim_max), title='Mean wind speed')
axes[0].plot(U_inlet, Z_inlet, 'k', marker = 'o', ms=3, lw=0.75)
axes[0].plot(U_incident, Z_incident, 'r', marker='o', ms=3, lw=0.75)
if err:
    axes[0].text(0.70, 0.05, '$e$ = {}%'.format(calcError(U_incident, U_inlet)), color='k', transform=axes[0].transAxes)

# Turbulence kinetic energy
axes[1].set(xlabel='$k$ $[m^2/s^2]$', xlim=(0.5*k_inlet.min(), 1.5*k_inlet.max()), ylim=(ylim_min, ylim_max),title='Turbulence kinetic energy')
axes[1].plot(k_inlet, Z_inlet, 'k', marker = 'o', ms=3, lw=0.75)
axes[1].plot(k_incident, Z_incident, 'r', marker='o', ms=3, lw=0.75)
if err:
    axes[1].text(0.70, 0.05, '$e$ = {}%'.format(calcError(k_incident, k_inlet)), color='k', transform=axes[1].transAxes)

# Turbulence dissipation rate
axes[2].set(xlabel='$\epsilon$ $[m^2/s^3]$', xlim=(0.8*epsilon_inlet.min(), 1.2*epsilon_inlet.max()), ylim=(ylim_min, ylim_max),title='Turbulence dissipation rate')
axes[2].semilogx(epsilon_inlet, Z_inlet, 'k', marker = 'o', ms=3, lw=0.75)
axes[2].semilogx(epsilon_incident, Z_incident, 'r', marker='o', ms=3, lw=0.75)
if err:
    axes[2].text(0.70, 0.25, '$e$ = {}%'.format(calcError(epsilon_incident, epsilon_inlet)), color='k', transform=axes[2].transAxes)

# Turbulence viscosity
axes[3].set(xlabel='$\\nu_T$ $[m^2/s$]', xlim=(0.8*nut_inlet.min(), 1.2*nut_inlet.max()), ylim=(ylim_min, ylim_max),title='Turbulence viscosity')
axes[3].plot(nut_inlet, Z_inlet, 'k', marker = 'o', ms=3, lw=0.75)
axes[3].plot(nut_incident, Z_incident, 'r', marker='o', ms=3, lw=0.75)
if err:
    axes[3].text(0.65, 0.25, '$e$ = {}%'.format(calcError(nut_incident, nut_inlet)), color='k', transform=axes[3].transAxes)


axes[0].legend(['Inlet','Incident'],loc='upper left', framealpha=0)

for axis in axes:
    axis.tick_params(axis="y",direction="in")
    axis.tick_params(axis="x",direction="in")
    if ylog:
        axis.set_yscale('log')
        axis.set(ylim=(Z_inlet.min()*0.75, Z_inlet.max()*1.2))


# # Add yPlus values to image:
# axes[3].text(0.45, 0.05, 'max $y^+$  ={:.2e}'.format(yPlus[1]), color='k', transform=axes[3].transAxes)
# axes[3].text(0.45, 0.1, 'mean $y^+$={:.2e}'.format(yPlus[2]), color='k', transform=axes[3].transAxes)
# axes[3].text(0.65, 0.15, '$y_P$ = {}m'.format(round(Z_inlet.min(),2)), color='k', transform=axes[3].transAxes)


plt.tight_layout()
plt.savefig('plot_ABL.jpg', dpi=300)