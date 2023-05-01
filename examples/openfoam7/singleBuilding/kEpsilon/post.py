# %%
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
# %matplotlib inline


# %%
# Load experiment data
df_exp = pd.read_csv('exp.csv', index_col='Point ID', float_precision='round_trip')

# %%
# Load cfd data

Uref = 5.81 # Reference velocity at z=0.5m in wind tunnel test CEDVAL A1-1 (Parente et al 2011)

# Read CFD point coords
df_cfd = pd.read_csv('postProcessing/probes/0/Ucoords.csv', float_precision='round_trip')

# Read mean wind speed
df_uvw = pd.read_csv('postProcessing/probes/0/U.csv',header=[0,1], index_col=[0])
df_uvw = df_uvw.iloc[-1]

# Calculate normalised mean velocity
df_U = pd.DataFrame([])
probelist = df_uvw.index.get_level_values(0).unique()

for probe in probelist:
    # df_U[probe] = np.array([np.sqrt(df_uvw[probe]['U']**2 + df_uvw[probe]['V']**2 + df_uvw[probe]['W']**2)/Uref])
    df_U[probe] = np.array([df_uvw[probe]['U']/Uref])

df_U = df_U.T.values

# Read and normalise TKE
df_k = pd.read_csv('postProcessing/probes/0/k.csv', header=[0,1], index_col=[0])
df_k.columns = df_k.columns.droplevel(1)
df_k = df_k.iloc[-1].T
df_k = df_k/(Uref**2)

df_k = df_k.values

# Append to cfd dataframe
df_cfd['U/Uref'] = df_U
df_cfd['k/Uref2'] = df_k
df_cfd = df_cfd.set_index('Probe ID')

# %%
def plotdata(data, xline, var):
    '''
    Function to select location and variable of plot
    '''
    data1 = data[data['X']==xline]
    return data1[var].values, data1['Z'].values

# %%
# Plot U profiles
fig, axes = plt.subplots(2, 3, figsize=(8, 8), dpi=100)

var = 'U/Uref'

xline=-0.072
axes[0,0].set(ylabel='$z$ (m)', title='x='+str(xline)+'m')
axes[0,0].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[0,0].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=-0.04
axes[0,1].set(title='x='+str(xline)+'m')
axes[0,1].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[0,1].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0
axes[0,2].set(title='x='+str(xline)+'m')
axes[0,2].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[0,2].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0.04
axes[1,0].set(xlabel='$U/U_{ref}$',ylabel='$z$ (m)', title='x='+str(xline)+'m')
axes[1,0].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[1,0].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0.105
axes[1,1].set(xlabel='$U/U_{ref}$', title='x='+str(xline)+'m')
axes[1,1].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[1,1].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0.3
axes[1,2].set(xlabel='$U/U_{ref}$', title='x='+str(xline)+'m')
axes[1,2].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[1,2].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

axes[0,0].legend(['CFD','Exp'],loc='upper left', framealpha=1)

for axis1 in axes:
    for axis in axis1:
        axis.tick_params(axis="y",direction="in")
        axis.tick_params(axis="x",direction="in")
        axis.set(ylim=([0,0.25]), xlim=([-0.5,1]))

plt.tight_layout()
plt.savefig('plot_U.jpg', dpi=300)

# %%
# Plot K profiles
fig, axes = plt.subplots(2, 3, figsize=(8, 8), dpi=100)

var = 'k/Uref2'

xline=-0.072
axes[0,0].set(ylabel='$z$ (m)', title='x='+str(xline)+'m')
axes[0,0].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[0,0].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=-0.04
axes[0,1].set(title='x='+str(xline)+'m')
axes[0,1].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[0,1].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0
axes[0,2].set(title='x='+str(xline)+'m')
axes[0,2].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[0,2].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0.04
axes[1,0].set(xlabel='$k/U^2_{ref}$',ylabel='$z$ (m)', title='x='+str(xline)+'m')
axes[1,0].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[1,0].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0.105
axes[1,1].set(xlabel='$k/U^2_{ref}$', title='x='+str(xline)+'m')
axes[1,1].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[1,1].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

xline=0.3
axes[1,2].set(xlabel='$k/U^2_{ref}$', title='x='+str(xline)+'m')
axes[1,2].plot(*plotdata(df_cfd, xline, var), c='r', ls='--')
axes[1,2].plot(*plotdata(df_exp, xline, var), 'o', mfc='none', ms=4, c='k')

axes[0,0].legend(['CFD','Exp'],loc='upper center', framealpha=1)

for axis1 in axes:
    for axis in axis1:
        axis.tick_params(axis="y",direction="in")
        axis.tick_params(axis="x",direction="in")
        axis.set(ylim=([0,0.25]), xlim=([0,0.2]))

plt.tight_layout()
plt.savefig('plot_k.jpg', dpi=300)

