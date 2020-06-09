import uproot
import os.path
import numpy as np
import pandas as pd
from tqdm import tqdm

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches

from scipy.optimize import curve_fit

from fits import *
from plt_helper import *

plt.style.use(url)

energies = np.array([20, 30, 50, 80, 100, 120, 150, 200, 250, 300])
true_E   = np.array([20, 30, 49.99, 79.93, 99.83, 119.20, 149.14, 197.32, 243.61, 287.18])

obs_dir = '../data/'

tmax_dt = []; err_tmax_dt = []
tmax_mc = []; err_tmax_mc = []

cog_dt = []; err_cog_dt = []
cog_mc = []; err_cog_mc = []

for energy in tqdm(energies):
    fname = obs_dir + 'lp_data_%i.h5' %energy
    lp_data = pd.read_hdf(fname, 's')
    x = np.array(lp_data.index); y = lp_data.values
    tmax, er_tmax = showerMaxFit(x,y)
    tmax_dt.append(tmax); err_tmax_dt.append(er_tmax)

    cog, er_cog = cogFit(x,y)
    cog_dt.append(cog); err_cog_dt.append(er_cog)

    fname = obs_dir + 'lp_mc_%i.h5' %energy
    lp_mc = pd.read_hdf(fname, 's')
    x = np.array(lp_mc.index); y = lp_mc.values
    tmax, er_tmax = showerMaxFit(x, y)
    tmax_mc.append(tmax); err_tmax_mc.append(er_tmax)

    cog, er_cog = cogFit(x, y)
    cog_mc.append(cog); err_cog_mc.append(er_cog)

tmax_dt = np.array(tmax_dt); err_tmax_dt = np.array(err_tmax_dt)
tmax_mc = np.array(tmax_mc); err_tmax_mc = np.array(err_tmax_mc)
cog_dt = np.array(cog_dt); err_cog_dt = np.array(err_cog_dt)
cog_mc = np.array(cog_mc); err_cog_mc = np.array(err_cog_mc)

##############
## Plotting ##
##############

plots_dir = '../plots/'
isdir = os.path.isdir(plots_dir) 
if not isdir: 
  print('Directory {} does not exist. Creating it.' .format(plots_dir))  
    os.mkdir(plots_dir)

### Shower max from fit

plt.figure(figsize = (6, 4))
### Data
p_dt, l_dt = plotShowerMax(tmax_dt, err_tmax_dt)
### Simulation
p_mc, l_mc = plotShowerMax(tmax_mc, err_tmax_mc, isMC = True)
plt.legend(handles=[p_mc, l_mc, p_dt, l_dt], fontsize = 8)

plt.show()
plt.savefig(plots_dir + 'tmax.pdf', bbox_inches='tight')


### Average COGz from fit

plt.figure(figsize = (6, 4))
### Data
p_dt, l_dt = plotCOGFit(cog_dt, err_cog_dt)
### Simulation
p_mc, l_mc = plotCOGFit(cog_mc, err_cog_mc, isMC = True)
plt.legend(handles=[p_mc, l_mc, p_dt, l_dt], fontsize = 8)

plt.show()
plt.savefig(plots_dir + 'cogz_from_fit.pdf', bbox_inches='tight')
