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

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--inj", type=bool, help="", default=False, required=False)

args = parser.parse_args()
isInj = args.inj

plt.style.use(url)

energies = np.array([20, 30, 50, 80, 100, 120, 150, 200, 250, 300])
true_E   = np.array([20, 30, 49.99, 79.93, 99.83, 119.20, 149.14, 197.32, 243.61, 287.18])

obs_dir = '../data/'

reso_dt = []; reso_mc = []; reso_inj = []
err_reso_dt = []; err_reso_mc = []; err_reso_inj = []
mip_dt  = []; mip_mc  = []; mip_inj = []
err_mip_dt = []; err_mip_mc = []; err_mip_inj = []

for energy, tE in zip(tqdm(energies), true_E):
    fname = obs_dir + 'esum_data_%i.h5' %energy
    esum_data = pd.read_hdf(fname, 's')
    res = repeatedGausFit(esum_data, tE)
    reso_dt.append(res[0]); mip_dt.append(res[2])
    err_reso_dt.append(res[1]); err_mip_dt.append(res[3])

    if isInj:
        fname = obs_dir + 'esum_data_inj_%i.h5' %energy
        esum_inj = pd.read_hdf(fname, 's')
        res = repeatedGausFit(esum_inj, tE)
        reso_inj.append(res[0]); mip_inj.append(res[2])
        err_reso_inj.append(res[1]); err_mip_inj.append(res[3])

    fname = obs_dir + 'esum_mc_%i.h5' %energy
    esum_mc = pd.read_hdf(fname, 's')
    res = repeatedGausFit(esum_mc, tE)
    reso_mc.append(res[0]); mip_mc.append(res[2])
    err_reso_mc.append(res[1]); err_mip_mc.append(res[3])

reso_dt = np.array(reso_dt); reso_mc = np.array(reso_mc); reso_inj = np.array(reso_inj)
err_reso_dt = np.array(err_reso_dt); err_reso_mc = np.array(err_reso_mc); err_reso_inj = np.array(err_reso_inj)
mip_dt  = np.array(mip_dt); mip_mc  = np.array(mip_mc); mip_inj = np.array(mip_inj)
err_mip_dt = np.array(err_mip_dt); err_mip_mc = np.array(err_mip_mc); err_mip_inj = np.array(err_mip_inj)

##############
## Plotting ##
##############

plots_dir = '../plots/'
isdir = os.path.isdir(plots_dir) 
if not isdir: 
    print('Directory {} does not exist. Creating it.' .format(plots_dir))  
    os.mkdir(plots_dir)

#### EM Resolution ####

plt.figure(figsize = (6, 4))

### Data
p_dt, l_dt = plotReso(reso_dt, err_reso_dt)

### Injection data
if isInj:
    p_inj, l_inj = plotReso(reso_inj, err_reso_inj, isInj = isInj)

### Simulation
p_mc, l_mc = plotReso(reso_mc, err_reso_mc, isMC = True)

if isInj:
    plt.legend(handles=[p_mc, l_mc, p_dt, l_dt, p_inj, l_inj], fontsize = 8)
else:
    plt.legend(handles=[p_mc, l_mc, p_dt, l_dt], fontsize = 8)

plt.show()
plt.savefig(plots_dir + 'resolution.pdf', bbox_inches='tight')

#### EM Linearity ####
plt.figure(figsize = (6, 4))
### Data
p_dt, l_dt = plotLinearity(mip_dt, err_mip_dt)
### Simulation
p_mc, l_mc = plotLinearity(mip_mc, err_mip_mc, isMC = True)
plt.legend(handles=[p_mc, l_mc, p_dt, l_dt], fontsize = 8)

plt.show()
plt.savefig(plots_dir + 'linearity.pdf', bbox_inches='tight')

#### EM Linearity ratio ####
plotLinearitySF(mip_dt, err_mip_dt, mip_mc, err_mip_mc)
plt.show()
plt.savefig(plots_dir + 'linearitySF.pdf', bbox_inches='tight')

