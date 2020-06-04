import uproot
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

url = 'https://gist.githubusercontent.com/bonanomi/d14780f7562cb2a22fdd753a9d4459d4/raw/034716767493fcfb7852c0c0e4555b86eafbb901/MyMPLStyle'

plt.style.use(url)

energies = np.array([20, 30, 50, 80, 100, 120, 150, 200, 250, 300])
true_E   = np.array([20, 30, 49.99, 79.93, 99.83, 119.20, 149.14, 197.32, 243.61, 287.18])

obs_dir = '../data/'

cg_dt = []
cg_mc  = []

for energy in tqdm(energies):
	fname = obs_dir + 'cog_data_%i.h5' %energy
    cog_data = pd.read_hdf(fname, 's')
    cg_dt.append(np.median(cog_data))

    fname = obs_dir + 'cog_mc_%i.h5' %energy
    cog_mc = pd.read_hdf(fname, 's')
    cg_mc.append(np.median(cog_mc))


##############
## Plotting ##
##############

plt.figure(figsize = (6, 4))
plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)
plt.plot(true_E, cg_dt, 'k.', label = 'Data')

plt.plot(true_E, cg_mc, marker = 's', mfc = 'None',
             color = 'r', linestyle = 'None', label = 'MC')

plt.grid(b = None)
plt.legend(fontsize = 8)
plt.ylabel(r'$\left\langle COG_z\right\rangle$', ha='right', y=1.0, fontsize = 12)
plt.xlabel('E beam [GeV]', ha='right', x=1.0, fontsize = 12)
plt.show()
plt.savefig('cog_average.pdf', bbox_inches='tight')
