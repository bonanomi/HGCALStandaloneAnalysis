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

url = 'https://gist.githubusercontent.com/bonanomi/d14780f7562cb2a22fdd753a9d4459d4/raw/c8fd6811858458ebb3e0a34f5a3a1a9584bcd7ce/MyMPLStyle'

plt.style.use(url)

energies = np.array([20, 30, 50, 80, 100, 120, 150, 200, 250, 300])
true_E   = np.array([20, 30, 49.99, 79.93, 99.83, 119.20, 149.14, 197.32, 243.61, 287.18])

obs_dir = '../data/'
plots_dir = '../plots/'
isdir = os.path.isdir(plots_dir) 
if not isdir: 
  print('Directory {} does not exist. Creating it.' .format(plots_dir))  
    os.mkdir(plots_dir)

tmax_dt = []; err_tmax_dt = []
tmax_mc = []; err_tmax_mc = []

cog_dt = []; err_cog_dt = []
cog_mc = []; err_cog_mc = []

three_energies = [20, 100, 300]
dt_longprofs = {}
mc_longprofs = {}

for energy in tqdm(energies):
    plt.figure(figsize = (6, 4))
    plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)

    fname = obs_dir + 'lp_data_%i.h5' %energy
    lp_data = pd.read_hdf(fname, 's')
    x = np.array(lp_data.index); y = lp_data.values
    
    if energy in three_energies:
        dt_longprofs[energy] = [x, y]

    plt.plot(x, y, 'k.--', label = 'Data')

    fname = obs_dir + 'lp_mc_%i.h5' %energy
    lp_mc = pd.read_hdf(fname, 's')
    x = np.array(lp_mc.index); y = lp_mc.values
    if energy in three_energies:
        mc_longprofs[energy] = [x, y]

    plt.plot(x, y, marker = 's', color = 'r', mfc = 'None', linestyle = '--', label = 'MC')

    plt.grid(b = None)
    plt.legend(fontsize = 10, title = r'$e^+$ %i GeV' %energy, title_fontsize = 12)
    plt.ylabel('Energy deposit [MIP]', ha='right', y=1.0, fontsize = 12)
    plt.xlabel('Shower depth [X0]', ha='right', x=1.0, fontsize = 12)
    plt.show()
    plt.savefig(plots_dir + 'EM_longProf_%i.pdf' %energy, bbox_inches='tight')

plt.figure(figsize=(8, 5))
plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc='left', fontsize = 15)

for energy in three_energies:
    x_dt, y_dt = dt_longprofs[energy]
    x_mc, y_mc = mc_longprofs[energy]
    plt.plot(x_dt, y_dt, 'k.--')
    plt.plot(x_mc, y_mc, marker = 's', color = 'r', mfc = 'None', linestyle = '--')

plt.text(11, 250,'20 GeV', style = 'italic')
plt.text(11, 1000,'100 GeV', style = 'italic')
plt.text(11, 2550,'300 GeV', style = 'italic')
plt.ylabel('Energy deposit [MIP]', ha='right', y=1.0, fontsize = 13)
plt.xlabel('Shower depth [X0]', ha='right', x=1.0, fontsize = 13)
k_dot = mlines.Line2D([], [], color='k', marker='o', linestyle = '--', alpha = 0.8, 
                          markersize=3.5, label='Data')
r_dot = mlines.Line2D([], [], color='r', marker='s', linestyle = '--', mfc = 'None', 
                          markersize=3.5, label='MC')
plt.legend(handles=[k_dot, r_dot],fontsize = 12)
plt.grid(b = None)
plt.show()
plt.savefig(plots_dir + 'EM_three_longProf.pdf', bbox_inches='tight')