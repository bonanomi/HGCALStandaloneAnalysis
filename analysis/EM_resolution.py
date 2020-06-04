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

reso_dt = []; reso_mc = []
mip_dt  = []; mip_mc  = []

for energy, tE in zip(tqdm(energies), true_E):
    fname = obs_dir + 'esum_data_%i.h5' %energy
    esum_data = pd.read_hdf(fname, 's')
    res = repeatedGausFit(esum_data, tE)
    reso_dt.append(res[0]); mip_dt.append(res[2])
    
    fname = obs_dir + 'esum_mc_%i.h5' %energy
    esum_mc = pd.read_hdf(obs_dir, 's')
    res = repeatedGausFit(esum_mc, tE)
    reso_mc.append(res[0]); mip_mc.append(res[2])


##############
## Plotting ##
##############

draw_e = np.linspace(15, 300, 1000)

plt.figure(figsize = (6, 4))
plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)
plt.plot(true_E, reso_dt, 'k.', label = 'Data')
popt, pcov = curve_fit(resolution, true_E, reso_dt)
sc_dt = popt[0]*1e2; cons_dt = popt[1]*1e2
err_sc_dt = np.sqrt(np.diag(pcov))[0]*1e2; err_const_dt = np.sqrt(np.diag(pcov))[1]*1e2
plt.plot(draw_e, resolution(draw_e, *popt), '--', color = 'k')

plt.plot(true_E, reso_mc, marker = 's', mfc = 'None', 
             color = 'r', linestyle = 'None', label = 'MC')
popt, pcov = curve_fit(resolution, true_E, reso_mc)
sc_mc = popt[0]*1e2; cons_mc = popt[1]*1e2
err_sc_mc = np.sqrt(np.diag(pcov))[0]*1e2; err_const_mc = np.sqrt(np.diag(pcov))[1]*1e2
plt.plot(draw_e, resolution(draw_e, *popt), '--', color = 'r')

k_dot = mlines.Line2D([], [], color='k', marker='o',alpha = 0.8, 
                          markersize=3.5, label='Data')
k_line = mlines.Line2D([], [], color='k', marker='None', mfc = 'None', linestyle = '-.',
                          markersize=3.5, 
label = r'S = %.2f $\pm$ %.2f $\sqrt{GeV}$ %%' '\n' r'C = %.2f $\pm$ %.2f %%' %(sc_dt, err_sc_dt, cons_dt, err_const_dt))

red_dot = mlines.Line2D([], [], color='red', marker='o', mfc = 'None',
                          markersize=3.5, label= 'MC v3')
red_line = mlines.Line2D([], [], color='red', marker='None', mfc = 'None', linestyle = '-.',
                          markersize=3.5, 
label = r'S = %.2f $\pm$ %.2f $\sqrt{GeV}$ %%' '\n' r'C = %.2f $\pm$ %.2f %%' %(sc_mc, err_sc_mc, cons_mc, err_const_mc))

plt.legend(handles=[red_dot, red_line, k_dot, k_line], fontsize = 8)
plt.ylabel(r'Gaussian $\frac{\sigma_{E}}{\left\langle E\right\rangle}$', ha='right', y=1.0, fontsize = 12)
plt.xlabel('E beam [GeV]', ha='right', x=1.0, fontsize = 12)
plt.show()
plt.savefig('resolution.pdf', bbox_inches='tight')
