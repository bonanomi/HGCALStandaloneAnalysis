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

url = 'https://gist.githubusercontent.com/bonanomi/d14780f7562cb2a22fdd753a9d4459d4/raw/034716767493fcfb7852c0c0e4555b86eafbb901/MyMPLStyle'

plt.style.use(url)

energies = np.array([20, 30, 50, 80, 100, 120, 150, 200, 250, 300])
true_E   = np.array([20, 30, 49.99, 79.93, 99.83, 119.20, 149.14, 197.32, 243.61, 287.18])

obs_dir = '../data/'

reso_dt = []; reso_mc = []
mip_dt  = []; mip_mc  = []

for energy, tE in zip(tqdm(energies), true_E):
    fname = obs_dir + 'cog_data_%i.h5' %energy
    cog_data = pd.read_hdf(fname, 's')

    fname = obs_dir + 'cog_mc_%i.h5' %energy
    cog_mc = pd.read_hdf(fname, 's')

    binning = np.linspace(6, 15, 50)

    t = plt.hist(cog_data, binning, density = True, histtype = 'step', color = 'k', label = 'Data') 
    t = plt.hist(cog_mc, binning, density = True, histtype = 'step', color = 'r', label = 'MC') 

    plt.figure(figsize = (6, 4))
    plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)

	plt.grid(b = None)
	plt.legend(fontsize = 10, title = r'$e^+$ %i GeV' %energy, title_fontsize = 12)
	plt.ylabel(r'a.u.', ha='right', y=1.0, fontsize = 12)
	plt.xlabel('Shower depth [X0]', ha='right', x=1.0, fontsize = 12)
	plt.show()
	plt.savefig('EM_cog_%i.pdf' %energy, bbox_inches='tight')
