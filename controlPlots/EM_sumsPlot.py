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
    fname = obs_dir + 'esum_data_%i.h5' %energy
    esum_data = pd.read_hdf(fname, 's')

    fname = obs_dir + 'esum_mc_%i.h5' %energy
    esum_mc = pd.read_hdf(obs_dir, 's')

    binning = np.linspace(np.median(esum_mc)*0.4, np.median(esum_mc)*1.1, 100)

    sf = np.median(esum_mc)/np.median(esum_data)

    t = plt.hist(esum_data*sf, binning, density = True, histtype = 'step', color = 'k', label = 'Data') 
    t = plt.hist(esum_mc, binning, density = True, histtype = 'step', color = 'r', label = 'MC') 

    plt.figure(figsize = (6, 4))
    plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)

	plt.grid(b = None)
	plt.legend(fontsize = 10, title = r'$e^+$ %i GeV' %energy, title_fontsize = 12)
	plt.ylabel(r'a.u.', ha='right', y=1.0, fontsize = 12)
	plt.xlabel('Reconstructed energy [MIP]', ha='right', x=1.0, fontsize = 12)
	plt.show()
	plt.savefig('EM_sum_%i.pdf' %energy, bbox_inches='tight')
