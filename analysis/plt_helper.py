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

url = 'https://gist.githubusercontent.com/bonanomi/d14780f7562cb2a22fdd753a9d4459d4/raw/c8fd6811858458ebb3e0a34f5a3a1a9584bcd7ce/MyMPLStyle'

plt.style.use(url)

energies = np.array([20, 30, 50, 80, 100, 120, 150, 200, 250, 300])
true_E   = np.array([20, 30, 49.99, 79.93, 99.83, 119.20, 149.14, 197.32, 243.61, 287.18])
draw_e = np.linspace(15, 300, 1000)

obs_dir = '../data/'

def plotReso(res, err, isMC = False, isInj = False):
    plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)

    label = 'Data'; color = 'k'; marker = '.'
    if isMC:
        label = 'MC'; color = 'r'; marker = 's'
    if isInj:
        label = 'Data (injection)'; color = 'b'
    plt_args = {'marker': marker, 'color': color, 'linestyle':'None', 'label': label}
    if isMC: plt_args['mfc'] = 'None'

    plt.errorbar(true_E, res, yerr = err, **plt_args)
    popt, pcov = curve_fit(resolution, true_E, res, sigma =err, absolute_sigma = True)
    stoc = popt[0]*1e2; const = popt[1]*1e2
    err_stoc = np.sqrt(np.diag(pcov))[0]*1e2; err_const = np.sqrt(np.diag(pcov))[1]*1e2
    plt.plot(draw_e, resolution(draw_e, *popt), '--', color = color)

    dot_args = {'color':color, 'marker':marker, 'markersize':3.5, 'label':label}
    if isMC:
        dot_args['mfc'] = 'None'
    dot = mlines.Line2D([], [], **dot_args)

    line = mlines.Line2D([], [], color=color, marker='None', mfc = 'None', linestyle = '-.',
                          markersize=3.5,
             label = r'S = %.2f $\pm$ %.2f $\sqrt{GeV}$ %%' '\n' r'C = %.2f $\pm$ %.2f %%' %(stoc, err_stoc, const, err_const))
    
    plt.grid(b=None)
    plt.ylabel(r'Gaussian $\frac{\sigma_{E}}{\left\langle E\right\rangle}$', ha='right', y=1.0, fontsize = 12)
    plt.xlabel('E beam [GeV]', ha='right', x=1.0, fontsize = 12)

    return dot, line

def plotLinearity(mip, err, isMC = False, isInj = False):
    plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)

    label = 'Data'; color = 'k'; marker = '.'
    if isMC:
        label = 'MC'; color = 'r'; marker = 's'
    if isInj:
        label = 'Data (injection)'; color = 'b'
    plt_args = {'marker': marker, 'color': color, 'linestyle':'None', 'label': label}
    if isMC: plt_args['mfc'] = 'None'

    plt.errorbar(true_E, mip, yerr = err, **plt_args)
    popt, pcov = curve_fit(fit_line, true_E, mip, sigma =err, absolute_sigma = True)
    stoc = popt[0]*1e2; const = popt[1]*1e2
    err_stoc = np.sqrt(np.diag(pcov))[0]*1e2; err_const = np.sqrt(np.diag(pcov))[1]*1e2
    plt.plot(draw_e, fit_line(draw_e, *popt), '--', color = color)

    dot_args = {'color':color, 'marker':marker, 'markersize':3.5, 'label':label}
    if isMC:
        dot_args['mfc'] = 'None'
    dot = mlines.Line2D([], [], **dot_args)

    line = mlines.Line2D([], [], color=color, marker='None', mfc = 'None', linestyle = '-.',
                          markersize=3.5,
             label = r'E = (%.2f $\pm$ %.2f $\frac{MIP}{GeV}$) + (%.2f $\pm$ %.2f)' %(stoc, err_stoc, const, err_const))
    
    plt.grid(b=None)
    plt.ylabel(r'$\left\langle E\right\rangle$', ha='right', y=1.0, fontsize = 12)
    plt.xlabel('E beam [GeV]', ha='right', x=1.0, fontsize = 12)

    return dot, line


def plotShowerMax(tmax, err, isMC = False, isInj = False):
    plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)

    label = 'Data'; color = 'k'; marker = '.'
    if isMC:
        label = 'MC'; color = 'r'; marker = 's'
    if isInj:
        label = 'Data (injection)'; color = 'b'
    plt_args = {'marker': marker, 'color': color, 'linestyle':'None', 'label': label}
    if isMC: plt_args['mfc'] = 'None'

    plt.errorbar(true_E, tmax, yerr = err, **plt_args)
    popt, pcov = curve_fit(logE, true_E, tmax, sigma =err, absolute_sigma = True)
    Ec = popt*1e3; err_Ec = np.sqrt(np.diag(pcov*1e3))

    plt.plot(draw_e, logE(draw_e, *popt), '--', color = color)

    dot_args = {'color':color, 'marker':marker, 'markersize':3.5, 'label':label}
    if isMC:
        dot_args['mfc'] = 'None'
    dot = mlines.Line2D([], [], **dot_args)

    line = mlines.Line2D([], [], color=color, marker='None', mfc = 'None', linestyle = '-.',
                          markersize=3.5,
                  label = r'$\log\left(\frac{E}{%.2f \pm %.2f \, MeV}\right) - 0.5$' %(Ec, err_Ec))
    
    plt.grid(b=None)
    plt.ylabel(r'Shower max [X0]', ha='right', y=1.0, fontsize = 12)
    plt.xlabel('E beam [GeV]', ha='right', x=1.0, fontsize = 12)

    return dot, line

def plotCOGFit(cog, err, isMC = False, isInj = False):
    plt.title(r'$\bf{{CMS}}$' ' Preliminary', style = 'italic', loc = 'left', fontsize = 10)

    label = 'Data'; color = 'k'; marker = '.'
    if isMC:
        label = 'MC'; color = 'r'; marker = 's'
    if isInj:
        label = 'Data (injection)'; color = 'b'
    plt_args = {'marker': marker, 'color': color, 'linestyle':'None', 'label': label}
    if isMC: plt_args['mfc'] = 'None'

    plt.errorbar(true_E, cog, yerr = err, **plt_args)
    popt, pcov = curve_fit(logE_cog, true_E, cog, sigma =err, absolute_sigma = True)
    Ec = popt*1e3; err_Ec = np.sqrt(np.diag(pcov*1e3))

    plt.plot(draw_e, logE_cog(draw_e, *popt), '--', color = color)

    dot_args = {'color':color, 'marker':marker, 'markersize':3.5, 'label':label}
    if isMC:
        dot_args['mfc'] = 'None'
    dot = mlines.Line2D([], [], **dot_args)

    line = mlines.Line2D([], [], color=color, marker='None', mfc = 'None', linestyle = '-.',
                          markersize=3.5,
                  label = r'1.5 + $\log\left(\frac{E}{%.2f \pm %.2f \, MeV}\right)$' %(Ec, err_Ec))
    
    plt.grid(b=None)
    plt.ylabel(r'$\left\langle COG_z\right\rangle$ [X0]', ha='right', y=1.0, fontsize = 12)
    plt.xlabel('E beam [GeV]', ha='right', x=1.0, fontsize = 12)

    return dot, line