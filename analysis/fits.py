import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf, erfc, gamma

def gaussian(x, A, mu, sigma):
    '''A scaled normal distribution.
    '''
    return A*np.exp(-0.5*((x - mu)/sigma)**2)

def crystal_ball(x, a, n, xb, sig):
    '''Crystal ball function with a left tail.
    '''
    x = x+0j
    if a < 0:
        a = -a
    if n < 0:
        n = -n
    aa = abs(a)
    A = (n/aa)**n * np.exp(- aa**2 / 2)
    B = n/aa - aa
    C = n/aa*1/(n-1)*np.exp(- aa**2 / 2)
    D = np.sqrt(np.pi/2)*(1+erf(aa/np.sqrt(2)))
    N = 1./(sig*(C+D))
    total += ((x-xb)/sig  > -a) * N * np.exp(- (x-xb)**2/(2.*sig**2))
    total += ((x-xb)/sig <= -a) * N * A * (B - (x-xb)/sig)**(-n)

    try:
        return total.real
    except:
        return total
    return total

def err_prop(sigma, E, err_sigma, err_E):
    '''Propagation of error for sigma/E.
    '''
    N = sigma/E
    dval1 = err_sigma/sigma
    dval2 = err_E/E
    return N*np.sqrt((dval1)**2 + (dval2)**2)

def resolution(x, stoc, const):
    '''Energy resolution function, w/o noise term.
    '''
    return np.sqrt((stoc / np.sqrt(x))**2 + const**2)#+ (noise / x) ** 2)

def showerShape(t, alpha, beta, E):
    '''Gamma distribution to fit longitudial shower shapes in particles calorimeters.
    See equation 2 in https://arxiv.org/pdf/hep-ex/0001020v1.pdf
    '''
    return E * ((beta * t) ** (alpha - 1) * beta * np.exp(-beta * t)) / gamma(alpha)

def Longo_fcn(z, a, w, b):
    '''Longo parametrization to fit longitudial shower shapes in particles calorimeters.
    See http://inspirehep.net/record/105195/?ln=fr
    '''
    return a*z**w*np.exp(-b*z)

def line(x, m, q):
    return m*x + q

def repeatedGausFit(histogram, energy, weights = [], drawFit=False, c = 'k', rangeInSigmaLeft = 1., rangeInSigmaRight = 2.5):
    Sampl = 2.198e-01
    Noise = 8.210e-03
    Linear = 3e-6
    h_med = histogram.median()

    '''Check whether we're working with raw energies or visible energy in GeV.
        Binning and std_dev of the initial distribution are set accrodingly.
    '''
    if h_med > 1000:
        h_std = np.sqrt(Noise*Noise + Sampl*Sampl*energy + Linear*Linear*(energy**2))*100.0
        binning = np.linspace(np.median(histogram)*0.8, np.median(histogram)*1.2, 100)
    if h_med < 1000:
        h_std = np.sqrt(Noise*Noise + Sampl*Sampl*energy + Linear*Linear*(energy**2))
        binning = np.linspace(np.median(histogram)*0.8, np.median(histogram)*1.2, 80)

    n, bins = np.histogram(histogram, density = 1, bins=binning)
    if len(weights) > 1:
        n, bins = np.histogram(histogram, density = 1, bins=binning, weights = weights)
    bins = 0.5*(bins[1:]+bins[:-1])
    ymaximum = n.max()
    s_fit = (bins > h_med - rangeInSigmaLeft*h_std)& (bins < (h_med + rangeInSigmaRight*h_std))
    popt, pcov = curve_fit(gaussian, bins[s_fit], n[s_fit], p0=[ymaximum, h_med, h_std])

    for i in range(9):
        mu_tmp = popt[1]
        std_tmp = popt[2]
        ymaximum_tmp = popt[0]
        s_fit = (bins > mu_tmp- rangeInSigmaLeft*std_tmp) & (bins < mu_tmp + rangeInSigmaRight*std_tmp)

        popt, pcov = curve_fit(gaussian, bins[s_fit], n[s_fit], p0=[ymaximum_tmp, mu_tmp, std_tmp])

    e_reco = popt[1]
    sigma = popt[2]
    err = np.sqrt(np.diag(pcov))
    E_err = err[1]
    sigma_err = err[2]
    res = ((n[s_fit] - gaussian(bins[s_fit], *popt))/gaussian(bins[s_fit], *popt))**2
    chi2 = np.sum(res)/(len(bins[s_fit])-len(popt))
    resolution = sigma/e_reco
    error = err_prop(sigma, e_reco, sigma_err, E_err)

    if drawFit==True:
     xdraw = np.linspace(popt[1]-rangeInSigmaLeft*popt[2],popt[1]+rangeInSigmaRight*popt[2], 1000)
     plt.plot(xdraw, gaussian(xdraw, popt[0], popt[1],  popt[2]), '--', linewidth=2.0, color = c)

    return resolution, error, e_reco, E_err, sigma, chi2
