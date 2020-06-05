import numpy as np

## Not considering runs with ID < 424, as they are problematic.
## This affects only 250 and 300 GeV statistics.

runs = {
    20:  [436, 437, 439, 441, 442, 443, 444, 447, 450, 451, 452, 453, 455],
    30:  [594, 595, 596, 597, 599, 601, 603, 604, 606, 607],
    50:  [456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 608, 609, 610, 611, 613, 614, 616, 617, 618, 619],
    80:  [466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 655, 656, 657, 659, 661, 663],
    100: [477, 479, 480, 481, 482, 484, 486, 487, 489, 490, 491],
    120: [620, 621, 622, 635, 636, 637, 639, 640, 641, 642, 643, 644],
    150: [493, 494, 495, 496, 501, 502, 503, 504, 505, 506, 507, 508, 509],
    200: [664, 665, 666, 667, 671, 672, 673, 674, 675, 676],
    250: [645, 646, 647, 648, 649, 650, 652, 653, 654], 
    300: [425, 424, 425, 426, 429, 431, 435]
}

data_prefix = '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/'
sim_prefix  = '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v3/electrons/'
inj_prefix  = '/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/injectionCalib_v4_patch0/'
path_hdf = '/eos/user/m/mbonanom/SWAN_projects/HGC-PulseFit/HGCAL TB 2018 Analysis/HDF_Analysis/'

energies = np.array([20, 30, 50, 80, 100, 120, 150, 200, 250, 300])
true_E   = np.array([20, 30, 49.99, 79.93, 99.83, 119.20, 149.14, 197.32, 243.61, 287.18])

w = [10.199, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851,
    9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 9.851, 11.360, 11.360, 11.360, 11.360,
    10.995, 10.995, 11.153, 36.139, 57.496, 57.654, 57.654, 57.654, 56.884, 38.620, 39.390,
    57.654, 58.083, 57.757, 56.574, 37.755]

X0 = [0.943, 1.923, 2.858, 3.838, 4.773, 5.753, 6.688, 7.668, 8.603, 9.583,
 10.518, 11.498, 12.433, 13.413, 14.348, 15.328, 16.263, 17.243, 18.178,
 19.157, 20.092, 21.240, 22.174, 23.322, 24.256, 25.505, 26.440,
 27.725, 30.566, 33.501, 36.437, 39.373, 42.308, 45.159, 46.275, 49.211,
 52.146, 55.276, 58.322, 61.175]

lay_X0 = dict()
weights = dict()
MIPtoGeV = 84.9*1e-6

for i in range(len(w)):
    weights[i+1] = (1+(w[i] / (1000 * MIPtoGeV))) * MIPtoGeV
    lay_X0[i+1] = (X0[i])