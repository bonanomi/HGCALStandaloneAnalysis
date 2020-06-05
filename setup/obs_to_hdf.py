import uproot
import os.path
import numpy as np
import pandas as pd
from tqdm import tqdm
from config import *
from hdf_helpers import esum_to_hdf, cog_to_hdf, lp_to_hdf

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--obs", type=str, help="esum, cog, longprof", default='esum', required=False)

args = parser.parse_args()
observable = args.obs

outdir = '../data/'

isdir = os.path.isdir(outdir) 
if not isdir: 
	print('Directory {} does not exist. Creating it.' .format(outdir))  
    os.mkdir(outdir)

for energy in tqdm(energies):
    if observable == 'esum':
        esum_to_hdf(energy, outdir = outdir)
        esum_to_hdf(energy, isInj = True, outdir = outdir)
        esum_to_hdf(energy, isMC=True, outdir = outdir)
    elif observable == 'cog':
        cog_to_hdf(energy, outdir = outdir)
        cog_to_hdf(energy, isMC=True, outdir = outdir)
    elif observable == 'longprof':
        lp_to_hdf(energy, outdir = outdir)
        lp_to_hdf(energy, isMC=True, outdir = outdir)
    else:
        print('Observable not defined.')
        break