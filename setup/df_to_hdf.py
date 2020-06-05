import os.path
import numpy as np
import uproot
from tqdm import tqdm
import pandas as pd
from reader import prepareHDF
from config import *

isdir = os.path.isdir(path_hdf) 
if not isdir: 
    print('Directory {} does not exist. Creating it.' .format(path_hdf))  
    os.mkdir(path_hdf)

for energy in energies:
	prepareHDF(energy, outdir = path_hdf)
    prepareHDF(energy, inInj = True, outdir = path_hdf)
    prepareHDF(energy, isMC = True, outdir = path_hdf)