import numpy as np
import pandas as pd
from observables import *
from config import energies, runs, data_prefix, sim_prefix, path_hdf

def esum_to_hdf(energy, isMC=False, dEdx = False, outdir = './'):
    if isMC:
        hdf_name = 'esum_mc_%i.h5' %energy
        df = pd.read_hdf(path_hdf+'mc_%i.h5' %energy, 'df')
    else:
        hdf_name = 'esum_data_%i.h5' %energy
        df = pd.read_hdf(path_hdf+'data_%i.h5' %energy, 'df')

    fout = outdir + hdf_name

    if dEdx:
    	fout = outdir + 'dEdx_' + hdf_name
        esum = get_totE_GeV(df)
    else:
        esum = get_totE(df)

    esum.to_hdf(fout, key='s')

def cog_to_hdf(energy, isMC=False, outdir = './'):
    if isMC:
        hdf_name = 'cog_mc_%i.h5' %energy
        df = pd.read_hdf(path_hdf+'mc_%i.h5' %energy, 'df')
    else:
        hdf_name = 'cog_data_%i.h5' %energy
        df = pd.read_hdf(path_hdf+'data_%i.h5' %energy, 'df')

    fout = outdir + hdf_name

    cog = do_COG(df)
    cog.to_hdf(fout, key='s')

def lp_to_hdf(energy, isMC=False, dEdx = False, outdir = './'):
    if isMC:
        hdf_name = 'cog_mc_%i.h5' %energy
        df = pd.read_hdf(path_hdf+'mc_%i.h5' %energy, 'df')
    else:
        hdf_name = 'cog_data_%i.h5' %energy
        df = pd.read_hdf(path_hdf+'data_%i.h5' %energy, 'df')

    fout = outdir + hdf_name

    if dEdx:
    	fout = outdir + 'dEdx_' + hdf_name
        lp = get_longProf_dEdx(df)
    else:
        lp = get_longProf(df)

    lp.to_hdf(fout, key='s')