import uproot
import numpy as np
import pandas as pd
from tqdm import tqdm
from config import *

def do_df(run, location, mc_range = 0, isMC=False, withDWC=True, h3_esum = False):
    branches = [u'rechit_chip', 'rechit_module', 'rechit_channel',
                u'rechit_energy', 'rechit_layer']
    if isMC:
        branches.append(u'ahc_energySum')

    dwc_branches = [u'ntracks', 'dwcReferenceType', 'b_x', 'b_y'] 
    rh_branches = [branch for branch in branches if 'rechit' in branch]
    xcorr = - 3.6 + 2.7 #x MC
    ycorr = + 2.6 - 1.0 #y MC
    
    fname = location + 'ntuple_%i.root' %run
    if isMC:
        if 'CMSSW11_0_withAHCAL_newBeamline' in location:
            fname = location + 'ntuple_sim_config22_pdgID11_beamMomentum%i_listFTFP_BERT_EMN_0000_%i.root' %(run, mc_range)
        else:
            fname = location + 'ntuple_sim_config22_pdgID11_beamMomentum%i_listFTFP_BERT_EMN_0000_0.root' %run
    key='rechitntupler/hits'
    df = uproot.open(fname)[key].pandas.df(branches)
    df = df.reset_index(level=1,drop=True)
    df.index.names = ['event']

    if withDWC:
        key_dwc='trackimpactntupler/impactPoints'
        dwc = uproot.open(fname)[key_dwc].pandas.df(dwc_branches)
        dwc = dwc.reset_index(drop=True)
        dwc.index.names = ['event']
        df_dwc = df.join(dwc)    
        sel_dwc = (df_dwc.ntracks == 1) & (df_dwc.dwcReferenceType == 13)

    else:
        df_dwc = df.copy()

    sel = (df_dwc.rechit_energy > 0.5)
    sel &= ~((df_dwc.rechit_chip == 3) & (df_dwc.rechit_channel == 22))
    sel &= ~((df_dwc.rechit_module == 78) & (df_dwc.rechit_chip == 0))
    sel &= (((df_dwc.rechit_layer != 37) & (df_dwc.rechit_layer != 36)))
    
    if withDWC:
        sel &= sel_dwc
        
    if h3_esum:
        return df_dwc[sel]
    
    df_sel = df_dwc[sel].copy()
    
    for name in df_sel.columns:
        df_sel[name] = df_sel[name].astype('float16')

    if isMC:
        if 'CMSSW11_0_withAHCAL_newBeamline' in location:
            df_sel['b_x'] = (df_sel['b_x'] - np.median(df_sel.b_x))
            df_sel['b_y'] = (df_sel['b_y'] - np.median(df_sel.b_y))
        else:
            df_sel['b_x'] = -df_sel['b_x'] + xcorr
            df_sel['b_y'] = -df_sel['b_y'] + ycorr        
    else:
        df_sel['b_x'] = df_sel['b_x'] + 2.7
        df_sel['b_y'] = df_sel['b_y'] - 1.
    
    
    sel = abs(df_sel.b_x) < 1.0
    sel &= abs(df_sel.b_y) < 1.0
    sel &= (df_sel.rechit_layer < 29)
    if isMC:
        sel &= (df_sel.ahc_energySum == 0)

    sel &= (df_sel[(df_sel.rechit_layer>28) & (df_sel.rechit_layer<=40)].groupby('event').rechit_layer.size() < 50)

    esum_tot = df_sel.groupby('event').rechit_energy.sum()
    esum_EE = df_sel[df_sel.rechit_layer < 29].groupby('event').rechit_energy.sum()
    cut = esum_EE/esum_tot > 0.95

    df_sel = df_sel[sel & cut].copy()

    return df_sel

def prepareHDF(energy, isInj=False, isMC=False, outdir='./'):
    dfs = []
    if isMC:
        runs_id = range(5)
        prefix = sim_prefix
        hdf_name = 'mc_%i.h5' %energy
    else:
        runs_id = runs[energy]
        prefix = data_prefix
        hdf_name = 'data_%i.h5' %energy
        if isInj:
            prefix = inj_prefix
            hdf_name = 'data_inj_%i.h5' %energy

    fout = outdir + hdf_name

    for run in tqdm(runs_id):
        try:
            if isMC:
                df_tmp = do_df(energy, prefix, mc_range = run, isMC = isMC)
            else:
                df_tmp = do_df(run, prefix)
            dfs.append(df_tmp)
        except:
            continue
    df = pd.concat(dfs, keys = range(len(dfs)))
    df.index = ['{}_{}'.format(i, j) for i, j in df.index]
    df.index.names = ['event']
    df.to_hdf(fout, key='df', mode='w')