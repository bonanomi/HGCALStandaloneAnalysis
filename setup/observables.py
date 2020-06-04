import pandas as pd
import numpy as np
from scipy.stats import binned_statistic
from config import lay_X0, weights

'''
 Some get functions for:
    * E sums
    * barycenters (option 'EE' and 'FH' to get only EE/FH barys)
    * layers barycenter
    * longitudinal profiles
    * Nhits profiles
'''

def get_totE(df, event = 'event'):
    totE = df.groupby(event).rechit_energy.sum()
    return totE

def get_totE_GeV(df, event = 'event'):
    ''' Longitudinal profile using dEdx weights.
        The calculation of dEdx weighted energies per layer is done directly here.
    '''
    df = df.copy()
    df['rechit_energy_w'] = df['rechit_energy']*(df['rechit_layer'].map(weights))
    totE = df.groupby(event).rechit_energy_w.sum()
    return totE

def get_baryX(df, option = 'calo', event = 'event'):
    df['x_timesE'] = df.rechit_energy * df.rechit_x
    totE = get_totE(df, event)
    bary_x = df.groupby(event)['x_timesE'].sum()/totE
    if (option == 'EE'):
        sel = df.rechit_layer < 29
        df_sel = df[sel].copy()
        df_sel['x_timesE'] = df_sel.rechit_energy * df_sel.rechit_x
        totE = get_totE(df_sel, event)
        bary_x = df_sel.groupby(event)['x_timesE'].sum()/totE
    if (option == 'FH'):
        sel = df.rechit_layer > 28
        sel &= df.rechit_layer < 40
        df_sel = df[sel].copy()
        df_sel['x_timesE'] = df_sel.rechit_energy * df_sel.rechit_x
        totE = get_totE(df_sel, event)
        bary_x = df_sel.groupby(event)['x_timesE'].sum()/totE
    return bary_x

def get_baryY(df, option = 'calo', event = 'event'):
    df['y_timesE'] = df.rechit_energy * df.rechit_y
    totE = get_totE(df, event)
    bary_y = df.groupby(event)['y_timesE'].sum()/totE
    if option == 'EE':
        sel = df.rechit_layer < 29
        df_sel = df[sel].copy()
        df_sel['y_timesE'] = df_sel.rechit_energy * df_sel.rechit_y
        totE = get_totE(df_sel)
        bary_y = df_sel.groupby(event)['y_timesE'].sum()/totE
    if (option == 'FH'):
        sel = df.rechit_layer > 28
        sel &= df.rechit_layer < 40
        df_sel = df[sel].copy()
        df_sel['y_timesE'] = df_sel.rechit_energy * df_sel.rechit_y
        totE = get_totE(df_sel)
        bary_y = df_sel.groupby(event)['y_timesE'].sum()/totE
    return bary_y

def get_baryZ(df, option = 'calo', event = 'event'):
    df['z_timesE'] = df.rechit_energy * df.rechit_z
    totE = get_totE(df, event)
    bary_z = df.groupby(event)['z_timesE'].sum()/totE
    if option == 'EE':
        sel = df.rechit_layer < 29
        df_sel = df[sel].copy()
        df_sel['z_timesE'] = df_sel.rechit_energy * df_sel.rechit_z
        totE = get_totE(df_sel)
        bary_z = df_sel.groupby(event)['z_timesE'].sum()/totE
    if (option == 'FH'):
        sel = df.rechit_layer > 28
        sel &= df.rechit_layer < 40
        df_sel = df[sel].copy()
        df_sel['z_timesE'] = df_sel.rechit_energy * df_sel.rechit_z
        totE = get_totE(df_sel)
        bary_z = df_sel.groupby(event)['z_timesE'].sum()/totE
    return bary_z

def get_dr(df, option = 'calo', event = 'event'):
    b_x = get_baryX(df, event)
    b_y = get_baryY(df, event)
    dr = np.hypot(b_x, b_y) 
    if option == 'EE':
        b_x = get_baryX(df, 'EE')
        b_y = get_baryY(df, 'EE')
        dr = np.hypot(b_x, b_y)
    return dr 

def get_baryLay(df, event = 'event'):
    df['lay_timesE'] = df.rechit_energy * df.rechit_layer
    totE = get_totE(df, event)
    bary_lay = df.groupby(event)['lay_timesE'].sum()/totE
    return bary_lay

def get_drProf(df, option = 'calo', event = 'event'):
    ''' Radial spread profile per layer (in EE only if option in set to 'EE').
        It returns layer and dr (= rechit_r - bary_r) per layer.
        These two quantities are evaluated as the median over all the events.
    '''
    if option == 'EE':
        sel = df.rechit_layer < 29
        df = df[sel].copy()
        
    bary_x = get_baryX(df, option, event)
    bary_y = get_baryY(df, option, event)
    dx = df['rechit_x'] - bary_x
    dy = df['rechit_y'] - bary_y
    dr = np.hypot(dx,dy)
    x = df.rechit_layer
    y = dr
    x_med, y_med = binned_statistic(x, [x,y], bins = 40, statistic='median').statistic
    sel = ~np.isnan(x_med)
    y_med = y_med[sel]
    x_med = x_med[sel]

    return x_med, y_med

def do_Hits(df, event = 'event'):
    sel = df.rechit_energy > 5
    df_sel = df[sel].copy()
    df_nhits = df_sel.groupby([event,'rechit_layer']).size().reset_index()
    df_nhits.columns = [event,'layer','nhits']
    lay_med_sum = df_nhits.groupby('layer')['nhits'].mean()
    
    return lay_med_sum

def get_longProf(df, event = 'event'):
    #    tot_sum = get_totE(df, event)
    df = df.copy()
    df['X0'] = df['rechit_layer'].map(lay_X0)
    lay_sum = df.groupby([event,'X0'])['rechit_energy'].sum()
    lay_sum_avg = lay_sum.reset_index().groupby('X0')['rechit_energy'].median()
    #    lay_sum_avg = (lay_sum/tot_sum).reset_index().groupby('rechit_layer')['rechit_energy'].mean()
    return lay_sum_avg

def get_longProf_dEdx(df, event = 'event'):
    ''' Longitudinal profile using dEdx weights.
        The calculation of dEdx weighted energies per layer is done directly here.
    '''
    df = df.copy()
    df['rechit_energy_w'] = df['rechit_energy']*(df['rechit_layer'].map(weights))
    df['X0'] = df['rechit_layer'].map(lay_X0)
    lay_sum = df.groupby([event,'X0'])['rechit_energy_w'].sum()
    lay_sum_avg = lay_sum.reset_index().groupby('X0')['rechit_energy_w'].median()
    # layers = df.groupby('X0').X0.mean()
    # esums = [df[df.X0 == X0].groupby(event).rechit_energy_w.sum() for layer in layers]
    # errors = [esum.std()/np.sqrt(len(esum)) for esum in esums]

    return lay_sum_avg#, errors

def do_COG(df, event = 'event'):
    df = df.copy()
    df['xx0'] = df['rechit_layer'].map(lay_X0)
    df['rE'] = df.xx0*(df.rechit_energy)
    etot = get_totE(df)
    cg = df.groupby(event).rE.sum()/(etot)
    return cg
