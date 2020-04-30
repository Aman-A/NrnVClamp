#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 03:55:46 2018

@author: amanaberra
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 00:27:31 2018

@author: amanaberra
"""
import os, pickle
import numpy as np
#from scipy.optimize import curve_fit
from lmfit import Model,Parameters
import matplotlib.pyplot as plt
data_folder = 'Data/Vclamp'
fig_folder = 'Figures/Vclamp'
chan_names=['na6_gp']
#chan_names=['NaTa_t','na6_gp','na16','nax_kole2008','nafJonas','nax','NaV','MCna1']
T = 37
save_data=False
save_figs=False
over_write=True
# Fit function
def exp_curr(t,A,tau,td):
    return np.power(A*(1-np.exp(-(t-td)/tau)),3)
# Define model    
model = Model(exp_curr,nan_policy='omit')
params = Parameters()
# add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
params.add_many(('A',-1,True,None,0,None,None),
                ('tau',0.1,True,0,0.3,None,None),
                ('td',0.1,True,0,1,None,None))
min_curr = -1e-7

fig = plt.figure(); 
ax = fig.add_subplot(111)
for chan_name in chan_names:
    
    with open(os.path.join(data_folder,chan_name+'_act_data_T_{}.pkl'.format(T)),'rb') as pickle_file:
                data = pickle.load(pickle_file)
    
    if 'taum' not in data.keys() or over_write: # don't rerun unless over_write == True
        print "Fitting tau for channel: ", chan_name
        t_vec = data['t_vec']
        ind_start = np.nonzero(t_vec >=0)[0][0]
        currs = data['currs']
        vs = data['vs']
        v_steps = data['v_steps']    
        Ena = data['params']['Ex']['ena']
        T = data['params']['T']
        taum = np.zeros(len(v_steps))
    
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        for i,ci in enumerate(currs):
            if (ci.min() < min_curr and ci.argmin() > ind_start and v_steps[i] != Ena):
                ti_rise = t_vec[ind_start:ci.argmin()]
                ci_rise = ci[ind_start:ci.argmin()] # g after voltage step    
            #    fits,covs = curve_fit(exp_curr,ti_rise,ci_rise,p0=[-1,0.1,0],\
            #                          bounds=([-50,1e-3,0],[0,30,10]),\
            #                          method='dogbox',maxfev=2000)
                params['A'].init_value = ci_rise.min()
                result = model.fit(ci_rise,t=ti_rise,params=params)    
                taum[i] = result.params['tau'].value
            #    taum[i] = fits[1]
                ax2.plot(ti_rise,ci_rise,'k-')
            #    ax.plot(ti_rise,exp_curr(ti_rise,fits[0],fits[1],fits[2]),'r.')
                ax2.plot(ti_rise,model.eval(result.params,t=ti_rise),'.',label=v_steps[i])
            else:
                taum[i] = np.nan
        ax2.set_title(chan_name)        
        ax2.legend()
        if save_data:
            datai = data
            datai['taum'] = taum
            datai['fit_function'] = 'Time-shifted single exponential'
            with open(os.path.join(data_folder,chan_name+'_act_data_T_{}.pkl'.format(T)),'wb') as pickle_file:
                print "Saving to",chan_name+'_act_data.pkl'
                pickle.dump(datai,pickle_file)
    else:
        print "Loading tau for channel: ", chan_name
        taum = data['taum']
        v_steps = data['v_steps']        
    # plot tau vs v on same axis for all channels    
    ax.plot(v_steps,taum*1e3,label=chan_name)
    ax.set_ylabel('Time constant ($\\mu s$)',usetex=True)
    ax.set_xlabel('Vm (mV)');
ax.legend()

if save_figs:    
    fig_name = os.path.join(fig_folder,''.join(chan_namei + '_' for chan_namei in chan_names) + 'taum')
    fig.savefig(fig_name,dpi=200)
    fig2_name = os.path.join(fig_folder,''.join(chan_namei + '_' for chan_namei in chan_names) + 'act_curr_fits')
    fig2.savefig(fig2_name,dpi=200)
    