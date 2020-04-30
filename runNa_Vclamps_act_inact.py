#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 12:48:41 2020

@author: amanaberra
"""
# Activation protocol for Na channels/current
import os
#os.system('nrnivmodl mechanisms')
#os.system('cp x86_64/special .')
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import numpy as np
# VClamp testing code
from neuron import h #, gui # make sure default -NSTACK and -NFRAME are set to 100000,20000
from Vclamp import Vclamp_tester, run_VClamp_act
from plot_Vclamp import plot_recs,plot_gmax, plot_tau
import pickle # for saving
from scipy.io import savemat
from fit_taus import fit_double_exp,fit_single_exp, single_exp
fig_folder = 'Figures/Vclamp'
data_folder = 'Data/Vclamp'
calc_tau = True
save_figs = False
save_data = False
#chan_names = ['MCna1']
#gbar_names = ['gna1bar']
#g_names = ['gna1']
all_channels = {'NaTa_t':{'gbar':'gNaTa_tbar','g':'gNaTa_t','color':'b'}, # Blue Brain axonal (Colbert & Pan 2002) 
            'na6_gp':{'gbar':'gbar','g':'g','color':'r'}, # Globus Pallidus Nav1.6 (Mercer 2007)
            'na16':{'gbar':'gbar','g':'gna','color':'b'}, # AIS Nav1.6 (Hu 2009)
            'na12':{'gbar':'gbar','g':'gna','color':'g'}, # AIS/soma Nav1.2 (Hu 2009)
            'NaTs2_t':{'gbar':'gNaTs2_tbar','g':'gNaTs2_t','color':'m'}, # Blue Brain somatic (Colbert & Pan 2002) 
            'naf':{'gbar':'gbar','g':'gna','color':'k'}, # Fast Na (Traub 2003)
            'nax_kole2008':{'gbar':'gbar','g':'gna','color':'k'}, # Axonal Nav (Kole 2008)
            'na_kole2008':{'gbar':'gbar','g':'gna','color':'k'}, # Somatic Nav (Kole 2008)
            'nafJonas':{'gbar':'gbar','g':'gna','color':'r'}, # Fast Na in MFB (based on Engel 2005, impl in Schmidt-Hieber 2008)
            'NaV':{'gbar':'gbar','g':'gna','color':'r'}, # Mouse Na from Allen models, based on (Carter 2012) (37° rec in mouse CA1 pyramids)
            'nax':{'gbar':'gbar','g':'gna','color':'r'}, # Axonal Na 8st model from Schmidt-Hieber 2010
            'na8st':{'gbar':'gbar','g':'gna','color':'k'}, # Somatic Na 8st model from Schmidt-Hieber 2010
            'MCna1':{'gbar':'gna1bar','g':'gna1','color':'r'}} # 2-closed, 1 open state Na channel (Baranauskas 2006)
sim_channels = ['NaTs2_t','NaTa_t','na8st','nax'] # 'na16','na12','nax',
#sim_channels = ['nax']
channels = {k:all_channels[k] for k in sim_channels if k in all_channels}
#colors = ['k','r','b','g','m']
curr_name = 'ina'
clamp_params = {'amp1':-120, 
                'dur1':50, 
                'amp2':-120,
                'dur2': 20,
                'amp3':0,
                'dur3':0,
                'v_stepi':-60,
                'v_stepf':20,
                'dv_step': 5,
                'rs':1e-3,
                'x':0.5}
# NEURON settings
# Settings
dt = 0.001 # 1 µs
Ena = 50 # 55 for Kole 2008
T = 37
q10 = 2.3 # rate constant coefficient

# Plot activation curves
fig2 = None # g/gmax curves
fig3 = None # taum curves
min_curr = 0.001 # relative to cmin_global
# curr_names all the same
#for chan_name,gbar_name,g_name,colori in zip(chan_names,gbar_names,g_names,colors):
for chan_name,_ in channels.items():
    
    clamp_test = Vclamp_tester(chan_name=chan_name,curr_name=curr_name,gbar_name=channels[chan_name]['gbar'],\
                                   g_name=channels[chan_name]['g'],clamp_params=clamp_params,dt=dt,T=T,q10=q10,ena=Ena)

    # Run vclamp simulations
    v_steps,t_vec,vs,currs,gs = run_VClamp_act(clamp_test, clamp_params) 
    cmin_global = min([min(c) for c in currs])    
    t_vec = t_vec - clamp_params['dur1'] # V step starts at 0
    # get tau fits for currents 
    if calc_tau:
        tau1 = np.zeros(len(v_steps))
        tau2 = np.zeros(len(v_steps))
        t_acts = []
        c_actfs = [] # 
        for i,c in enumerate(currs):
            t_act = t_vec[np.int(clamp_params['dur1']/dt):c.argmin()]
            c_act = c[np.int(clamp_params['dur1']/dt):c.argmin()]
            if c.min() < min_curr*cmin_global:
                try:
                    tau1[i],fit = fit_single_exp(t_act,c_act)
                    tau2[i] = np.nan
                    c_actfi = single_exp(t_act,fit[0],fit[1],fit[2])
                    t_acts.append(t_act)
                    c_actfs.append(c_actfi)
#                    print("{}: single exp for v = {}".format(chan_name,v_steps[i]))
#                    tau1[i],tau2[i]  = fit_double_exp(t_vec[c.argmin():],c[c.argmin():])
                except RuntimeError:
                    tau1[i],tau2[i],fit  = fit_double_exp(t_vec[c.argmin():],c[c.argmin():])
#                    tau1[i] = fit_single_exp(t_vec[c.argmin():],c[c.argmin():])
#                    tau2[i] = np.nan
                    print("{}: double exp for v = {}".format(chan_name,v_steps[i]))
            else:
                print("{}: current smaller than {} for v = {}".format(chan_name, min_curr,v_steps[i]))
                tau1[i] = np.nan; tau2[i] = np.nan # ignore these values    
        print("Fit current curves to exponential")
    # normalize all g vecs to max g
    gmax = max([max(g) for g in gs])
    g_vecsn = [g/gmax for g in gs]
    
    fig1, ax1, ax2, ax3 = plot_recs(t_vec,vs,currs,gs) # plot v,curr, and g/gmax
    ax1.set_title(chan_name)
    # save recordings for this current
    if save_figs:
        fig1_name = os.path.join(fig_folder,chan_name + '_recs')
        fig1.savefig(fig1_name,dpi=200)
    # Add gmax curve to fig2
    if not fig2:
        fig2,ax4, gmaxs = plot_gmax(v_steps,gs,channels[chan_name]['color'],chan_name) # create fig
    else:
        _,_,gmaxs = plot_gmax(v_steps,gs,channels[chan_name]['color'],chan_name,fig2,ax4) # add to fig
    # Add tau curve to fig3
    if calc_tau:
        # add fits to current plot
        for tf,cf in zip(t_acts,c_actfs):
            ax2.plot(tf,cf,'k--')
        if not fig3:
            fig3,ax5 = plot_tau(v_steps,tau1,tau2,channels[chan_name]['color'],chan_name) # create fig
        else:
            plot_tau(v_steps,tau1,tau2,channels[chan_name]['color'],chan_name,fig3,ax5) # add to fig
    datai = {}    
    datai['v_steps'] = v_steps
    datai['gs'] = gs
    datai['currs'] = currs
    datai['vs'] = vs
    datai['t_vec'] = t_vec
    datai['gmaxs'] = gmaxs
    datai['clamp_params'] = clamp_params
    datai['params'] = clamp_test.params
    if save_data:
        with open(os.path.join(data_folder,chan_name+'_act_data_T_{}.pkl'.format(T)),'wb') as pickle_file:
            pickle.dump(datai,pickle_file)
        savemat(os.path.join(data_folder,chan_name+'_data.mat'),datai)
#    taus_all[chan_name] = {}
#    taus_all[chan_name]['tau1'] = tau1
#    taus_all[chan_name]['tau2'] = tau2
# Save single gmax figure
plt.show()
if save_figs:    
    fig2_name = os.path.join(fig_folder,''.join(chan_namei + '_' for chan_namei in chan_names) + 'gmax')
    fig2.savefig(fig2_name,dpi=200)