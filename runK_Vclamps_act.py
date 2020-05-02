#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 23:38:24 2020

@author: amanaberra
"""

# Activation protocol for Na channels/current
import os
# os.system('nrnivmodl mechanisms')
# os.system('cp x86_64/special .')
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
import numpy as np
from math import pi
# VClamp testing code
from neuron import h #, gui # make sure default -NSTACK and -NFRAME are set to 100000,20000
from Vclamp import Vclamp_tester
from plot_Vclamp import plot_recs,plot_gmax, plot_tau
import pickle # for saving
from scipy.io import savemat
from scipy.ndimage import gaussian_filter
from fit_taus import fit_double_exp,fit_single_exp, single_exp
fig_folder = 'Figures/Vclamp'
data_folder = 'Data/Vclamp'
calc_tau = False # TODO fix time constant fitting
fc = None # Hz (fc = 10 kHz from Schmidt-Hieber 2010 for gaussian filter)
save_figs = True
save_data = True
#chan_names = ['MCna1']
#gbar_names = ['gna1bar']
#g_names = ['gna1']
all_channels = {'K_Tst':{'gbar':'gK_Tstbar','g':'gK_Tst','color':'b'}, # Blue Brain Transient K (Korngreen and Sakmann 2000)
                'K_Pst':{'gbar':'gK_Pstbar','g':'gK_Pst','color':'m'}, # Blue Brain Persistent K (Korngreen and Sakmann 2000)
                'Kv1':{'gbar':'gbar','g':'gkv1','color':'r'}, # n^8 Kv1 channel from axons (Kv1_axonal from Hallerman 2012)
                'SKv3_1':{'gbar':'gSKv3_1bar','g':'gSKv3_1','color':'g'}, # Blue Brain Kv3.1 (Shaw related)
                'Kv7':{'gbar':'gbar','g':'g','color':'k'}, # Kv7 axonal m-current Hallerman 2012 
                'Im':{'gbar':'gImbar','g':'gIm','color':'y'}, # Blue Brain M-urrent (Adams 1982)
                } 

sim_channels = ['Im','Kv7']
# sim_channels = ['Kv7','Im']
# TODO figure out why first current of Kv7 always has NaNs

channels = {k:all_channels[k] for k in sim_channels if k in all_channels}
#colors = ['k','r','b','g','m']
clamp_params = {'amp1':-120, 
                'dur1':100, 
                'amp2':-90,
                'dur2': 200,
                'amp3':0,
                'dur3':0,
                'v_stepi':-100,
                'v_stepf':50,
                'dv_step': 5,
                'amp_to_step': 'amp2', # activation protocol steps 2nd phase
                'rs':1e-5,# MOhm
                'x':0.5}
# NEURON settings
# Settings
curr_name = 'ik'
dt = 0.025 # 1 µs
Ek = -100 # 55 for Kole 2008
T = 37
q10 = 2.3 # rate constant coefficient
# Plot activation curves
fig2 = None # g/gmax curves
fig3 = None # taum curves
min_curr = 0 # 0.1% relative to cmin_global
# curr_names all the same
#for chan_name,gbar_name,g_name,colori in zip(chan_names,gbar_names,g_names,colors):
for chan_name,chan in channels.items():
    
    clamp_test = Vclamp_tester(chan_name=chan_name,curr_name=curr_name,\
                               gbar_name=chan['gbar'],g_name=chan['g'],\
                            clamp_params=clamp_params,dt=dt,T=T,q10=q10,ek=Ek)

    # Run vclamp simulations
    v_steps,t_vec,vs,currs,gs = clamp_test.run_protocol(clamp_params)
    cmax_global = np.nanmax([np.nanmax(abs(c)) for c in currs]) # abs min
    t_vec = t_vec - clamp_params['dur1'] # V step starts at 0
    # get tau fits for currents 
    if calc_tau:
        tau1 = np.zeros(len(v_steps))
        tau2 = np.zeros(len(v_steps))
        t_acts = []
        c_actfs = [] # 
        for i,c in enumerate(currs):
            t_act = t_vec[np.int(clamp_params['dur1']/dt):c.argmax()]
            c_act = c[np.int(clamp_params['dur1']/dt):c.argmax()]
            if c.max() > min_curr*cmax_global:
                try:
                    if fc is not None:
                        Fs = 1e3/dt # sampling freq
#                        c = 2
#                        sigma = np.sqrt(2*np.log(c))*Fs/(2*pi*fc)
                        sigma = Fs/(2*pi*fc)
                        c_act = gaussian_filter(c_act,sigma)
                    p0 = (1,1,0.1) # start at amp = 1 (norm.), tau = 1ms, delay = 1
                    bounds = ([0.9,0.005,0],[1.5,5,1])
                    tau1[i],fit = fit_single_exp(t_act,c_act/np.max(np.abs(c_act)),p0=p0,bounds=bounds)
                    tau2[i] = np.nan
                    c_actfi = single_exp(t_act,fit[0]*np.max(np.abs(c_act)),fit[1],fit[2])
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
    gmax = np.nanmax([np.nanmax(g) for g in gs])
    g_vecsn = [g/gmax for g in gs]
    # g_vecsn = gs
    
    fig1, ax1, ax2, ax3 = plot_recs(t_vec,vs,currs,g_vecsn) # plot v,curr, and g/gmax
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
            fig3,ax5 = plot_tau(v_steps,tau1,channels[chan_name]['color'],chan_name) # create fig
        else:
            plot_tau(v_steps,tau1,channels[chan_name]['color'],chan_name,fig3,ax5) # add to fig
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
    if calc_tau:
        fig3_name = os.path.join(fig_folder,''.join(chan_namei + '_' for chan_namei in sim_channels) + 'taum')
        fig3.savefig(fig3_name,dpi=200)

    fig2_name = os.path.join(fig_folder,''.join(chan_namei + '_' for chan_namei in sim_channels) + 'gmax')
    fig2.savefig(fig2_name,dpi=200)