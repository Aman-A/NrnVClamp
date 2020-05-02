#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 15:42:25 2018

@author: amanaberra
"""
import os
import pickle 
import numpy as np 
from plot_Vclamp import plot_recs,plot_gmax#,plot_tau
#chan_names = ['nafJonas2']
# chan_names = ['na8st','nax','MCna1']
# chan_names = ['SKv3_1']
chan_names = ['Im','Kv7']
# chan_names = ['nax','nav6','NaTa_t']
T = 37
rec_x_lims=(0,100)
#chan_names = ['naf','na16']
#chan_names = ['NaTa_t','naf']
colors = [(0, 0.4470, 0.7410),(0.8500,0.3250,0.0980),(0.9290,0.6940, 0.1250),
          (0.494,0.184,0.556),(0.466,0.674,0.188)]
colors = colors[0:len(chan_names)]
fig_folder = 'Figures/Vclamp'
data_folder = 'Data/Vclamp'
save_figs = True
# Load and plot data
data = {}
fig2 = None
for chan_name, colori in zip(chan_names,colors):
    with open(os.path.join(data_folder,chan_name+'_act_data_T_{}.pkl'.format(T)),'rb') as pickle_file:
            data[chan_name] = pickle.load(pickle_file)
    
    datai = data[chan_name]
    gmax = np.nanmax([np.nanmax(g) for g in datai['gs']])
    g_vecsn = [g/gmax for g in datai['gs']]
    # x_lims = (-5,datai['clamp_params']['dur2'])
    # Fig 1
    fig1, ax1, ax2, ax3 = plot_recs(datai['t_vec'],datai['vs'],datai['currs'],g_vecsn,rec_x_lims) # plot v,curr, and g/gmax
    ax1.set_title(chan_name)
    if save_figs:
        fig1_name = os.path.join(fig_folder,chan_name + '_recs')
        fig1.savefig(fig1_name,dpi=200)
    # Fig 2
    if not fig2:
        fig2,ax4, _ = plot_gmax(datai['v_steps'],datai['gs'],colori,chan_name) # create fig
    else:
        plot_gmax(datai['v_steps'],datai['gs'],colori,chan_name,fig2,ax4) # add to fig
            
if save_figs:    
    fig2_name = os.path.join(fig_folder,''.join(chan_namei + '_' for chan_namei in chan_names) + 'gmax')
    fig2.savefig(fig2_name,dpi=200)            
