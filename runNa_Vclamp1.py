#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This code runs 
Created on Mon Jun 18 09:49:09 2018

@author: amanaberra
"""
import os
#os.system('nrnivmodl mechanisms')
#os.system('cp x86_64/special .')
import matplotlib.pyplot as plt
#matplotlib.use('qt5agg')
# VClamp testing code
#from neuron import h, gui # make sure default -NSTACK and -NFRAME are set to 100000,20000
from Vclamp import Vclamp_tester, run_VClamp_act
from plot_Vclamp import plot_recs,plot_gmax

fig_folder = 'Figures/Vclamp'
save_fig = False
chan_name = 'na16' # NaTa_t
curr_name = 'ina'
gbar_name = 'gbar' # gNaTa_tbar
g_name='gna' # gNaTa_t
clamp_params = {'amp1':-120, 
                'dur1':40, 
                'amp2':-80,
                'dur2': 30,
                'amp3':0,
                'dur3':0,
                'v_stepi':-120,
                'v_stepf':40,
                'dv_step': 5,
                'rs':1,
                'x':0.5}
# NEURON settings
# Settings
dt = 0.005 # 1 µs
Ena = 50
T = 37
q10 = 2.3 # rate constant coefficient

clamp_test = Vclamp_tester(chan_name=chan_name,curr_name=curr_name,gbar_name=gbar_name,\
                               g_name=g_name,clamp_params=clamp_params,dt=dt,T=T,q10=q10,ena=Ena)

# Run vclamp simulations
v_steps,t_vec,vs,currs,gs = run_VClamp_act(clamp_test, clamp_params) 


t_vec = t_vec - clamp_params['dur1'] # V step starts at 0
# normalize all g vecs to max g
gmax = max([max(g) for g in gs])
g_vecsn = [g/gmax for g in gs]
#g_vecsn = [g.to_python(np.zeros(len(g)))/gmax for g in gs] #g_vecn = g_vecn/g_vecn.max()
fig1, ax1, ax2, ax3 = plot_recs(t_vec,vs,currs,g_vecsn) # plot v,curr, and g/gmax
# Plot gmax vs v
fig2,ax4, gmaxs = plot_gmax(v_steps,gs,'k',chan_name)
plt.show()
# Save figures
if save_fig:
    fig1_name = os.path.join(fig_folder,chan_name + '_recs')
    fig2_name = os.path.join(fig_folder,chan_name + '_gmax')
    fig1.savefig(fig1_name,dpi=200)
    fig2.savefig(fig2_name,dpi=200)