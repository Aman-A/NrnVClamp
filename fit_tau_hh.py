#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 00:27:31 2018

@author: amanaberra
"""
import os, pickle
import numpy as np
#from scipy.optimize import curve_fit,leaste_squares
from lmfit import Model,Parameters
import matplotlib.pyplot as plt
data_folder = 'Data/Vclamp'
chan_name='NaTa_t'
T = 37
with open(os.path.join(data_folder,chan_name+'_act_data_T_{}.pkl'.format(T)),'rb') as pickle_file:
            data = pickle.load(pickle_file)

t_vec = data['t_vec']
gs = data['gs']
vs = data['vs']
v_steps = data['v_steps']
gmaxs = data['gmaxs']
# alpha/beta
def rate(v,R,vth,k):
    return R*(v-vth)/(1-np.exp(-(v-vth)/k))
# minf/hinf
def xinf(v,Ra,vtha,ka,Rb,vthb,kb,n):
    alpha_m = rate(v,Ra,vtha,ka)
    beta_m = rate(v,Rb,vthb,kb)
    return np.power(alpha_m/(alpha_m + beta_m),n)
params = Parameters()
# add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
params.add_many(('Ra',0.182,True,1e-5,1,None,None),
                ('vtha',-32,True,-100,0,None,None),
                ('ka',0,True,0.1,10,None,None),
                ('Rb',0.182,True,1e-5,1,None,None),
                ('vthb',-32,True,-100,0,None,None),
                ('kb',0,True,0.1,10,None,None),
                ('n',3,False,1,3,None,None))
def xtau(v,Ra,tha,ka,Rb,vthb,kb):
    alpha_m = rate(v,Ra,vtha,ka)
    beta_m = rate(v,Rb,vthb,kb)
    return 1/(alpha_m + beta_m)
# returns g(t) for voltage step (constant v) 
def gna(t,gbar,m0,minf,taum,h0,hinf,tauh):
    m_t = m0 - (m0 - minf)*(1 - np.exp(-t/taum))
    h_t = h0 - (h0 - hinf)*(1 - np.exp(-t/tauh))
    return gbar*np.power(m_t,3)*h_t

def gna_m3(t,gbar,m0,minf,taum):
    m_t = m0 - (m0 - minf)*(1 - np.exp(-t/taum))    
    return gbar*np.power(m_t,3)

taum = np.zeros(len(v_steps))
minf = np.zeros(len(v_steps))
tauh = np.zeros(len(v_steps))
hinf = np.zeros(len(v_steps))
gbars = np.zeros(len(v_steps))

gmodel = Model(gna)
params = Parameters()
# add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
#params.add_many(('gbar',0.012,True,1e-4,100,None,None),
#                ('m0',0,True,0,0.1,None,None),
#                ('minf',0.5,True,0,1,None,None),
#                ('taum',0.1,True,1e-3,0.5,None,None))
params.add_many(('gbar',0.012,True,1e-4,100,None,None),
                ('m0',0,True,0,0.1,None,None),
                ('minf',0.5,True,0,1,None,None),
                ('taum',0.1,True,1e-3,10,None,None),
                ('h0',1,True,0.8,1,None,None),
                ('hinf',0.2,True,0,1,None,None),
                ('tauh',1,True,1e-3,50,None,None))
#gmodel.set_param_hint('gbar',min=1e-4,max=100) # set bounds
#gmodel.set_param_hint('m0',min=0,max=1)
#gmodel.set_param_hint('minf',min=0,max=1)
#gmodel.set_param_hint('taum',min=1e-3,max=30)
#gmodel.set_param_hint('h0',min=0,max=1)
#gmodel.set_param_hint('hinf',min=0,max=1)
#gmodel.set_param_hint('tauh',min=1e-3,max=100)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for i,gi in enumerate(gs):
    ind_start = np.nonzero(t_vec >=0)[0][0]
    ti = t_vec[ind_start:gi.argmax()]
    gi = gi[ind_start:gi.argmax()] # g after voltage step    
#    fits,covs = curve_fit(gna,ti,gi,p0=[0.012,0.5,0.5,0.1,0.5,0.5,1],\
#                          bounds=([1e-4,0,0,1e-3,0,0,1e-3],[100,1,1,30,1,1,100]),\
#                          method='dogbox',maxfev=2000)
    result = gmodel.fit(gi,params,t=ti)
    gbars[i] = result.params['gbar'].value
    minf[i] = result.params['minf'].value
    taum[i] = result.params['taum'].value
#    hinf[i] = result.params['hinf'].value
#    tauh[i] = result.params['tauh'].value
    ax1.plot(ti,gi,'-')
    ax1.plot(ti,gmodel.eval(result.params,t=ti),'.')
    
fig = plt.figure(); 
ax1 = fig.add_subplot(211)
ax1.plot(v_steps,minf,label='$m_{\infty}$')
#ax1.plot(v_steps,hinf,label='$h_{\infty}$')
ax1.set_ylabel('State'); ax1.legend()
ax2 = fig.add_subplot(212)
ax2.plot(v_steps,taum,label='$\\tau_{m}$')
#ax2.plot(v_steps,tauh,label='$\\tau_{h}$')
ax2.set_ylabel('Time constant (ms)')
ax2.set_xlabel('Vm (mV)'); ax2.legend()

fig2 = plt.figure()
ax3 = fig2.add_subplot()
ax3.plot(v_steps,gbars)
ax3.set_xlabel('Vm (mV)')
ax3.set_ylabel('gbar (S/cm2)')
#fig = plt.figure(); ax = fig.add_subplot(111)
#ax.plot(ti,gi,label='Recording')
#ax.plot(ti,gna(ti,fits[0],fits[1],fits[2],fits[3],fits[4],fits[5],fits[6]),'o',label='Fit',markersize=2)
#ax.set_xlabel('time (ms)')
#ax.set_ylabel('gna (S/cm2)')
#plt.legend()