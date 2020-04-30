#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 06:21:34 2018

@author: amanaberra
"""
import numpy as np
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import pickle
plt.rc('font',size=10)
execfile("run_stim.py")
#mtau_scaling=[0.05,0.1,0.23,0.28,0.3,0.31,0.5,0.7,0.8,1]    
mtau_scaling = np.arange(0.1,1.1,0.1)
htau_scaling = np.arange(0.1,1.1,0.1)
dv_dt_thresh = 50


def run_scale_mtau(mtau_scaling,threshEs):            
    
    v_all = []
    h.tstop = 3 # reset to 1ms
    fig = plt.figure(); ax = fig.add_subplot(111)
    for i,mtau_s in enumerate(mtau_scaling):
        print "Scaling mtau to", mtau_s
        change_mtau(mtau_s,axonal)
        h.AMP = threshEs[i]
        print "Threshold E =", h.AMP
        stimul(plot_on=True,v_label="{:.2f}".format(mtau_s),fig=fig,ax=ax)        
        v_all.append(v_term.to_python())
    return fig, ax


def plot_mtau_thresh(mtau_scaling,threshEs,savefig=False):
    per_diff = 100*(threshEs - threshEs.max())/threshEs.max()
    fig = plt.figure(); 
    ax1 = fig.add_subplot(111); 
    ax2 = ax1.twinx()
    ax1.plot(mtau_scaling,threshEs,'-o',color='k')
    ax2.plot(mtau_scaling,per_diff,'-o',color='k')
    ax1.set_ylabel('Threshold E (V/m)')
    ax2.set_ylabel('Percent change (%)')    
    ax1.set_xlabel("$\\tau_{m}$ scaling",usetex=True)
    plt.show()
    fig.tight_layout()
    if savefig:
        fig.savefig("mtau_scaling",dpi=200)
    return fig,ax1,ax2

# input just AP waveform
def calc_Vm_thresh(t,v,dv_dt_thresh=50):
    dvdt = np.diff(v)/np.diff(t)
    td = t[:-1]
    th_ind = np.nonzero(dvdt >= dv_dt_thresh)[0][0]
    if dvdt[th_ind] > dv_dt_thresh: # linear interpolation (slightly above threshold)
        t1 = td[th_ind-1]; t2 = td[th_ind]
        dv1 = dvdt[th_ind-1]; dv2 = dvdt[th_ind]
        th = np.interp(dv_dt_thresh,[dv1, dv2],[t1,t2]) # threshold time point 
        v1 = v[th_ind-1]; v2 = v[th_ind]
        vth = np.interp(th,[t1,t2],[v1,v2])
    else: # data point is exactly at threshold
        th = td[th_ind]
        vth = v[th_ind]
    return vth, th # threshold Vm, time of threshold crossing

    
    
def scale_mtau(mtau_scaling):            
    threshEs = np.zeros(len(mtau_scaling))
    spktimes = np.zeros(len(mtau_scaling))
    vthresh = np.zeros(len(mtau_scaling))
    v_all = []
    fig = plt.figure(); ax = fig.add_subplot(111)
    plt.show()
    for i,mtau_s in enumerate(mtau_scaling):
        print "Scaling mtau to", mtau_s
        change_mtau(mtau_s,axonal)
        h.AMP = 100
        h.tstop = 1 # reset to 1ms
        h.find_thresh()
        # get threshold
        threshEs[i] = h.threshE # V/m        
        h.AMP = h.threshE
        print "Threshold E =", h.threshE
        h.tstop = 2.5
        stimul(v_label="{:.2f}".format(mtau_s),fig=fig,ax=ax)        
        # get Vm recording
        v_all.append(v_term.to_python())
        # get time of first threshold crossing
        spktimes[i] = h.ap_times_all.x[0] 
        ti = t.to_python(np.zeros(len(t)))
        #ti0 = ti[:int(np.nonzero(abs(ti-h.ap_times_all.x[0]) < 1e-6)[0][0])]
        # get vm before threshold crossing
        vi = v_term.to_python(np.zeros(len(v_term)))[:int(np.nonzero(abs(ti-h.ap_times_all.x[0]) < 1e-6)[0][0])]
        vthresh[i] = vi[argrelextrema(vi,np.greater)[0]].max() # get (largest) local max
                
    mtau_scaling = np.array(mtau_scaling)
    return threshEs,v_all,spktimes,vthresh

def plot_vall(mtau_scaling,t,v_all,savefig=False):
    fig = plt.figure(); ax = fig.add_subplot(111)
    for mtau_s,v in zip(mtau_scaling,v_all):
        ax.plot(t,v,label="{:.2f}".format(mtau_s))
    ax.set_xlabel('time (ms)')
    ax.set_ylabel('$V_{m}$ (mV)',usetex=True)
    ax.set_xlim(0,t.x[int(t.size()-1)])
    plt.legend()
    fig.tight_layout()
    if savefig:
        fig.savefig("Vm_mtau_scaling",dpi=200)
    return fig,ax
        

#threshEs, spktimes,vthresh,v_all = scale_mtau(mtau_scaling)
#fig,ax = plot_vall(mtau_scaling,t,v_all,savefig=True)
#fig2,ax2,ax3 = plot_mtau_thresh(mtau_scaling,threshEs,savefig=True)
    
def scale_mtau_htau(mtau_scaling,htau_scaling):            
    threshEs = np.zeros((len(htau_scaling),len(mtau_scaling)))    
    v_all = []
    for i,htau_s in enumerate(htau_scaling):
        change_htau(htau_s,axonal)
        v_all.append([])    # ith htau_scaling   
        v_alli = []
        for j,mtau_s in enumerate(mtau_scaling):
            change_mtau(mtau_s,axonal) 
            print "htau =", htau_s, "mtau =", mtau_s
            h.AMP = 100
            h.tstop = 1 # reset to 1ms
            h.find_thresh()
            # get threshold
            threshEs[i,j] = h.threshE # V/m 
            print "Threshold E =", h.threshE
            h.AMP = h.threshE
            h.tstop = 2.5            
            stimul(plot_on=False)            
            v_alli.append(v_term.to_python()) # append jth mtau_scaling
        v_all[i].append(v_alli) # append all to ith row of v_all
    return threshEs,v_all
                
def plotvall_mh_scaling(t,v_all,mtau_scaling,htau_scaling):
    fig = plt.figure(); 
    for i,htau_s in enumerate(htau_scaling):        
        for j,mtau_s in enumerate(mtau_scaling):
            axij = fig.add_subplot(len(htau_scaling),len(mtau_scaling),1+np.ravel_multi_index((i,j),\
                                      (len(htau_scaling),len(mtau_scaling))))
            axij.plot(t,v_all[i][0][j])
            axij.set_title("$\\tau_{h}=%.1f, \\tau_{m}=%.1f$" % (htau_s,mtau_s),usetex=True)
            axij.set_xlim(0,2)
            axij.set_ylim(-90,60)
            if i == len(htau_scaling)-1: # bottom row
                axij.set_xlabel('time (ms)')
            else:
                axij.set_xticklabels([])                
            if j == 0: # leftmost column
                axij.set_ylabel('$V_{m}$ (mV)',usetex=True)
            else:
                axij.set_yticklabels([])                
    return fig

def plotthresh_mh_scaling(threshEs,mtau_scaling,htau_scaling):
    m = mtau_scaling; h = htau_scaling
    t = threshEs;    
    cmap = "viridis"
    fig = plt.figure(); ax = fig.add_subplot(111)
    imj = ax.imshow(threshEs,interpolation='bilinear',extent=[0,mtau_scaling[-1],0,htau_scaling[-1]],cmap=cmap)
    fig.colorbar(imj)
    ax.set_aspect('auto')
    ax.set_xlabel('$\\tau_{m}$ scaling',usetex=True)
    ax.set_ylabel('$\\tau_{h}$ scaling',usetex=True)
    ax.set_title('Threshold E (V/m)')
    return fig
    
    
threshEs,v_all = scale_mtau_htau(mtau_scaling,htau_scaling)

data={}
data['t'] = t
data['mtau_scaling'] = mtau_scaling
data['htau_scaling'] = htau_scaling
data['v_all'] = v_all
data['threshEs'] = threshEs
with open('Data/v_all_mtau_htau','wb') as pickle_file:
    pickle.dump(data,pickle_file) 
    
fig = plotvall_mh_scaling(t,v_all,mtau_scaling,htau_scaling)
fig2 = plotthresh_mh_scaling(threshEs,mtau_scaling,htau_scaling)


#threshEsm,v_allm,spktimes,vthresh = scale_mtau(mtau_scaling)