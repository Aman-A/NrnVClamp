#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 09:38:12 2018

@author: amanaberra
"""
import matplotlib
#matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
from numpy import ones,array,nanmax
# centers time axis on t_start (sets to t = 0)
def plot_recs(t_vec,vs,currs,gs,x_lims=None):
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    for v_vec,i_vec,g_vec in zip(vs,currs,gs):   
        ax1.plot(t_vec,v_vec,linewidth=0.5); 
        ax2.plot(t_vec,i_vec,linewidth=0.5);         
        ax3.plot(t_vec,g_vec,linewidth=0.5); 
#        ax1.plot(t_vec,v_vec,color='k',linewidth=0.5); 
#        ax2.plot(t_vec,i_vec,color='k',linewidth=0.5);         
#        ax3.plot(t_vec,g_vec,color='k',linewidth=0.5); 
    ax1.set_ylabel('V (mV)')
    ax2.set_ylabel('i (mA/cm2)')
    ax3.set_ylabel('g/gmax')        
    ax3.set_xlabel('time (ms)')
    if x_lims is None:
        x_lims = (-0.1,5)        
    ax1.set_xlim(x_lims)
    ax2.set_xlim(x_lims)
    ax3.set_xlim(x_lims)
    plt.draw()
    return fig, ax1, ax2, ax3
    
def plot_gmax(v_steps,gs,color='k',label=None,fig=None,ax=None):
    gmaxs = [nanmax(g) for g in gs]
    gmaxs_max = nanmax(gmaxs) 
    if not fig:        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(array([v_steps[0],v_steps[-1]]),0.5*ones((2,)),'--k')
    ax.plot(v_steps,gmaxs/gmaxs_max,'-o',color=color,label=label)    
    ax.set_xlim(v_steps[0],v_steps[-1])
    ax.set_xlabel('V (mV)')
    ax.set_ylabel('g/gmax')
    ax.legend()
    plt.draw()
    return fig,ax, gmaxs

def plot_tau(v_steps,tau1,color='k',label=None,fig=None,ax=None,tau2=None):
    if not fig:        
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.plot(v_steps,tau1*1e3,'-o',color=color,markerfacecolor=color,label=label+": tau1")
    if tau2 is not None:
        ax.plot(v_steps,tau2*1e3,'-o',color=color,markerfacecolor='w',label=label+": tau2")
    ax.set_xlabel('V (mV)')
    ax.set_ylabel('tau (us)')
    ax.legend()
    plt.draw()
    return fig,ax