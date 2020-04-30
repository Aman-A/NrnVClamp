#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 13:58:16 2018

@author: amanaberra
"""

import numpy as np
import scipy
import scipy.io
# Outputs normalized TMS waveform E and time vector t
# [t,E] = TMWwave(dt,tstop,delay,mode,plot_wave)

def TMS_wave(dt=0.005,tstop=1,delay=0.005,mode=3,plot_wave=1):                  
    # Get data
    data = scipy.io.loadmat('TMSwaves.mat');
    if mode is 3: # MagproX100 Monophasic waveform - Recorded by Lari Koponen (6/12/18)            
        trec = data['tm']; Erec = data['Erec_m'];                        
        title_string = 'Normalized E - MagproX100 Monophasic, delay = {:.1f} ms'.format(delay);            
    elif mode is 4: # MagproX100 Biphasic waveform - Recorded by Lari Koponen (6/12/18)
        trec = data.tb; Erec = data['Erec_b'];            
        title_string = 'Normalized E - MagproX100 Biphasic, delay = {:.1f} ms'.format(delay);            
    elif mode is 5: # MagproX100 Half-sine waveform - Recorded by Lari Koponen (6/12/18)            
        trec = data.th; Erec = data['Erec_h'];            
        title_string = 'Normalized E - MagproX100 Half sine, delay = {:.1f} ms'.format(delay);            
    # Downsample for simulation specific time step
    DT = np.mean(np.diff(trec,axis=0));
    sample_factor = int(dt/DT);
    if sample_factor < 1: # dt chosen is smaller than DT of recording
       sample_factor = 1; # don't downsample 
       print 'dt chosen smaller than recording, set to {:f} ms\n'.format(DT) 
    trec = trec.flatten(); Erec = Erec.flatten()
    # construct final vectors            
    if delay >= dt:      
        tvec = np.append([0,delay-dt],np.arange(delay,tstop+dt,dt))# include delay           
        Evec = np.append([0,0],np.interp(np.arange(0,trec[-1],dt),trec,Erec),0);
    else:
        tvec = np.arange(0,tstop+dt,dt); # pulse starts at t = 0
        Evec = np.append(np.interp(np.arange(0,trec[-1],dt),trec,Erec),0);                
    
    tvec = np.append(tvec,np.arange(tvec[-1]+dt,tstop+dt,dt)) # at end point to end at 0
           
    if tvec.size > Evec.size:
        Evec = np.append(Evec,np.zeros((tvec.size-Evec.size,1))); # add trailing zeros
    else:
        Evec = Evec[0:tvec.size];
    
    if plot_wave:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='w');
        ax = plt.subplot(111)
        ax.plot(tvec,Evec,color='k',linewidth=2)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        plt.title(title_string);
        plt.xlabel('time (ms)'); plt.ylabel('E (normalized)')
        plt.xlim([0, tstop]);    
        plt.show()
    return tvec, Evec
