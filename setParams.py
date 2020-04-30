#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 13:18:35 2018

@author: amanaberra
"""
import numpy as np

def setParams(nrn_model_ver = "max4rP60c"):
    dt = 0.005
    record_dt = dt
    tstop = 1
    amp = 100 # E-field amp
    st_amp = 0 # current injection amp
    st_del = 0
    st_dur = tstop
    st_mode = 1 # 1 - DC, 2 - rect pulses, 3 - triangular pulses
    cell_model_id = 32 
    v_init = -70
    theta = 180
    phi = 0
    axon_num_nodes = 10
    ss_init = 1
    temp = 37 
    record_mode = 4 # 1 - record V and spike times from all x=0.5, 2 - record V and spike times from all x=0.5, x=1, 3 - record V from soma, spike times from all x=0.5, x = 1
    synapses_on = 0
    resave = 0
    load_potentials = 0
    scale_mainaxon_diam=1 
    if nrn_model_ver is "max4rP60c":
        # adult, rat with myelinated axon, using Zhu 2000 scaling
        myelinate_axon=1
        prune_axon=0
        replace_axon = 0 # 1 - just iseg, 2 - expand main ax terminal, 3- set terminals passive + expand main ax terminal 
        scale_axon_diam=1.322 # 1.322 using Zhu 2000 soma scaling
        scale_apic_diam=1.248 # 1.248
        scale_basal_diam=1.133 # 1.133
        scale_soma_area=1.322 # 1.322 using Zhu 2000 soma scaling
        scale_basal_L=1        
        min_myelinD=0.2 # um - matches paper
        max_myelin_order = 0
    elif nrn_model_ver is "maxH":
        # adult, human with myelinated axon using L2/3 PC rat BB: human Eyal 2018
        myelinate_axon=1
        prune_axon=0
        replace_axon=0
        scale_axon_diam=2.453 # use soma scaling
        scale_apic_diam=1.876
        scale_basal_diam=1.946
        scale_soma_area=2.453
        scale_basal_L=1.17 # using Romand 2011 L5 PC growth
        min_myelinD=0.2 # um - matches paper
        max_myelin_order=0
    elif nrn_model_ver is "umaxr":
        # P14 rat, unmyelinated, no scaling
        myelinate_axon=0
        prune_axon=0
        replace_axon=0
        scale_axon_diam=1
        scale_apic_diam=1
        scale_basal_diam=1
        scale_soma_area=1
        scale_basal_L=1
        min_myelinD=0.2 # um - matches paper
        max_myelin_order=0
    print "nrn_model_ver = ", nrn_model_ver
    Params = np.array([dt,tstop,v_init,amp,st_amp,st_del,st_dur,cell_model_id,\
            theta,phi,axon_num_nodes,ss_init,record_mode,\
            synapses_on,resave,st_mode,record_dt,load_potentials,replace_axon,\
            myelinate_axon,prune_axon,scale_axon_diam,scale_apic_diam,scale_basal_diam,\
            scale_soma_area,scale_basal_L,scale_mainaxon_diam,min_myelinD,\
            temp,max_myelin_order])
    return Params

    
    
    