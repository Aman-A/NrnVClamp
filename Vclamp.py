#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 12:01:36 2018

@author: amanaberra
"""

from neuron import h # make sure default -NSTACK and -NFRAME are set to 100000,20000
h.load_file('stdrun.hoc')
import numpy as np

class Vclamp_tester:
    def __init__(self,chan_name='NaTa_t',curr_name='ina',gmax=1,gbar_name = 'gNaTa_tbar', \
                 g_name='gNaTa_t', clamp_params=dict(), dt=0.001,T=37,q10=2.3,**kwargs):
        if kwargs:            
            Ex = kwargs # ena,ek, etc.
        else:
            Ex = dict(ena=50) # default Ena for NaTa_t"            
        self.chan_name = chan_name
        self.curr_name = curr_name
        self.gbar_name = gbar_name + '_' + chan_name # add suffix
        self.g_name = g_name + '_' + chan_name # add suffix
        self.gmax = gmax
        self.params = dict(dt=dt,T=T,q10=q10,Ex=Ex)
        # initialize to     
        self.cell = self.create_cell()    
        self.clamp = self.create_SEclamp(clamp_params)    
        self.set_nrn_params()
    # Create cell    
    def create_cell(self):
        #h("create soma")
        #cell = h.soma
        cell = h.Section(name='soma')
        cell.L = 4/np.pi # gives area of 4 um2
        cell.diam = 1
        cell.cm = 0
        cell.insert(self.chan_name)
        for ex_name,ex_val in self.params['Ex'].items():
            setattr(cell,ex_name,ex_val)
        setattr(cell,self.gbar_name,1)
        #cell.ena = self.params['Ena']
        #cell.ek = self.params['Ek']
        cell.nseg = 1
        cell.Ra = 1e22
#        self.cell = cell   
        return cell
    # Set up voltage clamp and recording
    def create_SEclamp(self,clamp_params):        
        if self.cell is None: self.create_cell() # create cell if doesn't exist 
        if not clamp_params: # empty, set to defaults
            clamp_params = {'v1':-120, 
                'dur1':200, 
                'v2':-80,
                'dur2':40,
                'v3': 10,
                'dur3': 20,
                'rs':1e-5,
                'x':0.5}
        clamp = h.SEClamp(self.cell(clamp_params['x']))
        clamp.rs = clamp_params['rs']
        clamp.dur1 = clamp_params['dur1']
        clamp.amp1 = clamp_params['amp1']
        clamp.dur2 = clamp_params['dur2']
        clamp.amp2 = clamp_params['amp2']
        clamp.dur3 = clamp_params['dur3']
        clamp.amp3 = clamp_params['amp3']
        self.params['tstop'] = clamp.dur1 + clamp.dur2 + clamp.dur3 # set tstop
        self.clamp = clamp
        return clamp # option to save externally 
    def set_nrn_params(self): # make sure to call before running (after create_SEclamp has been created)
        self.params
        h.dt = self.params['dt'] # set NEURON parameters
        h.celsius = self.params['T']
        h("q10 = {}".format(self.params['q10']))
        h.tstop = self.params['tstop']
        print("dt = {:.3f} ms, celsius = {}, q10 = {}, tstop = {}".format(h.dt,h.celsius,h.q10,self.params['tstop']))
    def create_recordings(self,cell):
        t_vec = h.Vector() # time vector
        t_vec.record(h._ref_t)        
        v_vec = h.Vector() # voltage vector
        v_vec.record(cell(0.5)._ref_v)
        i_vec = h.Vector() # current vector
        i_vec.record(getattr(cell(0.5),"_ref_{}".format(self.curr_name)))        
        g_vec = h.Vector() # conductance vector
        g_vec.record(getattr(cell(0.5),"_ref_{}".format(self.g_name)))
        return t_vec,v_vec,i_vec,g_vec
        

                

        
# Simulation loop

# Run activation protocol
def run_VClamp_act(vclamp_tester, clamp_params):        
    t_vec,v_vec,i_vec,g_vec = vclamp_tester.create_recordings(vclamp_tester.cell)
    
    v_steps = np.arange(clamp_params['v_stepi'], clamp_params['v_stepf'],clamp_params['dv_step'])
    vs= []
    currs = []
    gs = []    
    for vi in v_steps:        
        vclamp_tester.clamp.amp2 = vi # set new voltage step
        # run
        h.v_init = clamp_params['amp1']
        h.run()                        
        vs.append(v_vec.to_python(np.zeros(len(v_vec))))
        currs.append(i_vec.to_python(np.zeros(len(i_vec))))
        gs.append(g_vec.to_python(np.zeros(len(g_vec))))   
    t_vec = t_vec.to_python(np.zeros(len(t_vec)))
    return v_steps,t_vec,vs,currs,gs # outputs all numpy arrays
   



    
