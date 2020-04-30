import numpy as np 
import os
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
#os.system('nrnivmodl mechanisms')
#os.system('cp x86_64/special .')
# Set parameters
from setParams import setParams
from TMS_wave import TMS_wave
# Start NEURON
from neuron import gui
from neuron import h # make sure default -NSTACK and -NFRAME are set to 100000,20000
# Load cell information 
from cells import get_cell_files
from load_cell import create_cell
#*********************************#
# Get Params
params_vec = setParams(nrn_model_ver = "maxH")
hoc_dir = "nrn_stim_hoc"
cell_dir = "cells"
run_name = "test"
dt = params_vec[0]
tstop = params_vec[1]
# Generate TMS waveform
mode = 3 # MagProX100 Monophasic
delay = dt
plot_wave = 0
tvec, Evec = TMS_wave(dt,tstop,delay,mode,plot_wave)

# Initialize hoc code
h('objref tvec, Evec, params_vec')
h.tvec = h.Vector(tvec)
h.Evec = h.Vector(Evec)
h.params_vec = h.Vector(params_vec)
h.load_file(1,os.path.join(hoc_dir,"init_field_stim.hoc"))

# Load cell
cell_name = "L5_TTPC2_cADpyr232_1"
syn_enabled = False # synapses_enabled 
cell_files = get_cell_files.get_cell_files()[cell_name]
# create_cell
h.getParams(h.params_vec) # assign parameter values
h.replace_ax = 2
cell = create_cell(syn_enabled,cell_name,cell_files['morphology'],
                   cell_files['biophys_template'],cell_files['syn_template'])
t = h.Vector()
v_soma = h.Vector()
v_term = h.Vector()
t.record(h._ref_t)
v_soma.record(cell.soma[0](0.5)._ref_v)
v_term.record(h.Node[267](0.5)._ref_v) # for maxH
# Find threshold
h.getes()
h.init_record_spikes()

def stimul(plot_on=True,v_label="Node[373]",fig=None,ax=None):
    h.stimul()  
    #plt.cla()
   #ax.plot(t,v_soma,color='k',linewidth=1)
    if plot_on:
        if not fig:
            fig = plt.figure()
            ax = plt.subplot(111)
        ax.plot(t,v_term,linewidth=1,label=v_label)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xlabel('time (ms)'); ax.set_ylabel('Vm (mV)')
        plt.xlim([0, h.tstop]); plt.ylim([-90, 50])   
        if v_label: # not set to empty
            ax.legend()
        plt.show()
        return fig, ax

#h.find_thresh()
#threshE = h.threshE    
#h.AMP = threshE
#stimul()
#print "Threshold E =", threshE

Node_secList = h.Node_secList
axonal = h.cell.axonal
def change_mtau(mtau_scale=1,seclist=axonal):
    for sec in seclist:
        if h.ismembrane("NaTa_t"):
            sec.mtau_scale_NaTa_t = mtau_scale

def change_htau(htau_scale=1,seclist=axonal):
    for sec in seclist:
        if h.ismembrane("NaTa_t"):
            sec.htau_scale_NaTa_t = htau_scale


#plt.figure()
#ax = plt.subplot(111)
#ax.plot(mtau_scaling,threshEs)
