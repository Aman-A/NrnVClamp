#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 21:07:12 2018

@author: amanaberra
"""
import os
from neuron import h
cell_dir = "cells"
hoc_dir = "nrn_stim_hoc"
def create_cell(syn_enabled,cell_name,morph_file,biophys_temp,syn_temp):
    prune_meth1_cells = ['L23_PC_cADpyr229_3','L5_TTPC2_cADpyr232_5']
    if cell_name in prune_meth1_cells:
        h("prune_meth = 1") # define hoc variable
    else: 
        h("prune_meth = 2")
    h.load_file("import3d.hoc")
    h.load_file(os.path.join(cell_dir,cell_name,"biophysics.hoc")) # cell specific biophys template
    h.load_file(os.path.join(cell_dir,cell_name,"synapses.hoc")) # cell specific biophys template    
    h.load_file(os.path.join(hoc_dir,"morphology_temp.hoc")) # generic morph template
    h.load_file(os.path.join(hoc_dir,"template.hoc")) # generic Blue Brain cell template
    h("objref cell")    
    h('cell = new BlueBrainCell({},"{}","{}","{}","{}")'.format(1 if syn_enabled else 0,cell_name,
                                morph_file,biophys_temp,syn_temp))
    h("modify_cell(cell,prune_meth)")
    cell = h.cell
    print "Loaded cell", cell_name
    return cell