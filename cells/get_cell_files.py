# -*- coding: utf-8 -*-
import os

# Makes dictionary with morphology & biophysics file names for each cell
temp_name_line_morph = 34 # 34th line in morphology.hoc
temp_name_word_morph = 2 # 2nd word in morphology.hoc
temp_name_line_biophys = 33 # 33rd line in biophysics.hoc
temp_name_word_biophys = 2 # 2nd word in biophysics.hoc
temp_name_line_syn = 31 # 31st line in synapses.hoc
temp_name_word_syn = 2 # 2nd word in synapses.hoc

def get_cell_files():
    cell_dir = os.path.dirname(os.path.realpath(__file__))
    cell_files = {} 
    for celli in os.listdir(cell_dir):
        celli_path = os.path.join(cell_dir,celli)
        if os.path.isdir(celli_path): # make sure its a folder            
             # morphology template name in morphology.hoc             
            with open(os.path.join(celli_path,"morphology.hoc")) as morph_hoc:
                morph_tempi = morph_hoc.readlines()[temp_name_line_morph-1].split()[temp_name_word_morph-1]
            # biophysics template name in biophysics.hoc            
            with open(os.path.join(celli_path,"biophysics.hoc")) as biophys_hoc: 
                biophys_tempi = biophys_hoc.readlines()[temp_name_line_biophys-1].split()[temp_name_word_biophys-1] 
            # synapses template name in biophysics.hoc            
            with open(os.path.join(celli_path,"synapses.hoc")) as syn_hoc: 
                syn_tempi = syn_hoc.readlines()[temp_name_line_syn-1].split()[temp_name_word_syn-1]                
             # morphology asc file
            files = os.listdir(celli_path) 
            morph_asci = next(f for f in files if f.startswith('dend'))
            # append celli to cell_files 
            cell_files[celli] = {'morphology': morph_asci, # asc file    
                'morph_template': morph_tempi,'biophys_template': biophys_tempi,
                'syn_template' : syn_tempi} 
    return cell_files
    

#cell_files = get_cell_files()