function move_cell_files(cell_ids)
if nargin == 0
   cell_ids = [17:21,32:36];  
end
cell_dir = '/Users/amanaberra/Documents/NEURON/Parallel_plate_TAL/cells/blue_brain_cells'; 
for i = cell_ids
   cell_name = cellmodelnames(i); 
   if exist(cell_name,'dir') == 0
      mkdir(cell_name);  
   end
   copyfile(fullfile(cell_dir,cell_name,'biophysics.hoc'),cell_name);
   copyfile(fullfile(cell_dir,cell_name,'morphology/*'),cell_name);
   copyfile(fullfile(cell_dir,cell_name,'morphology.hoc'),cell_name);
   copyfile(fullfile(cell_dir,cell_name,'synapses/synapses.hoc'),cell_name);
end
