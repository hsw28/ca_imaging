function mutualinfo_openfield_shuff_4SLURM

  maxNumCompThreads(52);
  cd ('/projects/p32072/data_eyeblink/rat314/slurm_variables')
  outputFile = fopen('/projects/p32072/data_eyeblink/rat314/slurm_variables/mutualinfo_shuff_output.log', 'w');
  try


    addpath(pwd);
    addpath(genpath('/home/hsw967/Programming/ca_imaging'));
    addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));


%  pool = c.parpool(8);

%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos

parpool('p32072', 12)
allvariables = load('/projects/p32072/data_eyeblink/rat314/slurm_variables/allvariables.mat');
pos_structure = allvariables.pos;
spikes = allvariables.peaks;
ca_MI = allvariables.MI;

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 5, ca_MI)
MI_shuff = f;
% Save the output to a .mat file


try
    save('./mutualinfo_shuff_output.mat', 'MI_shuff');
catch e
    disp('Error occurred while saving the file:');
    disp(e.message);
end


fprintf(outputFile, 'Computation completed successfully.\n');
catch e
fprintf(outputFile, 'An error occurred: %s\n', e.message);
end
 fclose(outputFile);
end
