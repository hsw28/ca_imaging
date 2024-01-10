function mutualinfo_openfield_shuff_4SLURM


  %maxNumCompThreads(str2num(getenv('SLURM_NPROCS'))); %%%%%trying without this



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

fprintf('starting4')
allvariables = load('allvariables.mat');
pos_structure = allvariables.pos;
spikes = allvariables.peaks;
ca_MI = allvariables.MI;

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 3, ca_MI)
%f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 500, ca_MI)
MI_shuff = f;
% Save the output to a .mat file
save('mutualinfo_shuff_output.mat', 'MI_shuff');
end
