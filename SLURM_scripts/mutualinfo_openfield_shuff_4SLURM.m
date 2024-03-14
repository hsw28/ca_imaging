function mutualinfo_openfield_shuff_4SLURM


    addpath(pwd);
    addpath(genpath('/home/hsw967/Programming/ca_imaging'));
    addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));
    addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/include'));





%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos
  %Ca_ts

fprintf('loading variables')
slurm_var = load('slurm_var.mat')
pos_structure = slurm_var.pos;
spikes = slurm_var.Ca_peaks;
ca_MI = slurm_var.MI;
ca_ts = slurm_var.Ca_ts;
clearvars slurm_var

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, ca_ts, 500, ca_MI);

end
