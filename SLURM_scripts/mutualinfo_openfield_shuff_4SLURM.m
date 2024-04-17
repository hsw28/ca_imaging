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
load('slurm_var.mat');
fprintf('loading pos')
pos_structure = slurm_var.pos;
fprintf('loading spikes')
spikes = slurm_var.Ca_peaks;
fprintf('loading TS')
ca_ts = slurm_var.Ca_ts;
fprintf('done loading')


fprintf('loading MI')
%for regular MI
%ca_MI = slurm_var.MI;
%for MI 5
ca_MI = slurm_var.MI5;


clearvars slurm_var
fprintf('cleared excess variables, starting process')
%for regular
%f = mutualinfo_openfield_shuff(spikes, pos_structure, 4, 2.5, ca_ts, 500, ca_MI);
%for MI5

f = mutualinfo_openfield_shuff(spikes, pos_structure, 4, 5, ca_ts, 500, ca_MI);

end
