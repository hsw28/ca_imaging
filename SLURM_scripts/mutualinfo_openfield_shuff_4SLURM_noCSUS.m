function mutualinfo_openfield_shuff_4SLURM_noCSUS


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
fprintf('loading MI')
ca_MI = slurm_var.MI_noCSUS;
fprintf('loading CSU_id')
CSUS_id = slurm_var.CSUS_id;
fprintf('done loading')


clearvars slurm_var
fprintf('cleared excess variables, starting process')


f = mutualinfo_openfield_shuff_noCSUS(spikes, pos_structure, 4, 2.5, ca_ts, 500, ca_MI, CSUS_id)

end
