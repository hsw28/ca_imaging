function mutualinfo_openfield_trace_shuff_4SLURM_noCSUS

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
fprintf('loading traces')
calcium_traces = slurm_var.Ca_traces;
fprintf('loading times')
ca_ts = slurm_var.Ca_ts;
fprintf('loading MI')
ca_MI = slurm_var.MI_trace_noCSUS;
fprintf('loading CSU_id')
CSUS_id = slurm_var.CSUS_id;
fprintf('done loading')

fprintf('clearing excess')
clearvars slurm_var
fprintf('starting script')
%regular
%f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 4, 2.5, ca_ts, 500, ca_MI);

%5
f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 4, 2.5, ca_ts, 500, ca_MI);
f = mutualinfo_openfield_trace_shuff_noCSUS(calcium_traces, pos_structure, 4, 2.5, ca_ts, 500, ca_MI, CSUS_id)
%finds mutual info for a bunch of cells




end
