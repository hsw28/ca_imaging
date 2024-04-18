function mutualinfo_openfield_trace_shuff_4SLURM

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
%for regular MI
ca_MI = slurm_var.MI_trace;
%for MI5
%ca_MI = slurm_var.MI_trace5;

fprintf('clearing excess')
clearvars slurm_var
fprintf('starting script')
%regular
%f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 4, 2.5, ca_ts, 500, ca_MI);

%5
f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 4, 5, ca_ts, 500, ca_MI);




end
