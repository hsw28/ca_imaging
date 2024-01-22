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

fprintf('loading traces')
allvariables = load('allvariables.mat');
calcium_traces = allvariables.Ca_traces;
clearvars allvariables

fprintf('loading pos')
pos = load('pos.mat');
pos_structure = pos.pos;

fprintf('loading MI')
MI = load('MI_CSUS.mat');
ca_MI = MI.MI_trace;

fprintf('loading timestamps')
Ca_ts = load('Ca_ts.mat')
ca_ts = Ca_ts.Ca_ts;


    f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 2, 2.5, ca_ts, 500, ca_MI);






end
