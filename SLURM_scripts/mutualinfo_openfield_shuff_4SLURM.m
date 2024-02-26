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

fprintf('loading pos')
pos = load('pos.mat');
pos_structure = pos.pos;

fprintf('loading peaks')
peaks = load('peaks.mat');
spikes = peaks.peaks;

fprintf('loading MI')
MIs = load('MI.mat');
ca_MI = MIs.MI;

fprintf('loading timestamps')
Ca_ts = load('Ca_ts.mat')
ca_ts = Ca_ts.Ca_ts;

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, ca_ts, 500, ca_MI);

end
