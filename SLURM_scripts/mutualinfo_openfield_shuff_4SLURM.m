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

fprintf('loading pos')
pos = load('pos.mat');
pos_structure = pos.pos;

fprintf('loading peaks')
peaks = load('peaks.mat');
spikes = peaks.peaks;

fprintf('loading MI')
MI = load('MI.mat');
ca_MI = MI.MI;

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 500, ca_MI)
MI_shuff = f;

end
