function mutualinfoCSUS_shuff_4SLURM


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


peaks = load('peaks.mat');
spikes = peaks.peaks;

MI = load('MI_CSUS.mat');
ca_MI = MI.MI_CSUS2;
ca_MI = MI.MI_CSUS2_pretrial;

CSUS_id = load('CSUS_id.mat')
CSUS_id = CSUS_id.CSUS_id;



f = mutualinfo_CSUS_shuff(spikes, CSUS_id, 0, 2, 500, ca_MI)

MI_shuff = f;


end
