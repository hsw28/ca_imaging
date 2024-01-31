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



CSUS = load('CSUS_id.mat')
CSUS_id = CSUS.CSUS_id;

MI = load('MI_CSUS.mat');




%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS AND PRETRIAL%%%%%%%%
%ca_MI = MI.MI_CSUS5_pretrial;
%f = mutualinfo_CSUS_shuff(spikes, CSUS_id, 1, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS AND NO PRETRIAL%%%%%%%%
ca_MI = MI.MI_CSUS5;
f = mutualinfo_CSUS_shuff(spikes, CSUS_id, 0, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS AND NO PRETRIAL%%%%%%%
ca_MI = MI.MI_CSUS2;
f = mutualinfo_CSUS_shuff(spikes, CSUS_id, 0, 2, 500, ca_MI)



end
