function mutualinfoCSUS_trace_shuff_4SLURM


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


  fprintf('loading traces')
  allvariables = load('allvariables.mat');
  calcium_traces = allvariables.Ca_traces;
  clearvars allvariables

MI = load('MI_CSUS.mat');
%ca_MI = MI.MI_CSUS2_trace;
ca_MI = MI.MI_CSUS2_trace_pretrial;
ca_MI = MI.MI_CSUS5_trace_pretrial;

CSUS = load('CSUS_id.mat')
CSUS_id = CSUS.CSUS_id;

allvariables = load('allvariables.mat');
calcium_traces = allvariables.Ca_traces;
clearvars allvariables

f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 0, 2, 500, ca_MI)




end
