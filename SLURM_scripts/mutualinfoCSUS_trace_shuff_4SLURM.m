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


CSUS = load('CSUS_id.mat')
CSUS_id = CSUS.CSUS_id;

allvariables = load('allvariables.mat');
calcium_traces = allvariables.Ca_traces;
clearvars allvariables

MI = load('MI_CSUS.mat');



%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS AND PRETRIAL%%%%%%%%
%ca_MI = MI.MI_CSUS5_trace_pretrial;
%f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 1, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS NO PRETRIAL%%%%%%%%
ca_MI = MI.MI_CSUS5_trace;
f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 0, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS NO PRETRIAL%%%%%%%%
ca_MI = MI.MI_CSUS2_trace;
f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 0, 2, 500, ca_MI)


end
