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
%ca_MI = MI.MI_CSUS2_trace_pretrial;
%ca_MI = MI.MI_CSUS5_trace_pretrial;

CSUS = load('CSUS_id.mat')
CSUS_id = CSUS.CSUS_id;

allvariables = load('allvariables.mat');
calcium_traces = allvariables.Ca_traces;
clearvars allvariables

MI = load('MI_CSUS.mat');

%finds 'mutual info' for CS/US/ non CS/US
%CSUS_structure should come from BULKconverttoframe.m
%do_you_want_CSUS_or_CSUSnone: 1 for only cs us, 0 for cs us none
%how many divisions you wanted-- for ex,
    % do_you_want_CSUS_or_CSUSnone = 1
    % how_many_divisions = 2 will just split between cs and us
                        %= 10 will split CS and US each into 5
%right now because im lazy how_many_divisions must be a factor of 10




%ca_MI = MI.MI_CSUS2_trace;
%ca_MI = MI.MI_CSUS2_trace_pretrial;

%f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 0, 2, 500, ca_MI)


%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS AND PRETRIAL%%%%%%%%
%ca_MI = MI.MI_CSUS5_trace_pretrial;
%f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 0, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS NO PRETRIAL%%%%%%%%
ca_MI = MI.MI_CSUS5_trace;
f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 1, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS NO PRETRIAL%%%%%%%%
ca_MI = MI.MI_CSUS2_trace;
f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 1, 2, 500, ca_MI)


end
