function mutualinfoCSUS_trace_shuff_4SLURM8


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



fprintf('loading variables')
load('slurm_var.mat')
fprintf('loading traces')
calcium_traces = slurm_var.Ca_traces;
fprintf('loading IDs')
CSUS_id = slurm_var.CSUS_id;





%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS AND PRETRIAL%%%%%%%%
%ca_MI = slurm_var.MI_CSUS5_trace_pretrial;
%clearvars slurm_var
%f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 1, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 5 DIVISIONS NO PRETRIAL%%%%%%%%
%ca_MI = slurm_var.MI_CSUS5_trace;
%clearvars slurm_var
%f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 0, 5, 500, ca_MI)

%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS NO PRETRIAL%%%%%%%%
fprintf('loading MI')
ca_MI = slurm_var.MI_CSUS8_trace;
fprintf('clearing var')
clearvars slurm_var
fprintf('starting')
f = mutualinfo_CSUS_trace_shuff(calcium_traces, CSUS_id, 0, 500, ca_MI)


end
