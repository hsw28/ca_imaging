function mutualinfoCSUS_trace_shuff_4SLURM2


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





%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS NO PRETRIAL%%%%%%%%
fprintf('loading MI')
ca_MI = slurm_var.MI_CSUS2_trace;
f = mutualinfo_CSUS_trace_shuff2(calcium_traces, CSUS_id, 0, 500, ca_MI)

ca_MI = slurm_var.MI_CSUS2_trace_pt;
f = mutualinfo_CSUS_trace_shuff2(calcium_traces, CSUS_id, 1, 500, ca_MI)


end
