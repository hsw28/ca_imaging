function mutualinfoCSUS_trace_shuff_4SLURM2_10


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


%for 10 divisions
CSUS_id = slurm_var.CSUS_id10;
ca_MI1 = slurm_var.MI_CSUS2_trace_10;
ca_MI2 = slurm_var.MI_CSUS2_trace_pt_10;


%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS NO PRETRIAL%%%%%%%%
fprintf('loading MI')
f = mutualinfo_CSUS_trace_shuff2(calcium_traces, CSUS_id, 0, 500, ca_MI1)


f = mutualinfo_CSUS_trace_shuff2(calcium_traces, CSUS_id, 1, 500, ca_MI2)


end
