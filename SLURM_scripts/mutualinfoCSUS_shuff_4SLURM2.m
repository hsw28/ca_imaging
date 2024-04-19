function mutualinfoCSUS_shuff_4SLURM2


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
fprintf('loading peaks')
spikes = slurm_var.Ca_peaks;
fprintf('loading CSUS_id')
CSUS_id = slurm_var.CSUS_id;



%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS AND NO PRETRIAL%%%%%%%
fprintf('loading MI1')
ca_MI0 = slurm_var.MI_CSUS2;
fprintf('loading MI2')
ca_MI1 = slurm_var.MI_CSUS2_pt;
clearvars slurm_var

%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS AND NO PRETRIAL%%%%%%%
fprintf('starting script no pretrial')
f = mutualinfo_CSUS_shuff2(spikes, CSUS_id, 0, 500, ca_MI0)


%%%%%BELOW IS FOR CSUS WITH 2 DIVISIONS AND PRETRIAL%%%%%%%
fprintf('starting script pretrial')
f = mutualinfo_CSUS_shuff2(spikes, CSUS_id, 1, 500, ca_MI1)



end
