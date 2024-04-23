function mutualinfoCSUS_shuff_4SLURM5


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
spikes = slurm_var.Ca_peaks;
CSUS_id = slurm_var.CSUS_id10;




ca_MI = slurm_var.MI_CSUS5;
clearvars slurm_var
f = mutualinfo_CSUS_shuff5(spikes, CSUS_id, 0, 500, ca_MI)



end
