module purge
module load matlab/r2018b
matlab -singleCompThread -batch "addpath(genpath('/home/hsw967/Programming/ca_imaging')); submit_mutualinfo_openfield_trace_shuff"
