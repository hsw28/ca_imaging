module purge
module load matlab/r2018b
matlab -singleCompThread -batch "addpath('/home/hsw967/Programming/ca_imaging/SLURM_scripts'); submit_mutualinfo_openfield_shuff"
