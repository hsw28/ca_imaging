module purge
module load matlab/r2018b
matlab -singleCompThread -batch "addpath(genpath('/home/hsw967/Programming/ca_imaging')); addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/')); submit_mutualinfo_openfield_shuff" -nodisplay -nosplash -nodesktop
