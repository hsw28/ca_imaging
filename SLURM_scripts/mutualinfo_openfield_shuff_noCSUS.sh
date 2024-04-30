module purge
module load matlab/r2022b
matlab -singleCompThread -batch "addpath(genpath('/home/hsw967/Programming/ca_imaging')); addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/')); submit_mutualinfo_openfield_shuff_noCSUS" -nodisplay -nosplash -nodesktop
