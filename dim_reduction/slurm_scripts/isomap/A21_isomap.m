maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(genpath('~/Documents/ca_imaging/dim_reduction/Isomap-master'));
addpath(genpath('~/Documents/ca_imaging/eyeblink'),'-end');
addpath(('~/Documents/ca_imaging'),'-end');

cd ~/Documents
load('EVERYTHING2.mat', 'frame_ts521')
load('EVERYTHING2.mat', 'frame_ts522')
load('EVERYTHING2.mat', 'frame_ts524')
load('EVERYTHING2.mat', 'frame_ts525')
load('EVERYTHING2.mat', 'Ca_traces')
load('EVERYTHING2.mat', 'times_US')
load('EVERYTHING2.mat', 'times_CS')
load('EVERYTHING2.mat', 'pos')
load('EVERYTHING2.mat', 'alignment_medium')

options.dims = 1:20;
D21_all = L2_distance(Ca_traces.CA_traces_2023_05_21, Ca_traces.CA_traces_2023_05_21);
[Y21_all, R21_all, E21_all] = IsomapII(D21_all, 'k', 20, options);


[wanted21 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_21, times_CS.CS_2023_05_21, times_US.US_2023_05_21, frame_ts521, pos.pos_2023_05_21, 1, 0);
time22 = vel(2,:);
moving22 = (Ca_traces.CA_traces_2023_05_22(:,wanted21));
D21_moving = L2_distance(moving22, moving22);
[Y21_moving, R21_moving, E21_moving] = IsomapII(D21_moving, 'k', 20, options);

[training21 ts21] = converttoframe(times_CS.CS_2023_05_21, times_US.US_2023_05_21, frame_ts521);
wanted21 = find(training21>0);
ts21 = ts21(wanted21);
training21 = (Ca_traces.CA_traces_2023_05_21(:,wanted21));
D21_training = L2_distance(training21, training21);
[Y21_training, R21_training, E21_training] = IsomapII(D21_training, 'k', 20, options);


save result_A21_isomap.mat
