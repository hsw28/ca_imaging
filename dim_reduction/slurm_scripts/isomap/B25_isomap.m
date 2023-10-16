maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(genpath('~/Documents/ca_imaging/dim_reduction/Isomap-master'));
addpath(genpath('~/Documents/ca_imaging/eyeblink'),'-end');
addpath(('~/Documents/ca_imaging'),'-end');

cd ~/Documents
load('EVERYTHING2.mat', 'frame_ts525')
load('EVERYTHING2.mat', 'frame_ts525')
load('EVERYTHING2.mat', 'frame_ts525')
load('EVERYTHING2.mat', 'frame_ts525')
load('EVERYTHING2.mat', 'Ca_traces')
load('EVERYTHING2.mat', 'times_US')
load('EVERYTHING2.mat', 'times_CS')
load('EVERYTHING2.mat', 'pos')
load('EVERYTHING2.mat', 'alignment_medium')

options.dims = 1:20;
D25_all = L2_distance(Ca_traces.CA_traces_2023_05_25, Ca_traces.CA_traces_2023_05_25);
[Y25_all, R25_all, E25_all] = IsomapII(D25_all, 'k', 20, options);


[wanted25 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_25, times_CS.CS_2023_05_25, times_US.US_2023_05_25, frame_ts525, pos.pos_2023_05_25_oval, 1, 0);
time25 = vel(2,:);
moving25 = (Ca_traces.CA_traces_2023_05_25(:,wanted25));
D25_moving = L2_distance(moving25, moving25);
[Y25_moving, R25_moving, E25_moving] = IsomapII(D25_moving, 'k', 20, options);

[training25 ts25] = converttoframe(times_CS.CS_2023_05_25, times_US.US_2023_05_25, frame_ts525);
wanted25 = find(training25>0);
ts25 = ts25(wanted25);
training25 = (Ca_traces.CA_traces_2023_05_25(:,wanted25));
D25_training = L2_distance(training25, training25);
[Y25_training, R25_training, E25_training] = IsomapII(D25_training, 'k', 20, options);


save result_B25_isomap.mat -v7.3
