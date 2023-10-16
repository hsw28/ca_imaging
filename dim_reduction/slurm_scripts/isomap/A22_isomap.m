maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(genpath('~/Documents/ca_imaging/dim_reduction/Isomap-master'));
addpath(genpath('~/Documents/ca_imaging/eyeblink'),'-end');
addpath(('~/Documents/ca_imaging'),'-end');

cd ~/Documents
load('EVERYTHING2.mat', 'frame_ts522')
load('EVERYTHING2.mat', 'frame_ts522')
load('EVERYTHING2.mat', 'frame_ts524')
load('EVERYTHING2.mat', 'frame_ts525')
load('EVERYTHING2.mat', 'Ca_traces')
load('EVERYTHING2.mat', 'times_US')
load('EVERYTHING2.mat', 'times_CS')
load('EVERYTHING2.mat', 'pos')
load('EVERYTHING2.mat', 'alignment_medium')

options.dims = 1:20;
D22_all = L2_distance(Ca_traces.CA_traces_2023_05_22(:,647:end), Ca_traces.CA_traces_2023_05_22(:,647:end));
[Y22_all, R22_all, E22_all] = IsomapII(D22_all, 'k', 20, options);


[wanted22 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_22, times_CS.CS_2023_05_22, times_US.US_2023_05_22, frame_ts522, pos.pos_2023_05_22, 1, 0);
wanted22 = wanted22(245:end);
time22 = vel(2,245:end);
moving22 = (Ca_traces.CA_traces_2023_05_22(:,wanted22));
D22_moving = L2_distance(moving22, moving22);
[Y22_moving, R22_moving, E22_moving] = IsomapII(D22_moving, 'k', 20, options);

[training22 ts22] = converttoframe(times_CS.CS_2023_05_22, times_US.US_2023_05_22, frame_ts522);
wanted22 = find(training22>0);
wanted22 = wanted22(19:end);
ts22 = ts22(wanted22);
training22 = (Ca_traces.CA_traces_2023_05_22(:,wanted22));
D22_training = L2_distance(training22, training22);
[Y22_training, R22_training, E22_training] = IsomapII(D22_training, 'k', 20, options);


save result_A22_isomap.mat
