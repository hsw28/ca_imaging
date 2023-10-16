maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(genpath('~/Documents/ca_imaging/dim_reduction/Isomap-master'));
addpath(genpath('~/Documents/ca_imaging/eyeblink'),'-end');
addpath(('~/Documents/ca_imaging'),'-end');

cd ~/Documents
load('EVERYTHING2.mat', 'frame_ts524')
load('EVERYTHING2.mat', 'frame_ts524')
load('EVERYTHING2.mat', 'frame_ts524')
load('EVERYTHING2.mat', 'frame_ts525')
load('EVERYTHING2.mat', 'Ca_traces')
load('EVERYTHING2.mat', 'times_US')
load('EVERYTHING2.mat', 'times_CS')
load('EVERYTHING2.mat', 'pos')
load('EVERYTHING2.mat', 'alignment_medium')

options.dims = 1:20;
D24_all = L2_distance(Ca_traces.CA_traces_2023_05_24, Ca_traces.CA_traces_2023_05_24);
[Y24_all, R24_all, E24_all] = IsomapII(D24_all, 'k', 20, options);


[wanted24 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_24, times_CS.CS_2023_05_24, times_US.US_2023_05_24, frame_ts524, pos.pos_2023_05_24_oval, 1, 0);
time24 = vel(2,:);
moving24 = (Ca_traces.CA_traces_2023_05_24(:,wanted24));
D24_moving = L2_distance(moving24, moving24);
[Y24_moving, R24_moving, E24_moving] = IsomapII(D24_moving, 'k', 20, options);

[training24 ts24] = converttoframe(times_CS.CS_2023_05_24, times_US.US_2023_05_24, frame_ts524);
wanted24 = find(training24>0);
ts24 = ts24(wanted24);
training24 = (Ca_traces.CA_traces_2023_05_24(:,wanted24));
D24_training = L2_distance(training24, training24);
[Y24_training, R24_training, E24_training] = IsomapII(D24_training, 'k', 20, options);


save result_B24_isomap.mat -v7.3
