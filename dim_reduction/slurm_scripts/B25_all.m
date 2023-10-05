maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(genpath('~/Documents/ca_imaging/dim_reduction/mind-paper-bb'));
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
[wanted25 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_25, times_CS.CS_2023_05_25, times_US.US_2023_05_25, frame_ts525, pos.pos_2023_05_25_oval, 1, 0);
time25 = vel(2,:);
moving25 = (Ca_traces.CA_traces_2023_05_25(:,wanted25));
result_B25_all = runMIND(moving25, time25);

save result_B25_all.mat result_B25_all
