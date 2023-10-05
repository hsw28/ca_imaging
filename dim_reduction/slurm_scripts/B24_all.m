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
[wanted24 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_24, times_CS.CS_2023_05_24, times_US.US_2023_05_24, frame_ts524, pos.pos_2023_05_24_oval, 1, 0);
time24 = vel(2,:);
moving24 = (Ca_traces.CA_traces_2023_05_24(:,wanted24));
result_B24_all = runMIND(moving24, time24);

save result_B24_all.mat result_B24_all
