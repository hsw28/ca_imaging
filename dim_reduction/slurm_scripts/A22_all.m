maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(genpath('~/Documents/ca_imaging/dim_reduction/mind-paper-bb'));
addpath(genpath('~/Documents/ca_imaging/eyeblink'),'-end');
addpath(('~/Documents/ca_imaging'),'-end');

cd ~/Documents
load('EVERYTHING2.mat', 'frame_ts522')
load('EVERYTHING2.mat', 'Ca_traces')
load('EVERYTHING2.mat', 'times_US')
load('EVERYTHING2.mat', 'times_CS')
load('EVERYTHING2.mat', 'pos')
[wanted22 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_22, times_CS.CS_2023_05_22, times_US.US_2023_05_22, frame_ts522, pos.pos_2023_05_22, 1, 0);
wanted22 = wanted22(245:end);
time22 = vel(2,245:end);
moving22 = (Ca_traces.CA_traces_2023_05_22(:,wanted22));
result_A22_all = runMIND(moving22,time22)

save result_A22_all.mat result_A22_all
