maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(genpath('~/Documents/ca_imaging/dim_reduction/mind-paper-bb'));
addpath(genpath('~/Documents/ca_imaging/eyeblink'),'-end');
addpath(('~/Documents/ca_imaging'),'-end');

cd ~/Documents
load('EVERYTHING2.mat', 'frame_ts521')
load('EVERYTHING2.mat', 'Ca_traces')
load('EVERYTHING2.mat', 'times_US')
load('EVERYTHING2.mat', 'times_CS')
load('EVERYTHING2.mat', 'pos')
[wanted22 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_21, times_CS.CS_2023_05_21, times_US.US_2023_05_21, frame_ts521, pos.pos_2023_05_21, 1, 0);
time21 = vel(2,:);
moving21 = (Ca_traces.CA_traces_2023_05_21(:,wanted21));
result_A21_all = runMIND(moving21, time21);
save result_A21_all.mat result_A21_all
