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
load('EVERYTHING2.mat', 'alignment_medium')

%an and b1, sep and combined
vA22B24 = find(alignment_medium(:,15)>0 & alignment_medium(:,17)>0);
vA22B24_22 = alignment_medium(vA22B24,15);
[wanted22 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_22, times_CS.CS_2023_05_22, times_US.US_2023_05_22, frame_ts522, pos.pos_2023_05_22, 1, 0);
wanted22 = wanted22(245:end);
time22 = vel(2,245:end);
moving22 = (Ca_traces.CA_traces_2023_05_22(vA22B24_22,wanted22));
result_AA22B24_22 = runMIND(moving22, time22);

%
%A22B24 = find(alignment_medium(:,15)>0 & alignment_medium(:,17)>0);
%A22B24_24 = alignment_medium(A22B24,17);
%[wanted24 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_24, times_CS.CS_2023_05_24, times_US.US_2023_05_24, frame_ts524, pos.pos_2023_05_24_oval, 1, 0);
%time24 = vel(2,:);
%moving24 = (Ca_traces.CA_traces_2023_05_24(A22B24_24,wanted24));
%result_A22B24_24 = runMIND(moving24, time24);
%
%moving22_24 = [moving22 moving24];
%time22_24 = [time22 time24];
%result_A22B24_concat = runMIND(moving22_24, time22_24);
%

save result_AA22B24_22.mat result_AA22B24_22
