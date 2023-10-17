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

timestamps = frame_ts522;
if isa(timestamps,'table')
  timestamps = table2array(timestamps);
  timestamps = timestamps(:,2);
end

if size(timestamps,2)==3
  timestamps = timestamps(:,2);
end

if timestamps(5)>2
timestamps = timestamps./1000;
end
tsindex = 2:2:length(timestamps);

time22 = timestamps(tsindex);
time22 = time22(647:end);

%an and b1, sep and combined
vA22B24 = find(alignment_medium(:,15)>0 & alignment_medium(:,17)>0);
vA22B24_22 = alignment_medium(vA22B24,15);
moving22 = (Ca_traces.CA_traces_2023_05_22(vA22B24_22,647:end));
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
