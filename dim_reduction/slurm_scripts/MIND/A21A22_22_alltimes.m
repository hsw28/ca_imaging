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

%an-1 and an, sep and combined

%A21A22 = find(alignment_medium(:,15)>0 & alignment_medium(:,14)>0);
%A21A22_21 = alignment_medium(A21A22,14);
%[wanted22 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_21, times_CS.CS_2023_05_21, times_US.US_2023_05_21, frame_ts521, pos.pos_2023_05_21, 1, 0);
%time21 = vel(2,:);
%moving21 = (Ca_traces.CA_traces_2023_05_21(A21A22_21,wanted21));
%result_AA21A22_21 = runMIND(moving21, time21);


vA21A22 = find(alignment_medium(:,15)>0 & alignment_medium(:,14)>0);
vA21A22_22 = alignment_medium(vA21A22,15);
moving22 = (Ca_traces.CA_traces_2023_05_22(vA21A22_22,647:end));
result_A21A22_22 = runMIND(moving22, time22);

%moving21_22 = [moving21 moving22];
%time21_22 = [time21 time22];
%result_A21A22_concat = runMIND(moving21_22, time21_22);


save result_A21A22_22_alltimes.mat result_A21A22_22
