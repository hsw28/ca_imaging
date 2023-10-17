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

timestamps = frame_ts524;
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

time24 = timestamps(tsindex);

%b1 and b2, sep and combined
vB24B25 = find(alignment_medium(:,17)>0 & alignment_medium(:,18)>0);
vB24B25_24 = alignment_medium(vB24B25,17);
moving24 = (Ca_traces.CA_traces_2023_05_24(vB24B25_24,:));
result_B24B25_24 = runMIND(moving24, time24)

%B24B25 = find(alignment_medium(:,17)>0 & alignment_medium(:,18)>0);
%B24B25_25 = alignment_medium(B24B25,18);
%[wanted25 temp temp1 vel post] = movingtimetraining(Ca_traces.CA_traces_2023_05_25, times_CS.CS_2023_05_25, times_US.US_2023_05_25, frame_ts525, pos.pos_2023_05_25_oval, 1, 0);
%time25 = vel(2,:);
%moving25 = (Ca_traces.CA_traces_2023_05_25(B24B25_25,wanted25));
%result_B24B25_25 = runMIND(moving25, time25)

%moving24_25 = [moving24 moving25];
%time24_25 = [time24 time25];
%result_B24_B25_concat = runMIND(moving24_25, time24_25);
%
save result_B24B25_24_alltimes.mat result_B24B25_24
