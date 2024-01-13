function mutualinfoCSUS_trace_shuff_4SLURM


    addpath(pwd);
    addpath(genpath('/home/hsw967/Programming/ca_imaging'));
    addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));
    addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/include'));





%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos


peaks = load('peaks.mat');
spikes = peaks.peaks;

MI = load('MI_CSUS.mat');
ca_MI = MI.MI_CSUS2_trace;

CSUS_id = load('CSUS_id.mat')
CSUS_id = CSUS_id.CSUS_id;

allvariables = load('allvariables.mat');
calcium_traces = allvariables.all_traces;
clearvars allvariables

f = mutualinfo_CSUS_shuff(calcium_traces, CSUS_id, 1, 2, 500, ca_MI)

MI_shuff = f;

% Get the current date and time as a string
currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
% Create a filename with the timestamp
filename = ['mutualinfo_shuff_output_', currentDateTime, '.mat'];
% Save the output to the .mat file with the timestamped filename
save(filename, 'MI_shuff');


end
