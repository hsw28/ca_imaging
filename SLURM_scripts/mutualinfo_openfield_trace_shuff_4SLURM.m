function mutualinfo_openfield_trace_shuff_4SLURM
  maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

addpath(pwd);
addpath(genpath('/home/hsw967/Programming/ca_imaging'));
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));


%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos
  %Ca_ts

allvariables = load('allvariables.mat');
calcium_traces = allvariables.all_traces;

pos = load('pos.mat');
pos_structure = pos.pos;

MI_trace = load('MI_trace.mat');
ca_MI = MI_trace.MI_trace;

Ca_ts = load('Ca_ts.mat')
ca_ts = Ca_ts.Ca_ts;

f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 2, 2.5, 500, ca_MI, ca_ts);
% Save the output to a .mat file
MI_trace_shuff = f;


% Get the current date and time as a string
currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
% Create a filename with the timestamp
filename = ['mutualinfo_trace_shuff_output_', currentDateTime, '.mat'];
% Save the output to the .mat file with the timestamped filename
save(filename, 'MI_trace_shuff');

end
