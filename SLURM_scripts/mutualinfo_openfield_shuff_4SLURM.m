function mutualinfo_openfield_shuff_4SLURM
  maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));


    addpath(pwd);
    addpath(genpath('/home/hsw967/Programming/ca_imaging'));
    addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));



    % open the parallel pool, recording the time it takes
  %  p = parcluster('local');
  %  tic;
  %  parpool(p, 52); % open the pool using 28 workers
  %  fprintf('Opening the parallel pool took %g seconds.\n', toc)




%    c = parcluster;
%    c.AdditionalProperties.WallTime = '03:00:00';
%    c.AdditionalProperties.AccountName = 'p32072';
%    c.AdditionalProperties.QueueName = 'short';
%    c.parpool(4)
%    tic;
%    fprintf('Opening the parallel pool took %g seconds.\n', toc)


%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos
pos = load('pos.mat');
pos_structure = pos.pos;

peaks = load('peaks.mat');
spikes = peaks.peaks;

MI = load('MI.mat');
ca_MI = MI.MI;

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 500, ca_MI)
MI_shuff = f;

% Get the current date and time as a string
currentDateTime = datestr(now, 'yyyymmdd_HHMMSS');
% Create a filename with the timestamp
filename = ['mutualinfo_shuff_output_', currentDateTime, '.mat'];
% Save the output to the .mat file with the timestamped filename
save(filename, 'MI_shuff');


end
