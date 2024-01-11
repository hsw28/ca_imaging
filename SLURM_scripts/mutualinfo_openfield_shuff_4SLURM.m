function mutualinfo_openfield_shuff_4SLURM
  maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));


    addpath(pwd);
    addpath(genpath('/home/hsw967/Programming/ca_imaging'));
    addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));

c = parcluster;
c.AdditionalProperties.WallTime = '03:00:00';
c.AdditionalProperties.AccountName = 'p32072';
c.AdditionalProperties.QueueName = 'short';

%  pool = c.parpool(8);

%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos
pos = load('pos.mat');
pos_structure = pos.pos;

peaks = load('pos.mat');
spikes = peaks.peaks;

MI = load('MI.mat');
ca_MI = MI.MI;

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 3, ca_MI)
%f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 500, ca_MI)
MI_shuff = f;
% Save the output to a .mat file
save('mutualinfo_shuff_output.mat', 'MI_shuff');

end
