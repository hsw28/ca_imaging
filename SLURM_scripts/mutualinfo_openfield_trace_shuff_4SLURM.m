
function mutualinfo_openfield_trace_shuff_4SLURM
% Refresh MATLAB's toolbox cache
rehash toolboxcache;

% Custom configuration for the cluster (if required)
configCluster;

c = parcluster;
c.AdditionalProperties.AccountName = 'p32072'; % Replace with your actual account name
c.AdditionalProperties.WallTime = '24:00:00'; % Set to match the wall time in your SLURM script
c.AdditionalProperties.QueueName = 'normal';
c.AdditionalProperties.MemUsage = '64gb';

c.saveProfile;

addpath(pwd);
addpath(genpath('/home/hsw967/Programming/ca_imaging'));
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));

% Set the number of computational threads to the number of allocated CPUs
if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
    maxNumCompThreads(str2num(getenv('SLURM_CPUS_PER_TASK')));
end

pool = c.parpool(8);

pos_structure = load('pos.mat');
calcium_traces = load('traces.mat')
ca_MI = load('MI_trace.mat')

f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 2, 2.5, 500, ca_MI)
% Save the output to a .mat file
MI_trace_shuff = f;
save('mutualinfo_trace_shuff_output.mat', 'MI_trace_shuff');
end
