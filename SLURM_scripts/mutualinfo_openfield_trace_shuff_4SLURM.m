% Get a handle to the cluster
c = parcluster;

%% Required arguments in order to submit a MATLAB GPU job
% Specify the walltime (e.g. 4 hours)
c.AdditionalProperties.WallTime = '024:00:00';
% Specify an account to use for MATLAB jobs (e.g. pXXXX, bXXXX, etc)
c.AdditionalProperties.AccountName = 'p32072';
% Specify a queue/partition to use for MATLAB jobs (e.g. short, normal, long)
c.AdditionalProperties.QueueName = 'normal';
% Specify number of GPUs
c.AdditionalProperties.GpusPerNode = 1;
% Specify type of GPU card to use (e.g. a100)
c.AdditionalProperties.GpuCard = 'a100';
%% optional arguments but worth considering
% Specify memory to use for MATLAB jobs, per core (default: 4gb)
c.AdditionalProperties.MemUsage = '64gb';
% Specify number of nodes to use
c.AdditionalProperties.Nodes = 1;
% Specify e-mail address to receive notifications about your job
c.AdditionalProperties.EmailAddress = 'hsw967@northwestern.edu';

rehash toolboxcache
configCluster

% The script that you want to run through SLURM needs to be in the MATLAB PATH
% Here we assume that quest_gpu_example.m lives in the same folder as submit_matlab_job.m
addpath(pwd);
addpath(genpath('/home/hsw967/Programming/ca_imaging'));
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));


% Finally we will submit the MATLAB script quest_gpu_example to SLURM such that MATLAB
job = c.batch('mutualinfo_openfield_trace_shuff_4SLURM', 'CurrentFolder', '.');


function mutualinfo_openfield_trace_shuff_4SLURM
  % Refresh MATLAB's toolbox cache
rehash toolboxcache;

% Custom configuration for the cluster (if required)
% configCluster;

% Set the number of computational threads to the number of allocated CPUs
if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
    maxNumCompThreads(str2num(getenv('SLURM_CPUS_PER_TASK')));
end

% Start a parallel pool with the specified number of workers
numWorkers = 8;
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool(numWorkers);
end

pos_structure = load('pos.mat');
calcium_traces = load('traces.mat')
ca_MI = load('MI_trace.mat')

f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 2, 2.5, 500, ca_MI)
% Save the output to a .mat file
MI_trace_shuff = f;
save('mutualinfo_trace_shuff_output.mat', 'MI_trace_shuff');
end
