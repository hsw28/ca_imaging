% Get a handle to the cluster
c = parcluster;

% Specify the walltime (e.g. 24 hours)
c.AdditionalProperties.WallTime = '24:00:00';
% Specify an account to use for MATLAB jobs (e.g. pXXXX, bXXXX, etc)
c.AdditionalProperties.AccountName = 'p32072';
% Specify a queue/partition to use for MATLAB jobs (e.g. short, normal, long)
c.AdditionalProperties.QueueName = 'normal';
% Specify memory to use for MATLAB jobs
c.AdditionalProperties.MemUsage = '64gb';
% Specify number of nodes to use
c.AdditionalProperties.Nodes = 1;
% Specify number of CPUs per node (optional, depending on your job requirements)
c.AdditionalProperties.CpusPerNode = 8;
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
job = c.batch('trace_shuff_4SLURM', 'CurrentFolder', '.');
