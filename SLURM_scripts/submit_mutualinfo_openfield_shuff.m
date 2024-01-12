% Get a handle to the cluster
c = parcluster;

%% Required arguments in order to submit MATLAB job

% Specify the walltime (e.g. 4 hours)
c.AdditionalProperties.WallTime = '12:00:00';

% Specify an account to use for MATLAB jobs (e.g. pXXXX, bXXXX, etc)
c.AdditionalProperties.AccountName = 'p32072';

% Specify a queue/partition to use for MATLAB jobs (e.g. short, normal, long)
c.AdditionalProperties.QueueName = 'normal';

%% optional arguments but worth considering

% Specify memory to use for MATLAB jobs, per core (default: 4gb)
c.AdditionalProperties.MemUsage = '16gb';

% Specify number of nodes to use
c.AdditionalProperties.Nodes = 1;

% Specify e-mail address to receive notifications about your job
c.AdditionalProperties.EmailAddress = 'hsw@northwestern.edu';

% The script that you want to run through SLURM needs to be in the MATLAB PATH
% Here we assume that quest_parallel_example.m lives in the same folder as submit_matlab_job.m

% Finally we will submit the MATLAB script quest_parallel_example to SLURM such that MATLAB
% will request enough resources to run a parallel pool of size 52 (i.e. parallelize over 52 CPUs).,
addpath(pwd)
addpath(genpath('/home/hsw967/Programming/ca_imaging'))
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/include'));
job = c.batch('mutualinfo_openfield_shuff_4SLURM', 'Pool', 4, 'CurrentFolder', '.');
