c = parcluster;

% Specify the walltime (e.g. 4 hours)
c.AdditionalProperties.WallTime = '02:00:00';

% Specify an account to use for MATLAB jobs (e.g. pXXXX, bXXXX, etc)
c.AdditionalProperties.AccountName = 'p32072';

% Specify a queue/partition to use for MATLAB jobs (e.g. short, normal, long)
c.AdditionalProperties.QueueName = 'short';

%% optional arguments but worth considering
% Specify memory to use for MATLAB jobs, per core (default: 4gb)
c.AdditionalProperties.MemUsage = '32gb';

% Specify number of nodes to use
c.AdditionalProperties.Nodes = 12;

% Specify e-mail address to receive notifications about your job
c.AdditionalProperties.EmailAddress = 'hsw@northwestern.edu';

% Require exclusive node
c.AdditionalProperties.RequireExclusiveNode = false;


job = c.batch('/home/hsw967/Programming/ca_imaging/SLURM_scripts/mutualinfo_openfield_shuff_4SLURM', 'CurrentFolder', '.');
