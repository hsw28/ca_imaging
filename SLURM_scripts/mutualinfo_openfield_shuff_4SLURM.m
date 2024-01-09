function mutualinfo_openfield_shuff_4SLURM

  % Refresh MATLAB's toolbox cache
  rehash toolboxcache;

  % Custom configuration for the cluster (if required)
  configCluster;

  c = parcluster;
  c.AdditionalProperties.AccountName = 'p32072'; % Replace with your actual account name
  c.AdditionalProperties.WallTime = '24:00:00'; % Set to match the wall time in your SLURM script
  c.AdditionalProperties.QueueName = 'normal';
  c.AdditionalProperties.MemUsage = '64gb';

  addpath(pwd);
  addpath(genpath('/home/hsw967/Programming/ca_imaging'));
  addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));


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
spikes = load('spikes.mat')
ca_MI = load('MI.mat')

f = mutualinfo_openfield_shuff(spikes, pos_structure, 2, 2.5, 500, ca_MI)
MI_shuff = f;
% Save the output to a .mat file
save('mutualinfo_shuff_output.mat', 'MI_shuff');
end
