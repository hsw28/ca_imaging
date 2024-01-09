function trace_shuff_4SLURM
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
