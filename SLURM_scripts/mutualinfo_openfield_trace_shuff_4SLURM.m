function mutualinfo_openfield_trace_shuff_4SLURM

% Start a parallel pool with the specified number of workers
numWorkers = 8; % This should match the number of CPUs requested in the SLURM script
poolobj = gcp('nocreate'); % Check if the pool already exists
if isempty(poolobj)
    poolobj = parpool(numWorkers);
end

pos_structure = load('pos.mat');
calcium_traces = load('traces.mat')
ca_MI = load('MI_trace.mat')

f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 2, 2.5, 500, ca_MI)
% Save the output to a .mat file
save('mutualinfo_shuff_output.mat', 'MI_shuff');
