
function mutualinfo_openfield_trace_shuff_4SLURM
  outputFile = fopen('mutualinfo_shuff_output.log', 'w');
  try
  maxNumCompThreads(str2num(getenv('SLURM_NPROCS')));

  % Refresh MATLAB's toolbox cache
  rehash toolboxcache;

  % Custom configuration for the cluster (if required)
  configCluster;

%{
  c = parcluster;
  c.AdditionalProperties.AccountName = 'p32072'; % Replace with your actual account name
  c.AdditionalProperties.WallTime = '24:00:00'; % Set to match the wall time in your SLURM script
  c.AdditionalProperties.QueueName = 'normal';
  c.AdditionalProperties.MemUsage = '64gb';


  % Set the number of computational threads to the number of allocated CPUs
  if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
      maxNumCompThreads(str2num(getenv('SLURM_CPUS_PER_TASK')));
  end
%}

addpath(pwd);
addpath(genpath('/home/hsw967/Programming/ca_imaging'));
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));


%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos
  %ca_ts

allvariables = load('allvariables.mat');
pos_structure = allvariables.pos;
calcium_traces = allvariables.all_traces;

MI_trace = load('MI_trace.mat');
ca_MI = MI_trace.MI_trace;

Ca_ts = load('Ca_ts.mat')
ca_ts = Ca_ts.Ca_ts;

f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 2, 2.5, 500, ca_MI, ca_ts);
% Save the output to a .mat file
MI_trace_shuff = f;

try
    save('./mutualinfo_trace_shuff_output.mat', 'MI_trace_shuff');
catch e
    disp('Error occurred while saving the file:');
    disp(e.message);
end

fprintf(outputFile, 'Computation completed successfully.\n');
catch e
fprintf(outputFile, 'An error occurred: %s\n', e.message);
end
 fclose(outputFile);
end
