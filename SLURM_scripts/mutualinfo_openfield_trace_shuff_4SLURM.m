function mutualinfo_openfield_trace_shuff_4SLURM

addpath(pwd);
addpath(genpath('/home/hsw967/Programming/ca_imaging'));
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/matlab/'));
addpath(genpath('/home/hsw967/Programming/data_analysis/hannah-in-use/include'));



%file allvariables.mat should contain
  %all_traces
  %MI
  %MI_trace
  %peaks
  %pos
  %Ca_ts

fprintf('loading traces')
allvariables = load('allvariables.mat');
calcium_traces = allvariables.Ca_traces;
clearvars allvariables

fprintf('loading pos')
pos = load('pos.mat');
pos_structure = pos.pos;

fprintf('loading MI')
MI = load('MI_CSUS.mat');
ca_MI = MI.MI_trace;

fprintf('loading timestamps')
Ca_ts = load('Ca_ts.mat')
ca_ts = Ca_ts.Ca_ts;


% Calculate the size of the structure
structure_info = whos('calcium_traces');
structure_size_gb = structure_info.bytes / (1024^3);  % Convert to GB
% Check if the structure size exceeds 4 GB

%{
if structure_size_gb > 4
    % Determine the number of fields in the structure
    num_fields = numel(fieldnames(calcium_traces));
    % Calcilate the number of fields for each new structure (approximately half)
    num_fields_per_structure = ceil(num_fields / 2);

    % Split the structure into two new structures
    fprinf('splitting structure')
    field_namesCT = fieldnames(calcium_traces);
    field_namesPS = fieldnames(pos_structure);
    field_namesMI = fieldnames(ca_MI);
    field_namesTS = fieldnames(ca_ts);

    % Create the first new structure
    calcium_traces1 = struct();
    for i = 1:num_fields_per_structure
        field_name = field_namesCT{i};
        calcium_traces1.(field_name) = calcium_traces.(field_name);
        field_name = field_namesPS{i};
        pos_structure1.(field_name) = pos_structure.(field_name);
        field_name = field_namesMI{i};
        ca_MI1.(field_name) = ca_MI.(field_name);
        field_name = field_namesTS{i};
        ca_ts1.(field_name) = ca_ts.(field_name);
    end

    % Create the second new structure


    calcium_traces2 = struct();
    for i = (num_fields_per_structure + 1):num_fields
      field_name = field_namesCT{i};
      calcium_traces2.(field_name) = calcium_traces.(field_name);
      field_name = field_namesPS{i};
      pos_structure2.(field_name) = pos_structure.(field_name);
      field_name = field_namesMI{i};
      ca_MI2.(field_name) = ca_MI.(field_name);
      field_name = field_namesTS{i};
      ca_ts2.(field_name) = ca_ts.(field_name);
    end


    % Clear the original structure to free up memory
    clear calcium_traces;
    clear pos_structure
    clear ca_MI
    clear ca_ts

    f1 = mutualinfo_openfield_trace_shuff(calcium_traces1, pos_structure1, 2, 2.5, 500, ca_MI1, ca_ts1);
    % Save the output to a .mat file
    MI_trace_shuff = f1;



    f2 = mutualinfo_openfield_trace_shuff(calcium_traces1, pos_structure1, 2, 2.5, 500, ca_MI1, ca_ts1);
    % Save the output to a .mat file
    MI_trace_shuff = f2;


  else
%}
    f = mutualinfo_openfield_trace_shuff(calcium_traces, pos_structure, 2, 2.5, ca_ts, 500, ca_MI);
    % Save the output to a .mat file
    MI_trace_shuff = f;


%end





end
