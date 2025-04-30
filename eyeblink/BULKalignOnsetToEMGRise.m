function [corrected_US corrected_CS] = BULKalignOnsetToEMGRise(EMG_struct, EMG_ts_struct, us_times_struct)

% CORRECTPREDICTEDONSETSSTRUCT
% Adjusts predicted onset times based on EMG by aligning to last negative slope.
%
% Inputs:
%   EMG_struct             - Struct of EMG signals (filtered recommended)
%   EMG_ts_struct          - Struct of timestamps (in seconds)
%   predicted_onsets_struct - Struct of predicted onset times (e.g., US times)
%
% Output:
%   corrected_onsets_struct - Struct of corrected onset times

predicted_onsets_struct = us_times_struct;

fields_EMG = fieldnames(EMG_struct);
fields_TS = fieldnames(EMG_ts_struct);
fields_US = fieldnames(predicted_onsets_struct);


if length(fields_EMG) ~= length(fields_TS) | length(fields_EMG) ~= length(fields_US)
  error('your structure lengths do night align')
end


corrected_onsets_struct = struct();

for i = 1:numel(fields_EMG)
    field_EMG = fields_EMG{i};
    field_TS = fields_TS{i};
    field_US = fields_US{i};

    emg_data = EMG_struct.(field_EMG);
    emg_ts = EMG_ts_struct.(field_TS);
    predicted_onsets = predicted_onsets_struct.(field_US);


    if length(emg_data)<4 || length(predicted_onsets)<4
        corrected_US.(field_US) = [];

        newname = strrep(field_US, 'US_', 'CS_');

        corrected_CS.(newname) = [];
        continue;
    end

    emg_data = filter_emg(emg_data);

    fs = 1 / nanmean(diff(emg_ts));
    search_pre = 0.15;  % seconds to search back
    delta_onsets = nan(size(predicted_onsets));
    true_onsets = nan(size(predicted_onsets));

    for t = 1:numel(predicted_onsets)
        pred_time = predicted_onsets(t);
        idx_pred = find(emg_ts >= pred_time, 1);
        if isempty(idx_pred) || idx_pred <= 1
            continue;
        end

        searcharea = round(fs * search_pre);
        if (idx_pred - searcharea) < 1 || idx_pred > numel(emg_data)
            continue;
        end

        segment = emg_data((idx_pred - searcharea):idx_pred);
        slopes = diff(segment);
        neg = find(slopes < 0);

        if isempty(neg)
            continue;
        end

        want = max(neg) + idx_pred - searcharea;
        true_onsets(t) = emg_ts(want);
        delta_onsets(t) = pred_time - emg_ts(want);
    end

    avg_shift = nanmean(delta_onsets);
    corrected_US.(field_US) = predicted_onsets + avg_shift;

    newname = strrep(field_US, 'US_', 'CS_');

    corrected_CS.(newname) = predicted_onsets + avg_shift -.5 - .010;
end
end
