function fixed_US = alignOnsetToEMGRise(emg_data, emg_ts, predicted_onsets, show)
  % ALIGNONSETSTRICTMINTOMAX_SAFE
  % Detect true EMG onset using last local min to peak monotonic rise constraint.
  % Now includes timing safeguards and rise duration limit.



  emg_data = filter_emg(emg_data);
  
  if nargin < 4
      show = false;
  end

  % Parameters
  fs = 1 / nanmean(diff(emg_ts));
  search_pre  = 0.15;   % 100 ms backward window

  true_onsets = nan(size(predicted_onsets));
  delta_onsets = nan(size(predicted_onsets));


  for i = 1:numel(predicted_onsets)

      pred_time = predicted_onsets(i);
      idx_pred = find(emg_ts >= pred_time, 1)

      searcharea = round(fs*search_pre);

      slopes = diff(emg_data((idx_pred-searcharea):(idx_pred)))

      neg = find(slopes<0);
      want = max(neg)+idx_pred-searcharea;

      true_onsets(i) = emg_ts(want);
      delta_onsets(i) = pred_time-emg_ts(want);
    end

    av = nanmean(delta_onsets)
    fixed_US = predicted_onsets+av;
