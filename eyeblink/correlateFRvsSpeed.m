function corrVals = correlateFRvsSpeed(animal, neuronVec, dateStr, win)
% Correlate firing rate vs speed across trials during conditioning period
% win = [0 2] by default (post-CS window)
%ex:   correlateFRvsSpeed('rat0314', 'allcells', '2023_05_22', [0,2])

if ischar(animal)
    animal = evalin('base', animal);
end

if nargin < 4, win = [0 2]; end

  spikeMat = animal.Ca_peaks.(['CA_peaks_' dateStr]);
  if isequal(neuronVec, 'allcells')
      neuronVec = 1:size(spikeMat, 1);
  end
  corrVals = nan(numel(neuronVec), 1);


% Extract data
ts = animal.pos.(['pos_' dateStr])(:,1);
xy = animal.pos.(['pos_' dateStr])(:,2:3);
spikeMat = animal.Ca_peaks.(['CA_peaks_' dateStr]);
cs_times = animal.CS_times.(['CS_' dateStr]);

% Compute running speed
dt = diff(ts);
dx = diff(xy);
speed = [0; sqrt(sum(dx.^2,2)) ./ dt];  % cm/s
speed(isinf(speed) | isnan(speed)) = 0;

% Interpolate speed for arbitrary timestamps
speed_interp = @(t) interp1(ts, speed, t, 'linear', 'extrap');

for ni = 1:numel(neuronVec)
    neuron = neuronVec(ni);
    spikes = spikeMat(neuron,:); spikes = spikes(~isnan(spikes));
    if numel(spikes) < 3, continue; end

    fr_per_trial = nan(size(cs_times));
    speed_per_trial = nan(size(cs_times));

    for t = 1:numel(cs_times)
        t_start = cs_times(t) + win(1);
        t_end   = cs_times(t) + win(2);

        % Firing rate in window
        nSpikes = sum(spikes >= t_start & spikes <= t_end);
        fr_per_trial(t) = nSpikes / (t_end - t_start);

        % Mean speed in window
        mask = ts >= t_start & ts <= t_end;
        if any(mask)
            speed_per_trial(t) = mean(speed(mask));
        end
    end

    % Correlate if enough valid trials
    valid = ~isnan(fr_per_trial) & ~isnan(speed_per_trial);
    if sum(valid) >= 5
        corrVals(ni) = corr(fr_per_trial(valid)', speed_per_trial(valid)', 'type', 'Pearson');
    end
end
end
