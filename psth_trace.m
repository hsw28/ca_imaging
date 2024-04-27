function f = psth_trace(triggerTimes, trace, Ca_time, sp)

% Assuming `triggerTimes` is a list of trigger times in seconds
% `LFP` is your list/array of LFP values
% `Fs` is the sampling rate of your LFP signal
%ex: psth_trace(rat0314.CS_times.CS_2023_05_25, rat0314.Ca_traces.CA_traces_2023_05_25(1,:), rat0314.Ca_ts.CA_time_2023_05_25)


Fs = 7.5; % Replace with your actual sampling rate
preTriggerTime = 1; % 1 second before the trigger
postTriggerTime = 2; % 1 second after the trigger


% Number of samples before and after trigger to include
samplesBefore = 8;
samplesAfter = 12;

% Initialize matrix to hold the LFP segments
numTriggers = length(triggerTimes);
LFPsegments = NaN(numTriggers, samplesBefore + samplesAfter + 1);
Ca_time(:,2) = Ca_time(:,2)./1000;


% Extract LFP segments around the trigger points
for i = 1:numTriggers

    curtrig = triggerTimes(i);
    [c index] = (min(abs(curtrig-Ca_time(:,2))));
    index = ceil(index./2);

    % Get the range of indices for the current trigger
    startIndex = (index) - samplesBefore;
    endIndex = (index) + samplesAfter;

    % Check for bounds, in case the trigger is at the edges
    if startIndex < 1
        continue
    end
    if endIndex > length(trace)
        continue
    end

    % Extract LFP segment and store it
    LFPsegments(i, :) = trace(startIndex:endIndex);
end

% Compute the average LFP signal across all triggers
averageLFP = nanmean(LFPsegments, 1);
f = averageLFP;

% Time vector for plotting
timeVector = (-samplesBefore:samplesAfter) / Fs;

% Plot the average LFP
subplot(4,2,sp)
plot(timeVector, averageLFP);
vline(0)
vline(.25)
vline(.75)
vline(.85)
