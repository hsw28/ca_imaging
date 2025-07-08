function plotPlaceAndTrialTracesTogether(animal, neuronVec, dateStr)
DEPPPPPPP
% For each neuron, generate a single figure with place maps and trial-aligned traces

for ni = 1:numel(neuronVec)
    neuronIdx = neuronVec(ni);

    figure('Color','w', 'Position', [300 300 1400 800]);

    % Create 2 axes for the place map
    ax1 = subplot(2, 3, 1);
    ax2 = subplot(2, 3, 2);

    % Plot place maps into those axes
    plot_place_vs_CSUS(animal, neuronIdx, dateStr, [ax1, ax2]);

    % Plot trial traces (raster, PSTH, calcium)
    subplot(2, 3, 4); hold on;
    subplot(2, 3, 5); hold on;
    subplot(2, 3, 6); hold on;

    % Call trial-aligned plotter with `part_of_diff_fig` flag to avoid new figure
    plotTraceRasterAndCalcium(...
        animal.Ca_peaks.(['CA_peaks_' dateStr]), ...
        animal.CS_times.(['CS_' dateStr]), ...
        [-1, 2], ...
        animal.Ca_traces.(['CA_traces_' dateStr]), ...
        animal.Ca_ts.(['CA_time_' dateStr]), ...
        neuronIdx, true ...
    );
end
end
