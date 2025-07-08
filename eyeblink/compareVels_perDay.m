function f = compareVels_perDay(perDayStruct)
% compares velocities in trace and non trace period
%% use this after running perDayStruct = alignEMGandVelocity_perDay(ratName)

    if isempty(perDayStruct)
        warning('No data to plot.');
        return;
    end

    % Organize data
    trace_vals_all = {};
    nontrace_vals_all = {};
    labels = {};
    rats = {};
    ps = [];

    for i = 1:length(perDayStruct)
        trace_vals_all{end+1} = perDayStruct(i).trace_vel_vals;
        nontrace_vals_all{end+1} = perDayStruct(i).nontrace_vel_vals;
        labels{end+1} = sprintf('%s\n%s', perDayStruct(i).rat, perDayStruct(i).date);
        ps(end+1) = perDayStruct(i).p;
    end



    % Flatten for boxplot
    all_vals = [];
    group = [];
    tick_labels = {};
    xtick_positions = [];

    for i = 1:length(trace_vals_all)
        ti = (i-1)*2 + 1;
        all_vals = [all_vals; nontrace_vals_all{i}(:); trace_vals_all{i}(:)];

        group = [group; repmat(ti, numel(nontrace_vals_all{i}(:)), 1); ...
                       repmat(ti+1, numel(trace_vals_all{i}(:)), 1)];

        tick_labels{ti} = '';  % blank space for alignment
        tick_labels{ti+1} = labels{i};
        xtick_positions(end+1) = ti + 0.5;
    end

    % Colors
    color_non = [0.6 0.6 0.6];
    color_trace = [0.2 0.5 0.9];

    % Plot
    figure('Color','w', 'Position', [300 300 100+90*length(ps), 400]);
    boxplot(all_vals, group, 'Colors', 'k', 'Symbol', '.k', 'Widths', 0.6);
    hold on;

    % Color fill boxes
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        if mod(length(h) - j + 1, 2) == 0
            patch(get(h(j), 'XData'), get(h(j), 'YData'), color_trace, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        else
            patch(get(h(j), 'XData'), get(h(j), 'YData'), color_non, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
    end

    set(gca, 'XTick', xtick_positions, 'XTickLabel', tick_labels, ...
        'XTickLabelRotation', 30, 'FontSize', 10);
    ylabel('Velocity (cm/s)');
    title('Trace vs. Non-Trace Velocity Per Day');
    box off;

    % Add p-value annotations
    ylim_buffer = 1.1 * max(all_vals);
    for i = 1:length(ps)
        x = xtick_positions(i);
        y = ylim_buffer;
        pval = ps(i);
        if pval < 0.001
            sig = '***';
        elseif pval < 0.01
            sig = '**';
        elseif pval < 0.05
            sig = '*';
        else
            sig = 'n.s.';
        end
        text(x, y, sig, 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
end
