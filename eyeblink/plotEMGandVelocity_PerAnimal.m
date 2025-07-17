function plotEMGandVelocity_PerAnimal(t_interp, mean_emg, sem_or_std_emg, mean_vel, sem_or_std_vel, ratName, spreadType)
    %takes input from alignEMGandVelocity 
    %spreadType: 'sem' or 'std' (default = 'sem')

    if nargin < 7
        spreadType = 'sem';
    end


    % Label for y-axis and figure title
    switch lower(spreadType)
        case 'sem'
            labelSuffix = '± SEM';
        case 'std'
            labelSuffix = '± STD';
        otherwise
            error('spreadType must be "sem" or "std"');
    end

    % Colors
    emgColor = [0.6 0.1 0.8];
    velColor = [0.1 0.6 0.9];
    traceColor = [0.9 0.9 0.9];

%    figure('Color', 'w', 'Position', [300 300 600 400]);


    % === Panel 1: EMG ===
    hold on;
    yyaxis right
    size([t_interp fliplr(t_interp)])
    size([mean_emg + sem_or_std_emg, fliplr(mean_emg - sem_or_std_emg)])
    fill([t_interp fliplr(t_interp)], [mean_emg + sem_or_std_emg, fliplr(mean_emg - sem_or_std_emg)], emgColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(t_interp, mean_emg, 'Color', emgColor, 'LineWidth', 1);
    yLimits = ylim;
    fill([0 0.75 0.75 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
         traceColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    ylabel(sprintf('EMG (a.u.) %s', labelSuffix));
    title(sprintf('%s: EMG and Speed Aligned to CS (%s)', ratName, upper(spreadType)));
    xlim([t_interp(1), t_interp(end)]);
    ymin = nanmin(mean_emg - sem_or_std_emg);
    ymax = nanmax(mean_emg + sem_or_std_emg);
    if isnan(ymin) || isnan(ymax) || ymin == ymax
        ylim([-1 1]);
    else
        ylim([ymin*0.9, ymax*1.1]);
    end
    ylim([-30000,30000])
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
    set(gca, 'Box', 'off', 'FontSize', 12);

    % === Panel 2: Velocity ===
    yyaxis left
    ylim([-5 90])
    fill([t_interp fliplr(t_interp)], ...
         [mean_vel + sem_or_std_vel, fliplr(mean_vel - sem_or_std_vel)], ...
         velColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(t_interp, mean_vel, 'Color', velColor, 'LineWidth', 1);
    yLimits = ylim;
    fill([0 0.75 0.75 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
         traceColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    ylabel(sprintf('Speed (cm/s) %s', labelSuffix));
    xlabel('Time from CS onset (s)');
    xlim([t_interp(1), t_interp(end)]);
    ymin = nanmin(mean_vel - sem_or_std_vel);
    ymax = nanmax(mean_vel + sem_or_std_vel);
    if isnan(ymin) || isnan(ymax) || ymin == ymax
        ylim([-1 1]);
    else
        ylim([ymin*0.9, ymax*1.1]);
    end
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
    set(gca, 'Box', 'off', 'FontSize', 12);
        ylim([-5 90])

    linkaxes(findall(gcf,'Type','axes'), 'x');
end
