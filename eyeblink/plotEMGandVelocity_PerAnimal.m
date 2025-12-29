function plotEMGandVelocity_PerAnimal(t_interp, mean_emg, sem_or_std_emg, mean_vel, sem_or_std_vel, ratName, spreadType, varargin)
% plotEMGandVelocity_PerAnimal
%   Plot EMG (right y-axis) and speed (left y-axis) aligned to CS.
%
% Supports two modes:
%   (A) Pooled mode (original-style):
%       plotEMGandVelocity_PerAnimal(t, mE, sE, mV, sV, ratName, spreadType)
%
%   (B) Per-day mode:
%       plotEMGandVelocity_PerAnimal(t, perDayOut, [], [], [], ratName, spreadType)
%       where perDayOut is struct array from alignEMGandVelocity(..., perDay=true)
%
% Name-Value options:
%   'TraceShade'   [0 0.75]    (seconds relative to CS)
%   'Axes'         []          (pooled: 1 axis handle; per-day: 1×nDays axes handles)
%   'MakeFigure'   true        (if no Axes provided, create figure)
%   'YLimVel'      []          (e.g., [-5 90])
%   'YLimEMG'      []          (e.g., [-30000 30000])

    if nargin < 7 || isempty(spreadType)
        spreadType = 'sem';
    end

    % ---- parse options ----
    p = inputParser;
    p.addParameter('TraceShade', [0 0.75], @(x)isnumeric(x)&&numel(x)==2);
    p.addParameter('Axes', [], @(x) isempty(x) || all(ishghandle(x)));
    p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
    p.addParameter('YLimVel', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    p.addParameter('YLimEMG', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
    p.parse(varargin{:});
    traceShade = p.Results.TraceShade;
    axIn       = p.Results.Axes;
    makeFig    = p.Results.MakeFigure;
    yLimVel    = p.Results.YLimVel;
    yLimEMG    = p.Results.YLimEMG;

    % ---- label suffix ----
    switch lower(spreadType)
        case 'sem'
            labelSuffix = '± SEM';
        case 'std'
            labelSuffix = '± STD';
        otherwise
            error('spreadType must be "sem" or "std"');
    end

    % ---- colors ----
    emgColor   = [0.6 0.1 0.8];
    velColor   = [0.1 0.6 0.9];
    traceColor = [0.9 0.9 0.9];

    % ---- detect per-day mode ----
    isPerDay = isstruct(mean_emg) && ~isempty(mean_emg) && ...
               isfield(mean_emg,'mean_emg') && isfield(mean_emg,'mean_vel');

    if isPerDay
        perDayOut = mean_emg;
        nDays = numel(perDayOut);

        if isempty(axIn)
            if makeFig
                figure('Color','w');
            end
            ax = gobjects(1,nDays);
            for d = 1:nDays
                ax(d) = subplot(1,nDays,d);
            end
        else
            ax = axIn(:).';
            if numel(ax) ~= nDays
                error('Axes must have length nDays (%d).', nDays);
            end
        end

        for d = 1:nDays
            axes(ax(d)); %#ok<LAXES>
            cla(ax(d));
            hold(ax(d),'on');

            if isempty(perDayOut(d).mean_emg) || isempty(perDayOut(d).mean_vel)
                title(ax(d), sprintf('%s: %s (no valid trials)', ratName, perDayOut(d).day), 'Interpreter','none');
                xlim(ax(d), [t_interp(1) t_interp(end)]);
                yyaxis(ax(d),'left');  line(ax(d), [0 0], ylim(ax(d)), 'Color','k','LineStyle','--');
                set(ax(d), 'Box','off', 'FontSize', 12);
                continue;
            end

            mE = perDayOut(d).mean_emg;
            sE = perDayOut(d).spread_emg;
            mV = perDayOut(d).mean_vel;
            sV = perDayOut(d).spread_vel;

            % -------- velocity (left) --------
            yyaxis(ax(d),'left');
            fill(ax(d), [t_interp fliplr(t_interp)], ...
                      [mV + sV, fliplr(mV - sV)], ...
                      velColor, 'FaceAlpha',0.3, 'EdgeColor','none');
            plot(ax(d), t_interp, mV, 'Color', velColor, 'LineWidth', 1);

            if ~isempty(yLimVel)
                ylim(ax(d), yLimVel);
            else
                ymin = nanmin(mV - sV);
                ymax = nanmax(mV + sV);
                if ~(isnan(ymin) || isnan(ymax) || ymin==ymax)
                    ylim(ax(d), [ymin*0.9, ymax*1.1]);
                end
            end

            addShade(ax(d), traceShade, traceColor, 0.2);
            ylabel(ax(d), sprintf('Speed (cm/s) %s', labelSuffix));

            % -------- EMG (right) --------
            yyaxis(ax(d),'right');
            fill(ax(d), [t_interp fliplr(t_interp)], ...
                      [mE + sE, fliplr(mE - sE)], ...
                      emgColor, 'FaceAlpha',0.3, 'EdgeColor','none');
            plot(ax(d), t_interp, mE, 'Color', emgColor, 'LineWidth', 1);

            if ~isempty(yLimEMG)
                ylim(ax(d), yLimEMG);
            else
                ymin = nanmin(mE - sE);
                ymax = nanmax(mE + sE);
                if ~(isnan(ymin) || isnan(ymax) || ymin==ymax)
                    ylim(ax(d), [ymin*0.9, ymax*1.1]);
                end
            end

            addShade(ax(d), traceShade, traceColor, 0.2);
            ylabel(ax(d), sprintf('EMG (a.u.) %s', labelSuffix));

            % -------- cosmetics --------
            xlim(ax(d), [t_interp(1) t_interp(end)]);
            yyaxis(ax(d),'left');
            line(ax(d), [0 0], ylim(ax(d)), 'Color','k','LineStyle','--');

            title(ax(d), sprintf('%s: %s', ratName, perDayOut(d).day), 'Interpreter','none');
            set(ax(d), 'Box','off', 'FontSize', 12);

            if d == 1
                % only left-most gets x-label if you prefer; comment out if you want all
                xlabel(ax(d), 'Time from CS onset (s)');
            else
                xlabel(ax(d), 'Time from CS onset (s)');
            end
        end

        linkaxes(ax,'x');
        return;
    end

    % ======================
    % pooled mode
    % ======================
    if isempty(axIn)
        if makeFig
            figure('Color','w');
        end
        ax = gca;
    else
        ax = axIn(1);
    end

    axes(ax); %#ok<LAXES>
    cla(ax);
    hold(ax,'on');

    % velocity left
    yyaxis(ax,'left');
    fill(ax, [t_interp fliplr(t_interp)], ...
              [mean_vel + sem_or_std_vel, fliplr(mean_vel - sem_or_std_vel)], ...
              velColor, 'FaceAlpha',0.3, 'EdgeColor','none');
    plot(ax, t_interp, mean_vel, 'Color', velColor, 'LineWidth', 1);

    if ~isempty(yLimVel)
        ylim(ax, yLimVel);
    else
        ymin = nanmin(mean_vel - sem_or_std_vel);
        ymax = nanmax(mean_vel + sem_or_std_vel);
        if ~(isnan(ymin) || isnan(ymax) || ymin==ymax)
            ylim(ax, [ymin*0.9, ymax*1.1]);
        end
    end
    addShade(ax, traceShade, traceColor, 0.2);
    ylabel(ax, sprintf('Speed (cm/s) %s', labelSuffix));

    % EMG right
    yyaxis(ax,'right');
    fill(ax, [t_interp fliplr(t_interp)], ...
              [mean_emg + sem_or_std_emg, fliplr(mean_emg - sem_or_std_emg)], ...
              emgColor, 'FaceAlpha',0.3, 'EdgeColor','none');
    plot(ax, t_interp, mean_emg, 'Color', emgColor, 'LineWidth', 1);

    if ~isempty(yLimEMG)
        ylim(ax, yLimEMG);
    else
        ymin = nanmin(mean_emg - sem_or_std_emg);
        ymax = nanmax(mean_emg + sem_or_std_emg);
        if ~(isnan(ymin) || isnan(ymax) || ymin==ymax)
            ylim(ax, [ymin*0.9, ymax*1.1]);
        end
    end
    addShade(ax, traceShade, traceColor, 0.2);
    ylabel(ax, sprintf('EMG (a.u.) %s', labelSuffix));

    % CS line + cosmetics
    xlim(ax, [t_interp(1) t_interp(end)]);
    yyaxis(ax,'left');
    line(ax, [0 0], ylim(ax), 'Color','k', 'LineStyle','--');
    xlabel(ax, 'Time from CS onset (s)');
    title(ax, sprintf('%s: EMG and Speed Aligned to CS (%s)', ratName, upper(spreadType)), 'Interpreter','none');
    set(ax, 'Box','off', 'FontSize', 12);
end

% ---- local helper: add shaded patch using current ylim ----
function addShade(ax, shadeWin, shadeColor, alphaVal)
    yL = ylim(ax);
    fill(ax, [shadeWin(1) shadeWin(2) shadeWin(2) shadeWin(1)], ...
              [yL(1) yL(1) yL(2) yL(2)], ...
              shadeColor, 'EdgeColor','none', 'FaceAlpha', alphaVal);

end
