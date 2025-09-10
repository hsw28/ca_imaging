function plotRateMaskSummary(N, normalized, MICutoff)
% Summary plot using rate masks defined as N std above mean
% For An-2, An-1, and An across all rats
% If you want data normalized by mean rate, pass normalized==1
%
% Optional:
%   MICutoff – [] (off) or scalar in [0,1]. When provided, ONLY INCLUDE cells
%               whose MI_noCSUS15_shuff.MI_<date>(:,3) >= MICutoff.

if nargin < 1 || isempty(N), N = 1; end
if nargin < 2 || isempty(normalized), normalized = 0; end
if nargin < 3, MICutoff = []; end

ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
nRats = numel(ratNames);

% Preallocate stats
rate_task_mask = nan(nRats,1);
rate_nontask_mask = nan(nRats,1);
rate_task_mask_non = nan(nRats,1);
rate_nontask_mask_non = nan(nRats,1);
sem_task_mask = nan(nRats,1);
sem_nontask_mask = nan(nRats,1);
sem_task_mask_non = nan(nRats,1);
sem_nontask_mask_non = nan(nRats,1);
pvals_taskmask = nan(nRats,1);
pvals_nonmask = nan(nRats,1);

figure('Color','w','Position',[200 300 1100 450]);

for r = 1:nRats
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    idx = find(strcmp(dateList, rat.An), 1);
    if isempty(idx) || idx < 3, warning('%s: not enough days around An.', ratNames{r}); continue; end
    theseDays = dateList(idx-2:idx);

    rt1_all = []; rt2_all = []; rp1_all = []; rp2_all = [];

    for d = 1:3
        dateStr  = theseDays{d};
        Sfld     = ['CA_peaks_'  dateStr];
        Mfld     = ['ratemask_'  dateStr];
        Rfld     = ['rates_'     dateStr];
        MIroot   = ['MI_'        dateStr];
        MIparent = 'MI_noCSUS15_shuff';  % <— source requested

        if ~isfield(rat.Ca_peaks, Sfld), continue; end
        spikeMat = rat.Ca_peaks.(Sfld);

        % rate mask (if missing, keep all)
        if isfield(rat, 'ratemask') && isfield(rat.ratemask, Mfld)
            rateMask = rat.ratemask.(Mfld);
        else
            rateMask = true(size(spikeMat,1),1);
        end

        nNeurons = size(spikeMat,1);

        % normalization vector
        if normalized==1 && isfield(rat, 'rates') && isfield(rat.rates, Rfld)
            meanrate = rat.rates.(Rfld);
            % ensure row vector of length nNeurons
            if iscolumn(meanrate), meanrate = meanrate.'; end
            if numel(meanrate) < nNeurons
                meanrate(end+1:nNeurons) = NaN;
            end
        else
            meanrate = ones(1,nNeurons);
        end



        % ----- OPTIONAL MI INCLUDE FILTER -----
        if ~isempty(MICutoff) && isfield(rat, MIparent) && isfield(rat.(MIparent), MIroot)
            MIvals = rat.(MIparent).(MIroot);   % n×4; use (:,3)
            miCol3 = MIvals(:,3);
            keepMI = miCol3 >= MICutoff;
            % combine with rateMask
            rateMask = rateMask(:) & keepMI(:);
        end


        % loop neurons
        for ni = 1:nNeurons
            if ~rateMask(ni), continue; end

            spikes = spikeMat(ni,:);
            spikes = spikes(~isnan(spikes));
            if isempty(spikes), continue; end

            [rt1, rt2, rp1, rp2] = RateMaskVsTask_summary(rat, ni, dateStr, N);
            if any(isnan([rt1, rt2, rp1, rp2])), continue; end

            rt1_all(end+1) = rt1 ./ meanrate(ni); %#ok<*AGROW>
            rt2_all(end+1) = rt2 ./ meanrate(ni);
            rp1_all(end+1) = rp1 ./ meanrate(ni);
            rp2_all(end+1) = rp2 ./ meanrate(ni);
        end
    end

    % means & SEMs
    rate_task_mask(r)        = mean(rt1_all, 'omitnan');
    sem_task_mask(r)         = std(rt1_all, 0, 'omitnan') / sqrt(max(1,numel(rt1_all)));
    rate_nontask_mask(r)     = mean(rt2_all, 'omitnan');
    sem_nontask_mask(r)      = std(rt2_all, 0, 'omitnan') / sqrt(max(1,numel(rt2_all)));
    rate_task_mask_non(r)    = mean(rp1_all, 'omitnan');
    sem_task_mask_non(r)     = std(rp1_all, 0, 'omitnan') / sqrt(max(1,numel(rp1_all)));
    rate_nontask_mask_non(r) = mean(rp2_all, 'omitnan');
    sem_nontask_mask_non(r)  = std(rp2_all, 0, 'omitnan') / sqrt(max(1,numel(rp2_all)));

    % paired t-tests (when possible)
    if numel(rt1_all) >= 3 && numel(rt2_all) >= 3
        [~, p] = ttest(rt1_all, rt2_all);
        length(rt1_all)
        pvals_taskmask(r) = p;
    end
    if numel(rp1_all) >= 3 && numel(rp2_all) >= 3
        [~, p] = ttest(rp1_all, rp2_all);
        pvals_nonmask(r) = p;
    end
end

figure
% --- Plot Task-based mask ---
subplot(1,2,1); cla; hold on;
barData = [rate_task_mask, rate_nontask_mask];
bh = bar(barData, 'grouped');
x1 = bh(1).XEndPoints;
x2 = bh(2).XEndPoints;
errorbar(x1, barData(:,1), sem_task_mask,        'k', 'LineStyle','none');
errorbar(x2, barData(:,2), sem_nontask_mask,     'k', 'LineStyle','none');
xticks(1:nRats); xticklabels(ratNames);
ylabel('Firing Rate (Hz)');
title('Task vs Non-task in task-defined mask');
legend({'In Mask','Out of Mask'}, 'Location', 'northwest');
set(gca,'Box','off');
bh(1).FaceColor = [0.2 0.4 0.8];
bh(2).FaceColor = [0.7 0.7 0.7];
pvals_taskmask
for r = 1:nRats
    if ~isnan(pvals_taskmask(r))
        starsAt(x1(r), x2(r), max(barData(r,:))*1.05, pvals_taskmask(r));
    end
end

% --- Plot Non-task-based mask ---
subplot(1,2,2); cla; hold on;
barData2 = [rate_task_mask_non, rate_nontask_mask_non];
bh2 = bar(barData2, 'grouped');
x1 = bh2(1).XEndPoints;
x2 = bh2(2).XEndPoints;
errorbar(x1, barData2(:,1), sem_task_mask_non,       'k', 'LineStyle','none');
errorbar(x2, barData2(:,2), sem_nontask_mask_non,    'k', 'LineStyle','none');
xticks(1:nRats); xticklabels(ratNames);
ylabel('Firing Rate (Hz)');
title('Task vs Non-task in non-task-defined mask');
legend({'Task','Non-task'}, 'Location', 'northwest');
set(gca,'Box','off');
bh2(1).FaceColor = [0.2 0.4 0.8];
bh2(2).FaceColor = [0.7 0.7 0.7];
for r = 1:nRats
    if ~isnan(pvals_nonmask(r))
        starsAt(x1(r), x2(r), max(barData2(r,:))*1.05, pvals_nonmask(r));
    end
end




% -------- Optional per-cell scatter (last computed arrays) --------
% If you want them across all rats/days, hoist collectors outside the per-rat loop.
figure

% --- Task mask scatter ---
subplot(1,2,1); hold on; axis square
title('Per-cell: Task mask');
xlabel('In-trial rate (task mask)'); ylabel('Out-of-trial rate (task mask)');

x = rt1_all(:); y = rt2_all(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);

scatter(x, y, 25, 'filled');

if ~isempty(x) && ~isempty(y)
    % unity line
    m = max([x; y]);
    if ~isfinite(m) || m<=0, m = 1; end
    plot([0 m],[0 m],'k--','LineWidth',1);

    % best-fit line
    b  = polyfit(x, y, 1);
    xx = linspace(0, m, 100);
    plot(xx, polyval(b, xx), 'r-', 'LineWidth', 1.5);

    % stats
    [R, P] = corr(x, y, 'rows','complete');  % Pearson by default
    text(0.05, 0.90, sprintf('r = %.3f, p = %.3g', R, P), ...
        'Units','normalized','FontSize',12,'VerticalAlignment','top');
end

% --- Spatial mask scatter ---
subplot(1,2,2); hold on; axis square
title('Per-cell: Spatial mask');
xlabel('In-trial rate (spatial mask)'); ylabel('Out-of-trial rate (spatial mask)');

x = rp1_all(:); y = rp2_all(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);

scatter(x, y, 25, 'filled');

if ~isempty(x) && ~isempty(y)
    % unity line
    m = max([x; y]);
    if ~isfinite(m) || m<=0, m = 1; end
    plot([0 m],[0 m],'k--','LineWidth',1);

    % best-fit line
    b  = polyfit(x, y, 1);
    xx = linspace(0, m, 100);
    plot(xx, polyval(b, xx), 'r-', 'LineWidth', 1.5);

    % stats
    [R, P] = corr(x, y, 'rows','complete');
    text(0.05, 0.90, sprintf('r = %.3f, p = %.3g', R, P), ...
        'Units','normalized','FontSize',12,'VerticalAlignment','top');
end


end % function

% --------------- helpers ---------------
function starsAt(x1, x2, y, p)
xc = mean([x1,x2]);
if p < 1e-3
    txt = '***';
elseif p < 1e-2
    txt = '**';
elseif p < 0.05
    txt = '*';
else
    txt = '';
end
if ~isempty(txt)
    text(xc, y, txt, 'FontSize', 14, 'HorizontalAlignment','center');
end
end





function [rate1, rate2, rate3, rate4] = RateMaskVsTask_summary(animal, neuronIdx, dateStr, N)
% RateMaskVsTask_summary
%   Returns four proportions for one neuron on one day:
%     rate1: fraction of TASK spikes that land in TASK-defined mask
%     rate2: fraction of NON-TASK (running>=4) spikes that land in TASK-defined mask
%     rate3: fraction of TASK spikes that land in NON-TASK-defined mask
%     rate4: fraction of NON-TASK (running>=4) spikes that land in NON-TASK-defined mask
%
% Notes
%   - Task window is [CS, CS+2] s (no speed filter).
%   - Non-task map/Spikes use only running samples/spikes (speed >= 4 cm/s)
%     and exclude all CS→CS+2 s windows.

rate1 = NaN; rate2 = NaN; rate3 = NaN; rate4 = NaN;

% --- Pull data
pos      = animal.pos.(['pos_' dateStr]);         % [t x y]
spikeRow = animal.Ca_peaks.(['CA_peaks_' dateStr])(neuronIdx, :);
csTimes  = animal.CS_times.(['CS_' dateStr])(:);

% sanitize spikes
spikes = spikeRow(~isnan(spikeRow) & spikeRow>0);
if isempty(spikes), return, end

ts     = pos(:,1);
pos_xy = pos(:,2:3);

% --- Build trial mask on the position timebase
taskMask = false(size(ts));
for k = 1:numel(csTimes)
    t0 = csTimes(k);
    taskMask = taskMask | (ts >= t0 & ts < t0 + 2);
end
taskTimes = ts(taskMask);
taskPos   = pos(taskMask, :);     % used to build the task map


% --- Split spikes into task vs non-task (no speed filter here)
isInTask_spk = false(size(spikes));
for k = 1:numel(csTimes)
    t0 = csTimes(k);
    isInTask_spk = isInTask_spk | (spikes >= t0 & spikes < t0 + 2);
end
taskSpikes    = spikes(isInTask_spk);
nonTaskSpikes = spikes(~isInTask_spk);

% If either group is empty we can't compute any of the four rates
if isempty(taskSpikes) || isempty(nonTaskSpikes), return, end

% --- Velocity on its own timebase, then restrict NON-TRIAL pos/spikes by speed >= 4
vel        = ca_velocity(pos);           % [speed; time]
vel_time   = vel(2,:)';
vel_mag    = vel(1,:)';

% interpolate x,y to velocity timestamps to form a pos stream on vel_time
interp_x = interp1(ts, pos_xy(:,1), vel_time, 'linear', NaN);
interp_y = interp1(ts, pos_xy(:,2), vel_time, 'linear', NaN);

% mark velocity samples that fall inside any CS→CS+2 window
inTrial_vel = false(size(vel_time));
for k = 1:numel(csTimes)
    t0 = csTimes(k);
    inTrial_vel = inTrial_vel | (vel_time >= t0 & vel_time < t0 + 2);
end

% keep only running & non-trial velocity samples with valid position
velKeep   = (vel_mag >= 4) & ~inTrial_vel & isfinite(interp_x) & isfinite(interp_y);
goodpos   = [vel_time(velKeep), interp_x(velKeep), interp_y(velKeep)];
if size(goodpos,1) < 5, return, end

% also restrict NON-TASK spikes to ones that were run>=4 at spike time
spike_vel = interp1(vel_time, vel_mag, nonTaskSpikes, 'linear', 'extrap');
nonTaskSpikes = nonTaskSpikes(spike_vel >= 4);
if isempty(nonTaskSpikes), return, end

N=1;
% ---------- Build rate maps and masks ----------
% Task map (no speed filter; uses only taskPos)
rate_task = CA_normalizePosData(taskSpikes(:), taskPos, 2.5, 1);
rate_task(rate_task==0) = NaN;
mu_task = nanmean(rate_task(:));
sd_task = nanstd(rate_task(:));
mask_task = rate_task > (mu_task + N*sd_task);
mask_task = rate_task > 0;

% Non-task map (run>=4 & non-trial; uses goodpos)
rate_non = CA_normalizePosData(nonTaskSpikes(:), goodpos, 2.5, 1);
rate_non(rate_non==0) = NaN;
mu_non = nanmean(rate_non(:));
sd_non = nanstd(rate_non(:));
mask_non = rate_non > (mu_non + N*sd_non);

% If either mask is all-NaN or empty, bail
if ~any(mask_task(:)) || ~any(mask_non(:)), return, end

% Bin edges that correspond to each map (so discretize matches CA_normalizePosData’s grid)
Xedges_task = linspace(min(taskPos(:,2)), max(taskPos(:,2)), size(rate_task,2)+1);
Yedges_task = linspace(min(taskPos(:,3)), max(taskPos(:,3)), size(rate_task,1)+1);
Xedges_non  = linspace(min(goodpos(:,2)), max(goodpos(:,2)), size(rate_non,  2)+1);
Yedges_non  = linspace(min(goodpos(:,3)), max(goodpos(:,3)), size(rate_non,  1)+1);

% Helper to get spike (x,y) by interpolation onto the pos timebase
sx_task = interp1(ts, pos_xy(:,1), taskSpikes,    'linear', 'extrap');
sy_task = interp1(ts, pos_xy(:,2), taskSpikes,    'linear', 'extrap');
sx_non  = interp1(ts, pos_xy(:,1), nonTaskSpikes, 'linear', 'extrap');
sy_non  = interp1(ts, pos_xy(:,2), nonTaskSpikes, 'linear', 'extrap');

% ---------- Map spikes to TASK mask ----------
bx = discretize(sx_task, Xedges_task);  by = discretize(sy_task, Yedges_task);
validT = bx>=1 & by>=1 & bx<=size(rate_task,2) & by<=size(rate_task,1);
if any(validT)
    in_task_T = mask_task(sub2ind(size(rate_task), by(validT), bx(validT)));
    rate1 = sum(in_task_T) / sum(validT);              % Task spikes in task-defined mask
end

bx = discretize(sx_non,  Xedges_task);  by = discretize(sy_non,  Yedges_task);
validN = bx>=1 & by>=1 & bx<=size(rate_task,2) & by<=size(rate_task,1);
if any(validN)
    in_task_N = mask_task(sub2ind(size(rate_task), by(validN), bx(validN)));
    rate2 = sum(in_task_N) / sum(validN);              % Non-task (run) spikes in task-defined mask
end

% ---------- Map spikes to NON-TASK mask ----------
bx = discretize(sx_task, Xedges_non);   by = discretize(sy_task, Yedges_non);
validT = bx>=1 & by>=1 & bx<=size(rate_non,2) & by<=size(rate_non,1);
if any(validT)
    in_non_T = mask_non(sub2ind(size(rate_non), by(validT), bx(validT)));
    rate3 = sum(in_non_T) / sum(validT);               % Task spikes in non-task-defined mask
end

bx = discretize(sx_non,  Xedges_non);   by = discretize(sy_non,  Yedges_non);
validN = bx>=1 & by>=1 & bx<=size(rate_non,2) & by<=size(rate_non,1);
if any(validN)
    in_non_N = mask_non(sub2ind(size(rate_non), by(validN), bx(validN)));
    rate4 = sum(in_non_N) / sum(validN);               % Non-task (run) spikes in non-task-defined mask
end
end
