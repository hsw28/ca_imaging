function f = run_speedBinMatched_epoch()

ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};

%R_CS    = run_speedBinMatched(ratNames, 'win',[0     0.25], 'binSize',.25/2, 'Plot', false, 'SpeedBinCount',1,'SpeedBinWidth',[]);   % CS
%R_trace = run_speedBinMatched(ratNames, 'win',[0.25  0.75], 'binSize',.1, 'Plot', false, 'SpeedBinCount',2,'SpeedBinWidth',[]);   % Trace
%R_US    = run_speedBinMatched(ratNames, 'win',[0.75  0.85], 'binSize',.1, 'Plot', false);   % US
R_CS_trace_US    = run_speedBinMatched(ratNames, 'win',[0 0.85], 'binSize',.1, 'Plot', false, 'ExcludeWin',[0 2]);   % CS
R_post  = run_speedBinMatched(ratNames, 'win',[0.85  2.00],'binSize',.1, 'Plot', false, 'ExcludeWin',[0 2]);   % Post
R_post2  = run_speedBinMatched(ratNames, 'win',[2  4],'binSize',.1, 'Plot', false, 'ExcludeWin',[0 4]);   % Post
R_post4  = run_speedBinMatched(ratNames, 'win',[4  6],'binSize',.1, 'Plot', false, 'ExcludeWin',[0 6);   % Post

epochRs    = {R_CS_trace_US, R_post,R_post2, R_post4};
epochNames = {'Tone+Trace+Shock','Post .85-2s','Post 2-4s','Post 4-6s'};

%epochRs    = {R_CS, R_trace, R_US, R_post,R_post2, R_post4};
%epochNames = {'CS','Trace','US','Post .85-2s','Post 2-4s','Post 4-6s'};

plot_epoch_foldchange_cellLevel(epochRs, epochNames, 'metric','fold');
plot_epoch_foldchange_POP(epochRs,       epochNames, 'metric','fold');

end

function plot_epoch_foldchange_cellLevel(epochRs, epochNames, varargin)
% plot_epoch_foldchange_cellLevel
%   Epoch-wise *cell-level* fold change bars.
%   Inputs:
%     epochRs    : cell array of R structs, one per epoch (output of run_speedBinMatched)
%     epochNames : cellstr of labels, e.g. {'CS','Trace','US','Post'}
%
%   For each epoch:
%     - Compute per-rat fold change (Trial / Non-trial) using tested cells.
%     - Average fold change across rats, plot mean ± SEM across rats.

p = inputParser;
p.addParameter('metric','fold',@(s) any(validatestring(s,{'fold','percent','delta'})));
p.addParameter('errType','sem',@(s) any(validatestring(s,{'sd','sem'})));
p.addParameter('barColor',[0.30 0.60 0.85]);
p.addParameter('epsDen',1e-12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
metric   = p.Results.metric;
errType  = p.Results.errType;
barColor = p.Results.barColor;
epsDen   = p.Results.epsDen;

nEpoch = numel(epochRs);
if numel(epochNames) ~= nEpoch
    error('epochRs and epochNames must have same length.');
end

meanAcrossRats = nan(nEpoch,1);
errAcrossRats  = nan(nEpoch,1);

for e = 1:nEpoch
    R = epochRs{e};
    if isempty(R) || ~isstruct(R)
        continue;
    end
    nR = numel(R);
    perRat_fc = nan(nR,1);

    for rr = 1:nR
        tri = []; non = [];
        for d = 1:numel(R(rr).perDay)
            S = R(rr).perDay{d};
            if isempty(S), continue; end
            tested = isfinite(S.pVal);
            tri = [tri; S.meanTrialRate(tested)];
            non = [non; S.meanNonTrialRate(tested)];
        end
        if isempty(tri) || isempty(non), continue; end
        K = min(numel(tri), numel(non));
        tri = tri(1:K);
        non = non(1:K);

        denom = max(non, epsDen);
        switch metric
            case 'fold'
                fc = tri ./ denom;
            case 'percent'
                fc = 100 * (tri - non) ./ denom;
            otherwise  % 'delta'
                fc = tri - non;
        end
        fc = fc(isfinite(fc));
        if isempty(fc), continue; end
        perRat_fc(rr) = mean(fc,'omitnan');
    end

    perRat_fc = perRat_fc(isfinite(perRat_fc));
    if isempty(perRat_fc), continue; end

    meanAcrossRats(e) = mean(perRat_fc,'omitnan');
    s                 = std(perRat_fc,'omitnan');
    if strcmpi(errType,'sem')
        errAcrossRats(e) = s / sqrt(numel(perRat_fc));
    else
        errAcrossRats(e) = s;
    end
end

figure('Color','w','Position',[420 320 540 420]); hold on
bar(1:nEpoch, meanAcrossRats, 'FaceColor',barColor, 'EdgeColor','none');
errorbar(1:nEpoch, meanAcrossRats, errAcrossRats, errAcrossRats, '.k', ...
    'LineWidth',1.2,'CapSize',10);

% unity line for fold, zero for percent/delta
switch metric
    case 'fold'
        yline(1,'r--','LineWidth',1.2);
    otherwise
        yline(0,'k--','LineWidth',1.2);
end

xticks(1:nEpoch);
xticklabels(epochNames);
xtickangle(0);

switch metric
    case 'fold'
        ylabel('Cell-level fold change (Trial / Non-trial)');
    case 'percent'
        ylabel('Cell-level change from non-trial (%)');
    otherwise
        ylabel('Cell-level \Delta rate (Hz)  Trial - Non-trial');
end
title('Cell-level epoch-wise change (averaged across rats)');
box on

yl = ylim;
pad = 0.08*(yl(2)-yl(1) + eps);
ylim([yl(1)-pad*0.4, yl(2)+pad]);
end

function plot_epoch_foldchange_POP(epochRs, epochNames, varargin)
% plot_epoch_foldchange_POP
%   Epoch-wise *POPULATION* fold change bars.
%   Uses S.pop.trialRatePerBin / nonRatePerBin (speed-matched POP bins).
%
%   Inputs:
%     epochRs    : cell array of R structs (output of run_speedBinMatched)
%     epochNames : labels, e.g. {'CS','Trace','US','Post'}

p = inputParser;
p.addParameter('metric','fold',@(s) any(validatestring(s,{'fold','percent','delta'})));
p.addParameter('errType','sem',@(s) any(validatestring(s,{'sd','sem'})));
p.addParameter('barColor',[0.30 0.60 0.85]);
p.addParameter('epsDen',1e-12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
metric   = p.Results.metric;
errType  = p.Results.errType;
barColor = p.Results.barColor;
epsDen   = p.Results.epsDen;

nEpoch = numel(epochRs);
if numel(epochNames) ~= nEpoch
    error('epochRs and epochNames must have same length.');
end

meanAcrossRats = nan(nEpoch,1);
errAcrossRats  = nan(nEpoch,1);

for e = 1:nEpoch
    R = epochRs{e};
    if isempty(R) || ~isstruct(R)
        continue;
    end
    nR = numel(R);
    perRat_fc = nan(nR,1);

    for rr = 1:nR
        tri = []; non = [];
        for d = 1:numel(R(rr).perDay)
            S = R(rr).perDay{d};
            if isempty(S) || ~isfield(S,'pop') || isempty(S.pop)
                continue;
            end
            tri = [tri; S.pop.trialRatePerBin(:)];
            non = [non; S.pop.nonRatePerBin(:)];
        end
        tri = tri(isfinite(tri));
        non = non(isfinite(non));
        if isempty(tri) || isempty(non), continue; end

        K = min(numel(tri), numel(non));
        tri = tri(1:K);
        non = non(1:K);

        denomMean = max(mean(non,'omitnan'), epsDen);
        deltaMean = mean(tri,'omitnan') - mean(non,'omitnan');

        switch metric
            case 'fold'
                fc = mean(tri,'omitnan') / denomMean;
            case 'percent'
                fc = 100 * (deltaMean / denomMean);
            otherwise
                fc = deltaMean;
        end
        perRat_fc(rr) = fc;
    end

    perRat_fc = perRat_fc(isfinite(perRat_fc));
    if isempty(perRat_fc), continue; end

    meanAcrossRats(e) = mean(perRat_fc,'omitnan');
    s                 = std(perRat_fc,'omitnan');
    if strcmpi(errType,'sem')
        errAcrossRats(e) = s / sqrt(numel(perRat_fc));
    else
        errAcrossRats(e) = s;
    end
end

figure('Color','w','Position',[440 340 540 420]); hold on
bar(1:nEpoch, meanAcrossRats, 'FaceColor',barColor, 'EdgeColor','none');
errorbar(1:nEpoch, meanAcrossRats, errAcrossRats, errAcrossRats, '.k', ...
    'LineWidth',1.2,'CapSize',10);

switch metric
    case 'fold'
        yline(1,'r--','LineWidth',1.2);
    otherwise
        yline(0,'k--','LineWidth',1.2);
end

xticks(1:nEpoch);
xticklabels(epochNames);
xtickangle(0);

switch metric
    case 'fold'
        ylabel('POPULATION fold change (Trial / Non-trial)');
    case 'percent'
        ylabel('POPULATION change from non-trial (%)');
    otherwise
        ylabel('POPULATION \Delta rate (Hz)  Trial - Non-trial');
end
title('POPULATION epoch-wise change (averaged across rats)');
box on

yl = ylim;
pad = 0.08*(yl(2)-yl(1) + eps);
ylim([yl(1)-pad*0.4, yl(2)+pad]);
end
