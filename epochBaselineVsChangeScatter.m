function R = epochBaselineVsChangeScatter(ratNames, varargin)
% epochBaselineVsChangeScatter
% For each epoch (CS/Trace/US/Post), scatter baseline FR (x) vs change in FR (y).
%
% Usage:
%   R = epochBaselineVsChangeScatter({'rat0222','rat0307','rat0313','rat0314','rat0816'});
%
% Options (name/value):
%   'EpochEdges'   : [-2 0 0.25 0.75 0.85 2.00]   % baseline + 4 epochs (s)
%   'ChangeMetric' : 'delta'    % 'abs' | 'frac' | 'delta'
%                     abs   = FR_e - FR_b
%                     frac  = (FR_e - FR_b) ./ (FR_b + eps0)
%                     delta = (FR_e - FR_b) ./ (FR_e + FR_b + eps0)  (your Δ)
%   'MinTrials'    : 20         % min trials required
%   'MinBaseSpk'   : 5          % min total baseline spikes across trials
%   'MinTrialSpk'  : 5          % min total spikes in 0–2 s across trials (stability)
%   'Eps'          : 1e-5       % stabilizer for frac/delta
%   'XScale'       : 'log'      % 'linear' | 'log'  (log handles zeros via +epsX)
%   'Plot'         : true
%
% Returns:
%   R.table : table with columns: rat, FR_base, FR_CS, FR_Trace, FR_US, FR_Post,
%             d_CS, d_Trace, d_US, d_Post (change per chosen metric)
%   R.rho   : Spearman rho per epoch
%   R.r     : Pearson r per epoch
%   R.p     : p-values for both correlations per epoch


p = inputParser;
addParameter(p,'EpochEdges',[-2 0 0.25 0.75 0.85 2.00]);
addParameter(p,'ChangeMetric','frac',@(s) any(validatestring(s,{'abs','frac','delta'})));
addParameter(p,'MinTrials',20,@isscalar);
addParameter(p,'MinBaseSpk',5,@isscalar);
addParameter(p,'MinTrialSpk',5,@isscalar);
addParameter(p,'Eps',1e-5,@isscalar);
addParameter(p,'XScale','linear',@(s) any(validatestring(s,{'linear','log'})));
addParameter(p,'MinBaselineHzForPlot',0,@isscalar);   % hide only in plot
addParameter(p,'PlotStyle','deciles',@(s) any(validatestring(s,{'scatter','deciles','hex'})));
addParameter(p,'Plot',true,@islogical);

parse(p,varargin{:});
E        = p.Results.EpochEdges;
cmode    = lower(p.Results.ChangeMetric);
minTr    = p.Results.MinTrials;
minBs    = p.Results.MinBaseSpk;
minTrSpk = p.Results.MinTrialSpk;
eps0     = p.Results.Eps;
xscale   = p.Results.XScale;
minBasePlot = p.Results.MinBaselineHzForPlot;
style    = p.Results.PlotStyle;
doplot   = p.Results.Plot;

labels = {'CS','Trace','US','Post'}; nE = 4;

FR_base = []; FR_e = []; d_all = []; rat_id = [];

for r = 1:numel(ratNames)
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    iAn   = find(strcmp(dates,rat.An),1);
    days  = dates(max(1,iAn-2):iAn);

    for d = 1:numel(days)
        D  = days{d};
        S  = rat.Ca_peaks.(sprintf('CA_peaks_%s',D));
        CS = rat.CS_times.(sprintf('CS_%s',D));
        M  = rat.ratemask.(sprintf('ratemask_%s',D))==1;

        if numel(CS) < minTr, continue, end

        nC = size(S,1);
        for c = 1:nC
            if ~M(c), continue, end
            st = S(c,:); st = st(~isnan(st) & st>0);
            if isempty(st), continue, end

            fr = zeros(1,numel(E)-1);
            cntBaseline = 0; cntTrial = 0;

            for ee = 1:numel(E)-1
                cnt = 0;
                for t = 1:numel(CS)
                    t0 = CS(t) + E(ee);
                    t1 = CS(t) + E(ee+1);
                    cnt = cnt + sum(st>=t0 & st<t1);
                end
                dur    = (E(ee+1)-E(ee)) * numel(CS);
                fr(ee) = cnt / max(dur, eps);
                if ee==1, cntBaseline = cntBaseline + cnt;
                else,      cntTrial    = cntTrial    + cnt;
                end
            end

            % Skip if BOTHs sides are too sparse
            if cntBaseline < minBs && cntTrial < minTrSpk
                continue
            end

            b    = fr(1);
            evec = fr(2:end);
            switch cmode
                case 'abs',   dvec = evec - b;                          % Hz
                case 'frac',  dvec = (evec - b) ./ (b + eps0);           % unitless
                case 'delta', dvec = (evec - b) ./ (evec + b + eps0);    % unitless
            end

            FR_base = [FR_base; b];
            FR_e    = [FR_e;    evec];
            d_all   = [d_all;   dvec];
            rat_id  = [rat_id;  r];
        end
    end
end

% ----- output table -----
T = table;
T.rat      = rat_id;
T.FR_base  = FR_base;
T.FR_CS    = FR_e(:,1);
T.FR_Trace = FR_e(:,2);
T.FR_US    = FR_e(:,3);
T.FR_Post  = FR_e(:,4);
T.d_CS     = d_all(:,1);
T.d_Trace  = d_all(:,2);
T.d_US     = d_all(:,3);
T.d_Post   = d_all(:,4);

% ----- correlations -----
rho = nan(1,nE); p_rho = nan(1,nE);
rlin= nan(1,nE); p_lin = nan(1,nE);
X = T.FR_base; Xplot = X; epsX = 1e-9;
if strcmp(xscale,'log'), Xplot = X + epsX; end
for e = 1:nE
    y = T{:,sprintf('d_%s',labels{e})};
    [rho(e), p_rho(e)] = corr(X, y, 'Type','Spearman','Rows','complete');
    [rlin(e), p_lin(e)] = corr(X, y, 'Type','Pearson' ,'Rows','complete');
end

% ----- plotting -----
if doplot
    figure('Color','w','Position',[100 100 1200 420]);
    tl = tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
    for e = 1:nE
        nexttile; hold on
        y = T{:,sprintf('d_%s',labels{e})};


        show = isfinite(y) & Xplot >= max(minBasePlot, epsX);

        switch style
            case 'scatter'
                scatter(Xplot(show), y(show), 7, 'k', 'filled', ...
                        'MarkerFaceAlpha',0.18, 'MarkerEdgeAlpha',0.18);

            case 'hex'
                % 2D density (hex-like via histcounts2 grid)
                if strcmp(xscale,'log')
                    xx = log10(Xplot(show));
                    xlab = 'Baseline FR (Hz, log_{10})';
                else
                    xx = Xplot(show); xlab = 'Baseline FR (Hz)';
                end
                yy = y(show);
                nx = 60; ny = 60;
                [N,Xedges,Yedges] = histcounts2(xx,yy,[nx ny]);
                imagesc([Xedges(1) Xedges(end)],[Yedges(1) Yedges(end)], N');
                axis xy
                colormap(parula); alpha 0.9
                if strcmp(xscale,'log'), xticks(get(gca,'XTick')); ...
                        xticklabels(arrayfun(@(v)sprintf('%.2g',10.^v), get(gca,'XTick'),'uni',0)); end
                xlabel(xlab);

            case 'deciles'
                % background faint scatter for context
                scatter(Xplot(show), y(show), 10, [0 0 0], 'filled', ...
                        'MarkerFaceAlpha',0.12, 'MarkerEdgeAlpha',0.12);
                % decile medians ± 95% CI
                qs = prctile(Xplot(show), 0:10:100);
                xc = nan(1,10); med = xc; lo = xc; hi = xc;
                for k = 1:10
                    idx = Xplot(show) >= qs(k) & Xplot(show) < qs(k+1);
                    yi = y(show); yi = yi(idx);
                    if isempty(yi), continue, end
                    med(k) = median(yi,'omitnan');
                    if numel(yi) >= 10
                        bs = bootstrp(1000,@median,yi);
                        lo(k) = prctile(bs,2.5); hi(k) = prctile(bs,97.5);
                    end
                    if strcmp(xscale,'log')
                        xc(k) = 10.^mean(log10([qs(k)+epsX, qs(k+1)+epsX]));
                    else
                        xc(k) = mean([qs(k), qs(k+1)]);
                    end
                end
                errorbar(xc, med, med-lo, hi-med, 'o-', 'LineWidth',1.5, ...
                         'MarkerFaceColor',[0.2 0.4 0.9], 'Color',[0.2 0.4 0.9]);
        end

        if strcmp(xscale,'log'), set(gca,'XScale','log'); end
        xlabel('Baseline FR (Hz)');
        switch cmode
            case 'abs',   yl = 'Change (Hz) = FR_{epoch} - FR_{base}';
            case 'frac',  yl = 'Fractional change = (FR_{e}-FR_{b})/(FR_{b}+\epsilon)';
            case 'delta', yl = '\Delta = (FR_{e}-FR_{b})/(FR_{e}+FR_{b}+\epsilon)';
        end
        ylabel(yl);
        title(sprintf('%s  (\\rho=%.2f, p=%.1g; r=%.2f, p=%.1g)', ...
              labels{e}, rho(e), p_rho(e), rlin(e), p_lin(e)));
        box off
    end
    title(tl, sprintf('Baseline vs Epoch Change (%s, %s)', upper(cmode), style));
end

R.table  = T;
R.rho    = rho;   R.p_rho = p_rho;
R.r      = rlin;  R.p     = p_lin;
R.params = p.Results;
end
