function F = MI_distribution2(RatList)
% MI_distribution_withSpeed  (BASELINE = WITH TASK / ALL CELLS)
% Same as MI_distribution, but ALSO includes SPEED-MATCHED REMOVED:
%   D_speed -> MI_noCSUS15_controlSpeed_all(:,2)
%
% Outputs:
%   F.equal_weight_pooled, F.perrat_overlay

if nargin < 1 || isempty(RatList)
    RatList = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end
RatNames = strrep(RatList,'rat','Rat ');

% -------- pull per-cell arrays using SAME mask/columns as your bar code --------
A_with  = cell(1,numel(RatList));   % WITH TASK      -> MI_wCSUS_all
B_no    = cell(1,numel(RatList));   % NO TASK        -> MI_noCSUS15_all
C_rand  = cell(1,numel(RatList));   % RAND REMOVED   -> MI_noCSUS15_control_all(:,2)
D_speed = cell(1,numel(RatList));   % SPEED MATCHED  -> MI_noCSUS15_controlSpeed_all(:,2)
Ns      = zeros(1,numel(RatList));

for i = 1:numel(RatList)
    varname = RatList{i};
    if ~evalin('base', sprintf('exist(''%s'',''var'')', varname))
        error('Variable %s not found in base workspace.', varname);
    end
    R = evalin('base', varname);

    need = {'MI_wCSUS_all','MI_noCSUS15_all','MI_noCSUS15_control_all','MI_noCSUS15_controlSpeed_all','ratemask_all'};
    missing = need(~cellfun(@(f)isfield(R,f), need));
    if ~isempty(missing)
        error('Fields missing in %s: %s', varname, strjoin(missing, ', '));
    end

    m_rows = (R.ratemask_all==1);                    % same mask as bar plot
    A  = R.MI_wCSUS_all(m_rows);                     % WITH TASK  (baseline)
    B  = R.MI_noCSUS15_all(m_rows);                  % NO TASK
    C  = R.MI_noCSUS15_control_all(m_rows,2);         % RAND REMOVED (your original)
    D  = R.MI_noCSUS15_controlSpeed_all(m_rows,2);    % SPEED-MATCHED REMOVED (new)

    A_with{i}  = A(:);
    B_no{i}    = B(:);
    C_rand{i}  = C(:);
    D_speed{i} = D(:);
    Ns(i)      = numel(A);
end

% -------- (A) sanity table: means and SEMs --------
fprintf('\n=== Sanity check (WITH TASK, NO TASK, RANDOM REMOVED, SPEED-MATCHED REMOVED) ===\n');
fprintf('Rat\t\t n\t mean(WITH)\t SEM(WITH)\t mean(NO)\t SEM(NO)\t mean(RAND)\t SEM(RAND)\t mean(SPEED)\t SEM(SPEED)\n');
for i = 1:numel(RatList)
    n  = Ns(i);
    mW = nanmean(A_with{i});  sW = nanstd(A_with{i})  / sqrt(n);
    mN = nanmean(B_no{i});    sN = nanstd(B_no{i})    / sqrt(n);
    mR = nanmean(C_rand{i});  sR = nanstd(C_rand{i})  / sqrt(n);
    mS = nanmean(D_speed{i}); sS = nanstd(D_speed{i}) / sqrt(n);
    fprintf('%-8s\t %4d\t % .6f\t % .6f\t % .6f\t % .6f\t % .6f\t % .6f\t % .6f\t % .6f\n', ...
        RatNames{i}, n, mW, sW, mN, sN, mR, sR, mS, sS);
end

% -------- Stats: paired t-tests & KS tests (per rat and pooled) --------
fprintf('\n=== Stats: WITH TASK vs NO TASK (paired across cells) ===\n');
print_pair_stats(RatNames, A_with, B_no, 'WithTask', 'NoTask');

fprintf('\n=== Stats: WITH TASK vs RANDOM REMOVED (paired across cells) ===\n');
print_pair_stats(RatNames, A_with, C_rand, 'WithTask', 'RandRemoved');

fprintf('\n=== Stats: WITH TASK vs SPEED-MATCHED REMOVED (paired across cells) ===\n');
print_pair_stats(RatNames, A_with, D_speed, 'WithTask', 'SpeedRemoved');

% -------- build %Δ (relative to WITH TASK) per rat --------
Pct_no_vs_with    = cellfun(@(B,A) 100*(B - A) ./ max(abs(A),eps), B_no,    A_with, 'uni',0);
Pct_rand_vs_with  = cellfun(@(C,A) 100*(C - A) ./ max(abs(A),eps), C_rand,  A_with, 'uni',0);
Pct_speed_vs_with = cellfun(@(D,A) 100*(D - A) ./ max(abs(A),eps), D_speed, A_with, 'uni',0);

for i = 1:numel(RatList)
    Pct_no_vs_with{i}    = Pct_no_vs_with{i}(isfinite(Pct_no_vs_with{i}));
    Pct_rand_vs_with{i}  = Pct_rand_vs_with{i}(isfinite(Pct_rand_vs_with{i}));
    Pct_speed_vs_with{i} = Pct_speed_vs_with{i}(isfinite(Pct_speed_vs_with{i}));
end

% shared binning from all rats (percent units)
all1 = vertcat(Pct_no_vs_with{:});
all2 = vertcat(Pct_rand_vs_with{:});
all3 = vertcat(Pct_speed_vs_with{:});
Edges    = pickEdgesPercent3(all1, all2, all3);
xcenters = (Edges(1:end-1)+Edges(2:end))/2;

% -------- (B) equal-weight pooled density with 95% CI + dashed mean lines --------
col_no    = [0.22 0.45 0.70];  % NO vs WITH (blue)
col_rand  = [0.85 0.33 0.10];  % RAND vs WITH (orange)
col_speed = [0.20 0.70 0.35];  % SPEED vs WITH (green)

H1 = cellfun(@(v) histcounts(v,Edges,'Normalization','probability'), Pct_no_vs_with,    'uni',0);
H2 = cellfun(@(v) histcounts(v,Edges,'Normalization','probability'), Pct_rand_vs_with,  'uni',0);
H3 = cellfun(@(v) histcounts(v,Edges,'Normalization','probability'), Pct_speed_vs_with, 'uni',0);
H1 = vertcat(H1{:});  H2 = vertcat(H2{:});  H3 = vertcat(H3{:});

mu1 = mean(H1,1,'omitnan');   % equal weight across rats
mu2 = mean(H2,1,'omitnan');
mu3 = mean(H3,1,'omitnan');

B = 1000;   % bootstrap over rats
boot1 = zeros(B,numel(mu1)); boot2 = zeros(B,numel(mu2)); boot3 = zeros(B,numel(mu3));
nR = size(H1,1);
for b = 1:B
    idx = randi(nR, nR, 1);
    boot1(b,:) = mean(H1(idx,:),1,'omitnan');
    boot2(b,:) = mean(H2(idx,:),1,'omitnan');
    boot3(b,:) = mean(H3(idx,:),1,'omitnan');
end
ci1 = prctile(boot1,[2.5 97.5],1);
ci2 = prctile(boot2,[2.5 97.5],1);
ci3 = prctile(boot3,[2.5 97.5],1);

% equal-weight MEANS of %Δ across rats (mean of per-rat means)
eqMean_no    = mean(cellfun(@(v) mean(v,'omitnan'), Pct_no_vs_with),    'omitnan');
eqMean_rand  = mean(cellfun(@(v) mean(v,'omitnan'), Pct_rand_vs_with),  'omitnan');
eqMean_speed = mean(cellfun(@(v) mean(v,'omitnan'), Pct_speed_vs_with), 'omitnan');

F.equal_weight_pooled = figure('Color','w','Name','Equal-weight pooled %Δ relative to WITH TASK');
hold on;
fill([xcenters fliplr(xcenters)], [ci1(1,:) fliplr(ci1(2,:))], col_no,    'FaceAlpha',0.15, 'EdgeColor','none');
fill([xcenters fliplr(xcenters)], [ci2(1,:) fliplr(ci2(2,:))], col_rand,  'FaceAlpha',0.15, 'EdgeColor','none');
fill([xcenters fliplr(xcenters)], [ci3(1,:) fliplr(ci3(2,:))], col_speed, 'FaceAlpha',0.15, 'EdgeColor','none');

plot(xcenters, mu1, 'LineWidth',2, 'Color', col_no);
plot(xcenters, mu2, 'LineWidth',2, 'Color', col_rand);
plot(xcenters, mu3, 'LineWidth',2, 'Color', col_speed);

xline(eqMean_no,    '--','LineWidth',1.8,'Color',col_no);
xline(eqMean_rand,  '--','LineWidth',1.8,'Color',col_rand);
xline(eqMean_speed, '--','LineWidth',1.8,'Color',col_speed);
xline(0,'k-','LineWidth',1);

ymax = max([mu1,mu2,mu3,ci1(:)',ci2(:)',ci3(:)']);
if ~isfinite(ymax) || ymax<=0, ymax=0.05; end
ylim([0, ymax*1.05]);
ytickformat('%.2f%%');
xlabel('% change in MI relative to WITH TASK (all cells)');
ylabel('Proportion of cells (equal-weight across rats)');
legend({'No-task 95% CI','Random-removed 95% CI','Speed-removed 95% CI', ...
        'No-task mean density','Random-removed mean density','Speed-removed mean density', ...
        'No-task mean %Δ','Random-removed mean %Δ','Speed-removed mean %Δ'}, 'Location','northwest');
title('Equal-weight pooled densities (baseline = WITH TASK / all cells)');
set(gca,'Box','off','LineWidth',1);
hold off;

% -------- (C) per-rat overlaid histograms (shared bins, proportion y-axis) --------
nCols = min(3, nR);
nRows = ceil(nR / nCols);
F.perrat_overlay = figure('Color','w','Name','Per-rat %Δ relative to WITH TASK (overlaid)');
tlo = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
for r = 1:nR
    nexttile(tlo,r);
    pn = Pct_no_vs_with{r};   pr = Pct_rand_vs_with{r};   ps = Pct_speed_vs_with{r};
    hold on;
    histogram(pn, Edges, 'Normalization','probability', 'FaceColor',col_no,   'EdgeColor','none','FaceAlpha',0.45);
    histogram(pr, Edges, 'Normalization','probability', 'FaceColor',col_rand, 'EdgeColor','none','FaceAlpha',0.45);
    histogram(ps, Edges, 'Normalization','probability', 'FaceColor',col_speed,'EdgeColor','none','FaceAlpha',0.45);
    xline(median(pn), '-','LineWidth',1.2,'Color',col_no);
    xline(median(pr), '-','LineWidth',1.2,'Color',col_rand);
    xline(median(ps), '-','LineWidth',1.2,'Color',col_speed);
    xline(0,'k-','LineWidth',1);
    xlim([Edges(1) Edges(end)]);
    ytickformat('%.2f%%');
    title(RatNames{r}, 'Interpreter','none');
    if r > (nRows-1)*nCols, xlabel('% Δ relative to WITH TASK'); end
    if mod(r-1,nCols)==0,   ylabel('Proportion of cells'); end
    set(gca,'Box','off','LineWidth',1);
    hold off;
end
sgtitle(tlo, 'Per-rat % change relative to WITH TASK (overlaid)');

end  % ===== end main =====


% ---------- helpers ----------
function Edges = pickEdgesPercent3(D1, D2, D3)
D = [D1(:); D2(:); D3(:)];
D = D(isfinite(D));
if isempty(D), Edges = -10:1:10; return; end
lohi = prctile(D,[1 99]); pad = 0.05*diff(lohi);
x1 = lohi(1)-pad; x2 = lohi(2)+pad;
nb = max(40, min(60, round((x2-x1)/max(0.5, 0.01*(x2-x1)))));
bw = (x2-x1)/nb; Edges = x1:bw:x2;
end

function print_pair_stats(RatNames, Xcell, Ycell, Xlab, Ylab)
R = numel(Xcell);
fprintf('Per-rat (paired t-test of %s - %s = 0; KS two-sample on raw arrays)\n', Xlab, Ylab);
fprintf('%-8s  %-6s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Rat','n_pairs','t(paired)','p_t','KS stat','p_KS','meanΔ');
allDiff = [];
for i = 1:R
    x = Xcell{i}; y = Ycell{i};
    m = isfinite(x) & isfinite(y);
    n = nnz(m);
    tstat = NaN; pt = NaN; D = NaN; pKS = NaN; md = NaN;
    if n >= 2
        d = x(m) - y(m);
        md = mean(d,'omitnan');
        [~, pt, ~, stats] = ttest(x(m), y(m));  % paired t-test
        tstat = stats.tstat;
        allDiff = [allDiff; d(:)]; %#ok<AGROW>
    end
    xf = x(isfinite(x)); yf = y(isfinite(y));
    if numel(xf) >= 2 && numel(yf) >= 2
        [~, pKS, ksstat] = kstest2(xf, yf);
        D = ksstat;
    end
    fprintf('%-8s  %6d  % -10.4f  % -10.3g  % -10.4f  % -10.3g  % -10.4f\n', ...
        RatNames{i}, n, tstat, pt, D, pKS, md);
end
tstat = NaN; pt = NaN; D = NaN; pKS = NaN; md = NaN;
if ~isempty(allDiff)
    md = mean(allDiff,'omitnan');
    [~, pt, ~, stats] = ttest(allDiff, 0);
    tstat = stats.tstat;
end
xf_all = []; yf_all = [];
for i = 1:R
    xf_all = [xf_all; Xcell{i}(isfinite(Xcell{i}))]; %#ok<AGROW>
    yf_all = [yf_all; Ycell{i}(isfinite(Ycell{i}))]; %#ok<AGROW>
end
if numel(xf_all) >= 2 && numel(yf_all) >= 2
    [~, pKS, ksstat] = kstest2(xf_all, yf_all);
    D = ksstat;
end
fprintf('POOLED    %6d  % -10.4f  % -10.3g  % -10.4f  % -10.3g  % -10.4f\n', ...
    numel(allDiff), tstat, pt, D, pKS, md);
end
