function OUT = compareSpatialVsTemporalLag(ratNames, varargin)
% compareSpatialVsTemporalLag (fixed)
% Overlays spatial (populationSpatialOrthogonality) and temporal
% (populationOrthogonality_group) lag curves.
%
% Example:
% OUT = compareSpatialVsTemporalLag( ...
%   {'rat0222','rat0307','rat0313','rat0314','rat0816'}, ...
%   'NBins',12,'Similarity','corr','BinMode','equal_occ','Mode','demean');

% ---------- parse args ----------
p = inputParser;  p.KeepUnmatched = true;
addParameter(p,'Mode','raw',@(s) any(strcmpi(s,{'raw','demean','zscore'})));
addParameter(p,'Similarity','corr',@(s) any(strcmpi(s,{'corr','pearson','cosine','spearman'})));
addParameter(p,'NBins',16,@(x) isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'GridRC',[],@(v) isempty(v) || (isnumeric(v)&&numel(v)==2&&all(v>=1)));
addParameter(p,'BinMode','equal_size',@(s) any(strcmpi(s,{'equal_occ','equal_size'})));
addParameter(p,'MinSpeed',4,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Adjacency','moore',@(s) any(strcmpi(s,{'moore','rook'})));
addParameter(p,'Days',[],@(d) ischar(d) || isstring(d) || iscellstr(d) || isempty(d));
addParameter(p,'NSplits',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=2));  % default later to NBins
addParameter(p,'WinSecs',[0 2],@(v) isnumeric(v)&&numel(v)==2&&v(2)>v(1));
addParameter(p,'MICutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'MIExcludeCutoff',[],@(x) isempty(x) || (isscalar(x)&&isfinite(x)));
addParameter(p,'Spatial',struct(),@(s) isstruct(s));
addParameter(p,'Temporal',struct(),@(s) isstruct(s));
parse(p,varargin{:});
C = p.Results;

% ---- build SPATIAL cfg ----
sp = struct( ...
   'NBins',C.NBins,'GridRC',C.GridRC,'BinMode',lower(C.BinMode), ...
   'MinSpeed',C.MinSpeed,'Mode',lower(C.Mode), ...
   'Similarity',lower(C.Similarity),'Adjacency',lower(C.Adjacency), ...
   'Days',C.Days,'MICutoff',C.MICutoff,'MIExcludeCutoff',C.MIExcludeCutoff, ...
   'PerRatFigures',false,'SaveFig',false);
if ~isempty(C.Spatial), sp = mergeStructs(sp, C.Spatial); end

% ---- build TEMPORAL cfg ----
temSim    = mapTemporalSimilarity(lower( pick(C.Temporal,'Similarity', C.Similarity) ));
temSplits = pick(C.Temporal,'NSplits', C.NSplits);
if isempty(temSplits), temSplits = sp.NBins; end   % match NBins by default

tp = struct( ...
   'NSplits',temSplits,'WinSecs',pick(C.Temporal,'WinSecs', C.WinSecs), ...
   'Mode',lower(pick(C.Temporal,'Mode', C.Mode)), ...
   'Similarity',temSim, ...
   'MICutoff',pick(C.Temporal,'MICutoff', C.MICutoff), ...
   'MIExcludeCutoff',pick(C.Temporal,'MIExcludeCutoff', C.MIExcludeCutoff), ...
   'DoStats',true,'GroupSaveFig',false);
if ~isempty(C.Temporal)
    tp = mergeStructs(tp, C.Temporal);
    tp.Similarity = mapTemporalSimilarity(tp.Similarity);
end

% ---------- run SPATIAL ----------
Rsp = populationSpatialOrthogonality(ratNames, ...
    'NBins',sp.NBins,'GridRC',sp.GridRC,'BinMode',sp.BinMode, ...
    'MinSpeed',sp.MinSpeed,'Mode',sp.Mode,'Similarity',sp.Similarity, ...
    'Adjacency',sp.Adjacency,'Days',sp.Days, ...
    'MICutoff',sp.MICutoff,'MIExcludeCutoff',sp.MIExcludeCutoff, ...
    'PerRatFigures',sp.PerRatFigures,'SaveFig',false);

% extract spatial lag (vector by integer lag)
sLag = [];
if isfield(Rsp,'pooledAcrossRats') && isfield(Rsp.pooledAcrossRats,'spatialStats') ...
        && isfield(Rsp.pooledAcrossRats.spatialStats,'lag_mean')
    sLag = Rsp.pooledAcrossRats.spatialStats.lag_mean(:);
end
if isempty(sLag)
    Cpool = Rsp.pooledAcrossRats.cosSimK;
    sLag  = lagMean_from_matrix(Cpool);
end

% ------- FORCE zero-lag at start for spatial -------
% If the provided lag vector starts at lag 1, prepend the lag-0 anchor.
% (Similarity at zero spatial distance is 1 by definition.)
if isempty(sLag) || ~isfinite(sLag(1)) || sLag(1) < 0.999
    sLag = [1; sLag(:)];
end
xS = (0:numel(sLag)-1).';

% ---------- run TEMPORAL ----------
Gtp = populationOrthogonality_group(ratNames, ...
    'NSplits',tp.NSplits,'WinSecs',tp.WinSecs, ...
    'Mode',tp.Mode,'Similarity',tp.Similarity, ...
    'MICutoff',tp.MICutoff,'MIExcludeCutoff',tp.MIExcludeCutoff, ...
    'DoStats',true,'NBoot',500,'NPerm',500, ...
    'GroupSaveFig',false);

% temporal lag (group mean ± SEM)
if isfield(Gtp,'lagMat') && ~isempty(Gtp.lagMat)
    tLag = Gtp.lagMat;                   % nRats × NSplits
    yT   = mean(tLag,1,'omitnan').';
    yTse = std(tLag,0,1,'omitnan').' ./ sqrt(max(1,sum(isfinite(tLag),1)).');
else
    Cg   = Gtp.simMean;
    yT   = lagMean_from_matrix(Cg);
    yTse = nan(size(yT));
end
xT = (0:numel(yT)-1).';

% ---------- plot ----------
F = figure('Color','w','Position',[120 120 720 440]); hold on
if any(isfinite(yTse))
    xx = [xT; flipud(xT)];
    yy = [yT - yTse; flipud(yT + yTse)];
    fill(xx, yy, [0 0 0], 'FaceAlpha', 0.1, 'EdgeColor','none');
end
p1 = plot(xT, yT, '-', 'LineWidth', 2.4);              % temporal
p2 = plot(xS, sLag, '-', 'LineWidth', 2.4);            % spatial
xlabel('Lag (bins)');
if any(strcmpi(C.Similarity,{'corr','pearson'})) || strcmpi(tp.Similarity,'pearson')
    ylabel('Mean Pearson similarity');
else
    ylabel(sprintf('Mean %s similarity', lower(tp.Similarity)));
end
legend([p1 p2], {'Temporal (populationOrthogonality\_group)', ...
                 'Spatial (populationSpatialOrthogonality)'}, ...
       'Location','northeast');
title('Lag curves: temporal vs spatial (group summaries)');
grid on; box on
xlim([0, max([xT(:); xS(:)])]);       % start at zero for both

% ---------- outputs ----------
OUT = struct();
OUT.spatial  = struct('R',Rsp);
OUT.temporal = struct('G',Gtp);
OUT.curves   = struct('spatial',struct('x',xS,'y',sLag), ...
                      'temporal',struct('x',xT,'y',yT,'y_sem',yTse));
OUT.figure = F;

end

% ===== helpers =====
function S = mergeStructs(A,B)
S = A; if isempty(B), return; end
f = fieldnames(B); for i=1:numel(f), S.(f{i}) = B.(f{i}); end
end

function v = pick(S, field, defaultVal)
if isfield(S,field) && ~isempty(S.(field)), v = S.(field); else, v = defaultVal; end
end

function s = mapTemporalSimilarity(s)
% temporal function accepts: 'cosine','pearson','spearman'
if any(strcmpi(s,{'corr','pearson'})), s = 'pearson'; end
if any(strcmpi(s,{'cos'})),            s = 'cosine';  end
% leave 'spearman' as-is
end

function y = lagMean_from_matrix(C)
K = size(C,1); y = nan(K,1);
for L = 0:(K-1)
    idx = diag(true(K-L,1), L) | diag(true(K-L,1), -L);
    y(L+1) = mean(C(idx), 'omitnan');
end
end
