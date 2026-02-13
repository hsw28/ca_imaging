function OUT = run_space_to_task_interference(ratNames, varargin)

OUT = struct();
OUT.PV = run_space_to_task_interference_PV_splitHalf_fast(ratNames, varargin{:});
fprintf('done with PV')
OUT.SC = run_space_to_task_interference_SC_splitHalf_fast(ratNames, varargin{:});

end

function ATTS = run_space_to_task_interference_PV_splitHalf_fast(ratNames, varargin)

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'TimeBins',15);
addParameter(p,'GridN',[2 2]);
addParameter(p,'UseSpeedMask',false);
addParameter(p,'VelThresh',4);
addParameter(p,'CellNorm','demean');
addParameter(p,'MinFrames',1);
addParameter(p,'MinTrialsPerK',2);
addParameter(p,'NResamp',50);
addParameter(p,'NPerm',100);
addParameter(p,'MaxPairsPV',100);
addParameter(p,'Seed',1);
parse(p,varargin{:});
o = p.Results;

rng(o.Seed);

if ischar(ratNames) || isstring(ratNames)
    ratNames = cellstr(ratNames);
end

nR = numel(ratNames);
B  = o.TimeBins;

ATTS = struct();
ATTS.params = o;
ATTS.rats   = repmat(struct('name','','days',[],'delta',struct()), nR, 1);

for ri = 1:nR
    rname = ratNames{ri}
    RAT   = evalin('base', rname);
    ATTS.rats(ri).name = rname;

    dates = autoDateList(RAT);
    idx = find(strcmp(dates, RAT.An), 1);
    if isempty(idx), dayIdx = max(1,numel(dates)-2):numel(dates);
    else, dayIdx = max(1,idx-2):idx;
    end
    days = dates(dayIdx);

    ATTS.rats(ri).days = repmat(struct('label','', ...
        'obs',struct('metricZ',nan(3,1)), ...
        'perm',struct('deltaZ',[])), numel(days),1);

    for di = 1:numel(days)
        dlabel = days{di}
        ATTS.rats(ri).days(di).label = dlabel;

        spk  = RAT.Ca_peaks.(sprintf('CA_peaks_%s',dlabel));
        posM = (RAT.pos.(sprintf('pos_%s',dlabel)));

        cs   = coerce_cs(RAT.CS_times.(sprintf('CS_%s',dlabel)));

        if isfield(RAT,'ratemask')
            rm = RAT.ratemask.(sprintf('ratemask_%s',dlabel));
            spk = spk(rm==1,:);
        end

        [ts,x,y] = coerce_pos(posM);
        if isempty(cs) || size(spk,1)<2, continue; end

        spikeCell = to_cell_spikes(spk);
        S = normalize_cells(spikes_to_matrix(spikeCell,ts), o.CellNorm);

        if o.UseSpeedMask
            v = local_speed(ts,x,y);
            speed_ok = v>=o.VelThresh;
        else
            speed_ok = true(size(ts));
        end

        nx=o.GridN(1); ny=o.GridN(2);
        bX = discretize(x, linspace(min(x),max(x),nx+1));
        bY = discretize(y, linspace(min(y),max(y),ny+1));
        bSpace = nan(size(x));
        m = isfinite(bX)&isfinite(bY);
        bSpace(m)=sub2ind([nx ny],bX(m),bY(m));
        K = nx*ny;

        % -------- BUILD PV(tr,b,k) --------
        PV = cell(numel(cs),B,K);
        PVmask = false(numel(cs),B,K);

        for tr=1:numel(cs)
            idx_all = find(ts>=cs(tr)+o.TraceWin(1) & ts<cs(tr)+o.TraceWin(2));
            if isempty(idx_all), continue; end
            e = round(linspace(1,numel(idx_all)+1,B+1));

            for b=1:B
                seg = idx_all(e(b):e(b+1)-1);
                if o.UseSpeedMask, seg=seg(speed_ok(seg)); end
                if isempty(seg), continue; end
                ks = unique(bSpace(seg));
                ks = ks(isfinite(ks));
                for k=ks'
                    f = seg(bSpace(seg)==k);
                    if numel(f)<o.MinFrames, continue; end
                    PV{tr,b,k}=mean(S(:,f),2,'omitnan');
                    PVmask(tr,b,k)=true;
                end
            end
        end

        [obsZ, permZ] = pv_metrics_fast(PV, PVmask, o);
        ATTS.rats(ri).days(di).obs.metricZ  = obsZ;
        ATTS.rats(ri).days(di).perm.deltaZ  = permZ;
    end

    ATTS.rats(ri).delta = aggregate_delta_across_days(ATTS.rats(ri).days);
end

plot_spaceTask_splitHalf_bars_group(ATTS,'PV');
plot_spaceTask_splitHalf_permHist_perRat(ATTS,'PV');

end

function ATTSsc = run_space_to_task_interference_SC_splitHalf_fast(ratNames, varargin)

p = inputParser;
addParameter(p,'TraceWin',[0 2]);
addParameter(p,'TimeBins',15);
addParameter(p,'GridN',[2 2]);
addParameter(p,'CellNorm','none');
addParameter(p,'MinFrames',2);
addParameter(p,'MinTrialsPerK',2);
addParameter(p,'MinBinsCorr',2);
addParameter(p,'NResamp',50);
addParameter(p,'NPerm',100);
addParameter(p,'MaxPairsSC',500);
addParameter(p,'Seed',1);
parse(p,varargin{:});
o = p.Results;

rng(o.Seed);

if ischar(ratNames)||isstring(ratNames)
    ratNames=cellstr(ratNames);
end

ATTSsc = struct();
ATTSsc.params=o;
ATTSsc.rats=repmat(struct('name','','days',[],'delta',struct()),numel(ratNames),1);

for ri=1:numel(ratNames)
    rname=ratNames{ri}
    RAT=evalin('base',rname);
    ATTSsc.rats(ri).name=rname;

    dates=autoDateList(RAT);
    idx=find(strcmp(dates,RAT.An),1);
    if isempty(idx), dayIdx=max(1,numel(dates)-2):numel(dates);
    else, dayIdx=max(1,idx-2):idx;
    end
    days=dates(dayIdx);

    ATTSsc.rats(ri).days=repmat(struct('label','', ...
        'obs',struct('metricZ',nan(3,1)), ...
        'perm',struct('deltaZ',[])), numel(days),1);

    for di=1:numel(days)
        dlabel=days{di}
        ATTSsc.rats(ri).days(di).label=dlabel;

        spk=RAT.Ca_peaks.(sprintf('CA_peaks_%s',dlabel));
        %posM=smoothpos(RAT.pos.(sprintf('pos_%s',dlabel)));
        posM=(RAT.pos.(sprintf('pos_%s',dlabel)));

        cs=coerce_cs(RAT.CS_times.(sprintf('CS_%s',dlabel)));

        if isfield(RAT,'ratemask')
            rm=RAT.ratemask.(sprintf('ratemask_%s',dlabel));
            spk=spk(rm==1,:);
        end

        [ts,x,y]=coerce_pos(posM);
        if isempty(cs)||size(spk,1)<2, continue; end

        S=normalize_cells(spikes_to_matrix(to_cell_spikes(spk),ts),o.CellNorm);

        nx=o.GridN(1); ny=o.GridN(2);
        bX=discretize(x,linspace(min(x),max(x),nx+1));
        bY=discretize(y,linspace(min(y),max(y),ny+1));
        bSpace=nan(size(x));
        m=isfinite(bX)&isfinite(bY);
        bSpace(m)=sub2ind([nx ny],bX(m),bY(m));
        K=nx*ny;

        B=o.TimeBins; nT=numel(cs); nC=size(S,1);
        V=cell(nT,K); Vmask=false(nT,K);

        for tr=1:nT
            idx_all=find(ts>=cs(tr)+o.TraceWin(1) & ts<cs(tr)+o.TraceWin(2));
            if isempty(idx_all), continue; end
            e=round(linspace(1,numel(idx_all)+1,B+1));
            for b=1:B
                seg=idx_all(e(b):e(b+1)-1);
                if isempty(seg), continue; end
                ks=unique(bSpace(seg));
                ks=ks(isfinite(ks));
                for k=ks'
                    f=seg(bSpace(seg)==k);
                    if numel(f)<o.MinFrames, continue; end
                    if isempty(V{tr,k}), V{tr,k}=nan(nC,B); end
                    V{tr,k}(:,b)=mean(S(:,f),2,'omitnan');
                    Vmask(tr,k)=true;
                end
            end
        end

        [obsZ,permZ]=sc_metrics_fast(V,Vmask,o);
        ATTSsc.rats(ri).days(di).obs.metricZ=obsZ;
        ATTSsc.rats(ri).days(di).perm.deltaZ=permZ;
    end

    ATTSsc.rats(ri).delta=aggregate_delta_across_days(ATTSsc.rats(ri).days);
end

plot_spaceTask_splitHalf_bars_group(ATTSsc,'SC');
plot_spaceTask_splitHalf_permHist_perRat(ATTSsc,'SC');

end


function [ts,x,y] = coerce_pos(pos)
if istable(pos)
    vn = lower(string(pos.Properties.VariableNames));
    tname = pick(vn, ["t","time","ts"]);
    xname = pick(vn, ["x","xpos","x_cm","xsmooth","posx","x_smooth","xcm"]);
    yname = pick(vn, ["y","ypos","y_cm","ysmooth","posy","y_smooth","ycm"]);
    ts = getcol(pos, tname); x = getcol(pos, xname); y = getcol(pos, yname);
else
    pos = double(pos);
    if size(pos,2) >= 3
        ts = pos(:,1); x = pos(:,2); y = pos(:,3);
    elseif size(pos,2) == 2
        x = pos(:,1); y = pos(:,2); ts = (1:size(pos,1))';
    else
        error('pos must be [n×3] or [n×2]');
    end
end
ts = double(ts(:)); x = double(x(:)); y = double(y(:));
end

function v = local_speed(ts,x,y)
dt = diff(ts); dt(~isfinite(dt) | dt<=0) = NaN;
dx = diff(x);  dy = diff(y);
spd = sqrt(dx.^2 + dy.^2) ./ dt;
v = [spd; spd(end)];
v(~isfinite(v)) = 0;
end

function cs = coerce_cs(csIn)
if istable(csIn)
    vn = lower(string(csIn.Properties.VariableNames));
    cname = pick(vn, ["cs_ms","cs_time","cstime","cs","onset","onsets","cs_onset","cs_time_ms","cue_onset","time","ts"]);
    cs = getcol(csIn, cname);
else
    cs = double(csIn(:));
end
if max(cs,[],'omitnan') > 1e4, cs = cs/1000; end
end

function cellSpk = to_cell_spikes(spk)
if iscell(spk), cellSpk = spk; return; end
[nc, ~] = size(spk);
cellSpk = cell(nc,1);
for c = 1:nc
    st = double(spk(c,:));
    st = st(isfinite(st) & st > 0);
    cellSpk{c} = st(:);
end
end

function S = spikes_to_matrix(spikeCell, ts)
ts = double(ts(:));
if numel(ts) < 2
    S = zeros(numel(spikeCell), numel(ts), 'single'); return;
end
dt = median(diff(ts),'omitnan');
if ~isfinite(dt) || dt<=0, dt = max(eps, mean(diff(ts),'omitnan')); end
edges = [ts - dt/2; ts(end) + dt/2];
nc = numel(spikeCell);
S  = zeros(nc, numel(ts), 'single');
for c = 1:nc
    st = spikeCell{c};
    if isempty(st), continue; end
    S(c,:) = histcounts(st, edges);
end
end

function S = normalize_cells(S, mode)
switch lower(mode)
    case 'none'
    case 'demean'
        mu = mean(S,2,'omitnan');
        S = S - mu;
    case 'zscore'
        mu = mean(S,2,'omitnan'); sd = std(S,0,2,'omitnan'); sd(sd==0|~isfinite(sd)) = 1; S = (S - mu) ./ sd;
    case 'meanrate'
      mu = mean(S,2,'omitnan');
      S = S ./ mu;
    otherwise
        error('CellNorm must be ''none'',''demean'',''zscore'',''meanrate''.');
end
end

function [r, nC] = safe_corr(a,b,minN,minStd)
a = a(:); b = b(:);
if numel(a) ~= numel(b)
    n = min(numel(a), numel(b)); a = a(1:n); b = b(1:n);
end
if nargin < 3 || isempty(minN),   minN = 3;       end
if nargin < 4 || isempty(minStd), minStd = 1e-12; end
m = isfinite(a) & isfinite(b);
nC = nnz(m);
if nC < minN, r = NaN; return; end
a = a(m); b = b(m);
sa = std(a); sb = std(b);
if ~isfinite(sa) || ~isfinite(sb) || sa < minStd || sb < minStd
    r = 0; return;
end
r = corr(a,b,'Type','Pearson');
end

function v = getcol(T, name)
if strlength(name)==0, v = []; return; end
v = T{:, find(strcmpi(T.Properties.VariableNames, char(name)),1)};
end

function name = pick(names, options)
name = "";
for k = 1:numel(options)
    idx = find(names == options(k), 1);
    if ~isempty(idx), name = names(idx); return; end
end
end

function [p, t, df, dz] = paired_t_z(z1, z2)
mask = isfinite(z1) & isfinite(z2);
z1 = z1(mask); z2 = z2(mask);
if numel(z1) < 2
    p = NaN; t = NaN; df = 0; dz = NaN; return;
end
[~, p, ~, st] = ttest(z1, z2);  % two-sided
t  = st.tstat;
df = st.df;
d  = z1 - z2;
dz = mean(d,'omitnan') / std(d, 0, 'omitnan');
end

function r = corr_temporal_vec(v1, v2, minBins)
v1 = v1(:); v2 = v2(:);
m = isfinite(v1) & isfinite(v2);
if nnz(m) < minBins, r = NaN; return; end
a = v1(m); b = v2(m);
if std(a) < 1e-12 || std(b) < 1e-12, r = 0; return; end
r = corr(a,b,'Type','Pearson');
end

% ======================================================================
% ======================= MISSING CORE ENGINES =========================
% ======================================================================

function z = pv_inin_z(A,B,ks,o)
% IN-IN per k
z = nan(numel(ks),1);
ii = 0;
for j = 1:numel(ks)
    k = ks(j);
    if isempty(A{k}) || isempty(B{k}), continue; end
    [r,~] = safe_corr(A{k},B{k},3,1e-12);
    if ~isfinite(r), continue; end
    ii = ii+1;
    z(ii) = atanh(clamp_r(r));
end
z = z(1:ii);
end


function z = pv_crossk_z(A,B,ks,nTarget,o,mode)
% OUT-OUT or IN-OUT: sample random (k1,k2), k1!=k2, matched count
% mode is only semantic label here; both use A{k1} vs B{k2}
MaxPairs = o.MaxPairsPV;
nK = numel(ks);
if nK < 2, z = []; return; end

nDraw = min(MaxPairs, max(10*nTarget, 200));
z = nan(min(nTarget,nDraw),1);

ii = 0;
for t = 1:nDraw
    k1 = ks(randi(nK));
    k2 = ks(randi(nK));
    if k2==k1, continue; end
    if isempty(A{k1}) || isempty(B{k2}), continue; end
    [r,~] = safe_corr(A{k1},B{k2},3,1e-12);
    if ~isfinite(r), continue; end
    ii = ii+1;
    z(ii) = atanh(clamp_r(r));
    if ii >= nTarget, break; end
end
z = z(1:ii);
end


function [obsZ, permDeltaZ] = sc_metrics_fast(V, Vmask, o)
% OPTION B IMPLEMENTATION (SC):
% Build TA{k} and TB{k} for each valid k (split halves).
% Create derangement mapping f over ks.
% For each k and cell c:
%   ININ   = corr(TA{k}(c,:),   TB{k}(c,:))
%   OUTOUT = corr(TA{f(k)}(c,:),TB{k}(c,:))
%   INOUT  = corr(TA{k}(c,:),   TB{f(k)}(c,:))
%
% obsZ: [3x1] mean Fisher-z across resamples for [ININ OUTOUT INOUT]
% permDeltaZ: [NPermx1] null for Δ = mean(ININ)-mean(OUTOUT) using label-shuffle of TB across ks.

[nT,K]   = size(Vmask);
NResamp  = o.NResamp;
NPerm    = o.NPerm;
MinT     = o.MinTrialsPerK;
MinBins  = o.MinBinsCorr;

% infer nC,B
tmp = [];
for tr=1:nT
    ks0 = find(Vmask(tr,:));
    if ~isempty(ks0)
        tmp = V{tr,ks0(1)};
        break;
    end
end
if isempty(tmp)
    obsZ = nan(3,1);
    permDeltaZ = [];
    return;
end
[nC,B] = size(tmp); %#ok<ASGLU>

% available trials per k
trsByK = cell(K,1);
for k=1:K
    trsByK{k} = find(Vmask(:,k));
end

% cache templates per resample: TA/TB and ks
CACHE = repmat(struct('ks',[],'TA',[],'TB',[],'ok',false), NResamp, 1);

for rr=1:NResamp
    TA = cell(K,1);
    TB = cell(K,1);
    okK = false(K,1);

    for k=1:K
        trs = trsByK{k};
        if numel(trs) < MinT, continue; end
        [hA,hB] = split_half(trs);
        if isempty(hA) || isempty(hB), continue; end

        TA{k} = mean_cat_V(V,hA,k,nC,B);
        TB{k} = mean_cat_V(V,hB,k,nC,B);
        if isempty(TA{k}) || isempty(TB{k}), continue; end
        okK(k)=true;
    end

    ks = find(okK);
    if numel(ks) < 2
        CACHE(rr).ok = false;
    else
        CACHE(rr).ok = true;
        CACHE(rr).ks = ks;
        CACHE(rr).TA = TA;
        CACHE(rr).TB = TB;
    end
end

% ---------- observed ----------
z_res = nan(NResamp,3);

for rr=1:NResamp
    if ~CACHE(rr).ok, continue; end
    ks = CACHE(rr).ks;
    TA = CACHE(rr).TA;
    TB = CACHE(rr).TB;

    fks = derangement_map(ks);

    [z_in, z_out, z_inout] = sc_map_metrics_z(TA, TB, ks, fks, MinBins);

    if isempty(z_in), continue; end
    z_res(rr,1)=mean(z_in,'omitnan');
    z_res(rr,2)=mean(z_out,'omitnan');
    z_res(rr,3)=mean(z_inout,'omitnan');
end

obsZ = mean(z_res,1,'omitnan')';
obsZ = obsZ(:);

% ---------- perm null ----------
permDeltaZ = nan(NPerm,1);
if NPerm < 5 || all(~isfinite(obsZ)), return; end

for pidx=1:NPerm
    deltas = nan(NResamp,1);

    for rr=1:NResamp
        if ~CACHE(rr).ok, continue; end
        ks = CACHE(rr).ks;
        TA = CACHE(rr).TA;
        TB = CACHE(rr).TB;

        % shuffle TB labels across ks
        ks_shuf = ks(randperm(numel(ks)));
        TBsh = cell(K,1);
        for ii=1:numel(ks)
            TBsh{ks(ii)} = TB{ks_shuf(ii)};
        end

        fks = derangement_map(ks);

        [z_in, z_out, ~] = sc_map_metrics_z(TA, TBsh, ks, fks, MinBins);
        if isempty(z_in) || isempty(z_out), continue; end

        deltas(rr) = mean(z_in,'omitnan') - mean(z_out,'omitnan');
    end

    permDeltaZ(pidx)=mean(deltas,'omitnan');
end

end

function [z_in, z_out, z_inout] = sc_map_metrics_z(TA, TB, ks, fks, MinBins)
% Compute Option-B SC metrics using a 1-to-1 mapping fks.
% Outputs are pooled across all (k,c) valid correlations.

z_in    = [];
z_out   = [];
z_inout = [];

for i=1:numel(ks)
    k  = ks(i);
    k2 = fks(i);

    A1 = TA{k};  B1 = TB{k};
    A2 = TA{k2}; B2 = TB{k2};

    if isempty(A1) || isempty(B1) || isempty(A2) || isempty(B2)
        continue;
    end

    for c=1:size(A1,1)
        r1 = corr_temporal_vec(A1(c,:), B1(c,:), MinBins); % ININ
        r2 = corr_temporal_vec(A2(c,:), B1(c,:), MinBins); % OUTOUT
        r3 = corr_temporal_vec(A1(c,:), B2(c,:), MinBins); % INOUT

        if isfinite(r1) && isfinite(r2) && isfinite(r3)
            z_in(end+1,1)    = atanh(clamp_r(r1)); %#ok<AGROW>
            z_out(end+1,1)   = atanh(clamp_r(r2)); %#ok<AGROW>
            z_inout(end+1,1) = atanh(clamp_r(r3)); %#ok<AGROW>
        end
    end
end

end

function fks = derangement_map(ks)
% Return a permutation of ks with no fixed points (fks(i) ~= ks(i)).
ks = ks(:);
n = numel(ks);

if n < 2
    fks = ks;
    return;
elseif n == 2
    fks = flipud(ks);
    return;
end

% rejection sampling (fine for small n)
for tries = 1:2000
    p = ks(randperm(n));
    if all(p ~= ks)
        fks = p;
        return;
    end
end

% deterministic fallback: cyclic shift until no fixed points
for shift = 1:n-1
    p = ks([shift+1:n, 1:shift]);
    if all(p ~= ks)
        fks = p;
        return;
    end
end

% last resort
fks = ks(randperm(n));
end


function z = sc_inin_z(TA,TB,ks,MinBins)
% pooled across (k,c): corr across bins for each cell
z = [];
for i=1:numel(ks)
    k = ks(i);
    A = TA{k}; B = TB{k};
    if isempty(A)||isempty(B), continue; end
    for c=1:size(A,1)
        r = corr_temporal_vec(A(c,:),B(c,:),MinBins);
        if isfinite(r)
            z(end+1,1) = atanh(clamp_r(r)); %#ok<AGROW>
        end
    end
end
end


function z = sc_crossk_z(TA,TB,ks,nTarget,MinBins,o,mode)
% sample cross-k correlations, each sample picks random (k1,k2,c)
MaxPairs = o.MaxPairsSC;
nK = numel(ks);
if nK<2, z=[]; return; end

% infer nC
tmp=[];
for i=1:numel(ks)
    if ~isempty(TA{ks(i)}), tmp=TA{ks(i)}; break; end
end
if isempty(tmp), z=[]; return; end
nC = size(tmp,1);

nDraw = min(MaxPairs, max(10*nTarget, 2000));
z = nan(min(nTarget,nDraw),1);
ii=0;

for t=1:nDraw
    k1 = ks(randi(nK));
    k2 = ks(randi(nK));
    if k2==k1, continue; end
    A = TA{k1}; B = TB{k2};
    if isempty(A)||isempty(B), continue; end
    c = randi(nC);
    r = corr_temporal_vec(A(c,:),B(c,:),MinBins);
    if ~isfinite(r), continue; end
    ii=ii+1;
    z(ii)=atanh(clamp_r(r));
    if ii>=nTarget, break; end
end
z = z(1:ii);
end


% ======================================================================
% =================== AGGREGATION + PLOTTING HELPERS ====================
% ======================================================================

function delta = aggregate_delta_across_days(daysStruct)
% daysStruct(di).obs.metricZ is [3x1]; Δ = ININ - OUTOUT.
% Combine perm nulls across days by taking aligned mean across days.

day_delta_obs = [];
day_perm = {};

for di = 1:numel(daysStruct)
    mz = daysStruct(di).obs.metricZ;
    if numel(mz)==3 && all(isfinite(mz))
        day_delta_obs(end+1,1) = mz(1)-mz(2); %#ok<AGROW>
    end

    pv = [];
    if isfield(daysStruct(di),'perm') && isfield(daysStruct(di).perm,'deltaZ')
        pv = daysStruct(di).perm.deltaZ;
        pv = pv(isfinite(pv));
    end
    if ~isempty(pv)
        day_perm{end+1,1} = pv(:); %#ok<AGROW>
    end
end

if isempty(day_delta_obs)
    delta = struct('obs_z',NaN,'obs_r',NaN,'perm_z',[],'p_left',NaN,'p_right',NaN,'p_two',NaN);
    return;
end

obs = mean(day_delta_obs,'omitnan');

perm_anim = [];
pL=NaN; pR=NaN; pT=NaN;

if ~isempty(day_perm)
    Pcommon = min(cellfun(@(v) nnz(isfinite(v)), day_perm));
    if ~isempty(Pcommon) && Pcommon>=5
        perm_anim = nan(Pcommon,1);
        for p=1:Pcommon
            vals = nan(numel(day_perm),1);
            for d=1:numel(day_perm)
                v = day_perm{d};
                if numel(v)>=p && isfinite(v(p)), vals(d)=v(p); end
            end
            perm_anim(p)=mean(vals(isfinite(vals)),'omitnan');
        end
        pR = mean(perm_anim >= obs);
        pL = mean(perm_anim <= obs);
        pT = 2*min(pR,pL);
    end
end

delta = struct('obs_z',obs,'obs_r',tanh(obs),'perm_z',perm_anim,'p_left',pL,'p_right',pR,'p_two',pT);
end


function plot_spaceTask_splitHalf_bars_group(ATTS, modeStr)
% 3-bar plot: ININ / OUTOUT / INOUT
% Thin lines: days; Thick: rat means; Bars: mean across rats.
% Prints paired t-tests:
%   - across rats (n=5) using per-rat means
%   - across days (n≈15) pooling all rat-days

Z_day = [];
aid = [];

nR = numel(ATTS.rats);
for ri=1:nR
    D = ATTS.rats(ri).days;
    for di=1:numel(D)
        mz = D(di).obs.metricZ;
        if numel(mz)==3 && all(isfinite(mz))
            Z_day(end+1,:) = mz(:)'; %#ok<AGROW>
            aid(end+1,1) = ri; %#ok<AGROW>
        end
    end
end
if isempty(Z_day)
    warning('No %s data to plot bars.', modeStr); return;
end

% per-rat mean (n=5)
ratZ = nan(nR,3);
for ri=1:nR
    m = (aid==ri);
    if any(m)
        ratZ(ri,:) = mean(Z_day(m,:),1,'omitnan');
    end
end
keepA = all(isfinite(ratZ),2);
ratZk = ratZ(keepA,:);
nA = size(ratZk,1);

% per-day pooled (n≈15)
dayZk = Z_day(all(isfinite(Z_day),2),:);
nD = size(dayZk,1);

% stats (rats)
[p12_r,t12_r,df12_r,dz12_r] = paired_t_z(ratZk(:,1), ratZk(:,2));
[p13_r,t13_r,df13_r,dz13_r] = paired_t_z(ratZk(:,1), ratZk(:,3));

% stats (days)
[p12_d,t12_d,df12_d,dz12_d] = paired_t_z(dayZk(:,1), dayZk(:,2));
[p13_d,t13_d,df13_d,dz13_d] = paired_t_z(dayZk(:,1), dayZk(:,3));

R_day = tanh(Z_day);
ratRk = tanh(ratZk);

figure('Color','w','Position',[140 140 980 560]); hold on
cmap = lines(nR);
make_light = @(c,frac) (1-frac)*c + frac*[1 1 1];

% day thin lines
for ri0=1:nR
    c_base=cmap(ri0,:);
    c_light=make_light(c_base,0.65);
    idx=find(aid==ri0);
    for j=1:numel(idx)
        rr = R_day(idx(j),1:2);
        plot(1:2, rr, '-', 'Color', c_light, 'LineWidth', 1.0);
        plot(1:2, rr, 'o', 'MarkerFaceColor', c_light, 'MarkerEdgeColor', c_light, 'MarkerSize', 4);
    end
end

% rat thick lines
kk=find(keepA);
for ii=1:nA
    ri0=kk(ii);
    rr=ratRk(ii,1:2);
    plot(1:2, rr, '-', 'Color', cmap(ri0,:), 'LineWidth', 2.7);
    plot(1:2, rr, 'o', 'MarkerFaceColor', cmap(ri0,:), 'MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerSize', 6);
end

bar(1, tanh(mean(ratZk(:,1),'omitnan')), 0.6, 'FaceColor',[0.30 0.60 1.00], 'EdgeColor','k');
bar(2, tanh(mean(ratZk(:,2),'omitnan')), 0.6, 'FaceColor',[0.60 0.60 0.60], 'EdgeColor','k');
%bar(3, tanh(mean(ratZk(:,3),'omitnan')), 0.6, 'FaceColor',[0.85 0.40 0.20], 'EdgeColor','k');

xlim([0.5 3.5]); xticks(1:3);
xticklabels({'IN–IN','OUT–OUT','IN–OUT'});
ylabel('Correlation (r)');
title(sprintf('%s split-half: IN/OUT identity', modeStr));
yline(0,'k:'); grid on; box on

yl=ylim;
dy = 0.06*range(yl);

% annotation: rats
y0 = yl(2) + dy;
line([1 2],[y0 y0],'Color','k','LineWidth',1.2);
text(1.5,y0+0.02*range(yl), ...
    sprintf('RATS (n=%d) IN–IN vs OUT–OUT: t(%d)=%.2f, p=%.3g, dz=%.2f', nA, df12_r,t12_r,p12_r,dz12_r), ...
    'HorizontalAlignment','center');

y1 = y0 + 1.1*dy;
line([1 3],[y1 y1],'Color','k','LineWidth',1.2);
text(2.0,y1+0.02*range(yl), ...
    sprintf('RATS (n=%d) IN–IN vs IN–OUT:  t(%d)=%.2f, p=%.3g, dz=%.2f', nA, df13_r,t13_r,p13_r,dz13_r), ...
    'HorizontalAlignment','center');

% annotation: days
y2 = y1 + 1.4*dy;
line([1 2],[y2 y2],'Color','k','LineWidth',1.2);
text(1.5,y2+0.02*range(yl), ...
    sprintf('DAYS (n=%d) IN–IN vs OUT–OUT: t(%d)=%.2f, p=%.3g, dz=%.2f', nD, df12_d,t12_d,p12_d,dz12_d), ...
    'HorizontalAlignment','center');

y3 = y2 + 1.1*dy;
line([1 3],[y3 y3],'Color','k','LineWidth',1.2);
text(2.0,y3+0.02*range(yl), ...
    sprintf('DAYS (n=%d) IN–IN vs IN–OUT:  t(%d)=%.2f, p=%.3g, dz=%.2f', nD, df13_d,t13_d,p13_d,dz13_d), ...
    'HorizontalAlignment','center');

ylim([yl(1) y3 + 0.12*range(yl)]);

fprintf('\n== %s split-half bars ==\n', modeStr);

fprintf('RATS (n=%d)\n', nA);
fprintf('  IN–IN vs OUT–OUT: t(%d)=%.2f, p=%.3g, dz=%.2f\n', df12_r,t12_r,p12_r,dz12_r);
fprintf('  IN–IN vs IN–OUT:  t(%d)=%.2f, p=%.3g, dz=%.2f\n', df13_r,t13_r,p13_r,dz13_r);

fprintf('DAYS (n=%d)\n', nD);
fprintf('  IN–IN vs OUT–OUT: t(%d)=%.2f, p=%.3g, dz=%.2f\n', df12_d,t12_d,p12_d,dz12_d);
fprintf('  IN–IN vs IN–OUT:  t(%d)=%.2f, p=%.3g, dz=%.2f\n', df13_d,t13_d,p13_d,dz13_d);

end


function plot_spaceTask_splitHalf_permHist_perRat(ATTS, modeStr)
% Per-rat histogram of animal-level perm null for Δz = ININ - OUTOUT.
% Fallback: if r.delta.perm_z is empty but per-day perms exist, rebuild perm_z.

nR=numel(ATTS.rats);
figure('Color','w','Position',[180 180 1100 700]);

for ri=1:nR
    r=ATTS.rats(ri);
    subplot(3,2,ri); hold on

    % pull obs
    obs = NaN;
    if isfield(r,'delta') && ~isempty(r.delta) && isfield(r.delta,'obs_z')
        obs = r.delta.obs_z;
    end

    % pull perm (preferred)
    perm = [];
    if isfield(r,'delta') && ~isempty(r.delta) && isfield(r.delta,'perm_z')
        perm = r.delta.perm_z;
    end
    perm = perm(isfinite(perm));

    % fallback: rebuild from day-level perms if needed
    if isempty(perm) && isfield(r,'days') && ~isempty(r.days)
        day_perm = {};
        for di=1:numel(r.days)
            if isfield(r.days(di),'perm') && isfield(r.days(di).perm,'deltaZ')
                v = r.days(di).perm.deltaZ;
                v = v(isfinite(v));
                if ~isempty(v), day_perm{end+1,1} = v(:); end %#ok<AGROW>
            end
        end
        if ~isempty(day_perm)
            Pcommon = min(cellfun(@(v) nnz(isfinite(v)), day_perm));
            if ~isempty(Pcommon) && Pcommon>=5
                perm = nan(Pcommon,1);
                for p=1:Pcommon
                    vals = nan(numel(day_perm),1);
                    for d=1:numel(day_perm)
                        v = day_perm{d};
                        if numel(v)>=p && isfinite(v(p)), vals(d)=v(p); end
                    end
                    perm(p)=mean(vals(isfinite(vals)),'omitnan');
                end
            end
        end
    end

    if isempty(perm) || ~isfinite(obs)
        title(sprintf('%s (no perm)', r.name));
        axis off;
        continue
    end

    histogram(perm,30,'Normalization','pdf','EdgeColor','none');
    yl=ylim;
    plot([obs obs],yl,'k-','LineWidth',2);

    pR=mean(perm>=obs);
    pL=mean(perm<=obs);
    pT=2*min(pR,pL);

    xlabel('\Delta z (IN–IN - OUT–OUT)');
    ylabel('Null density');
    title(sprintf('%s: %s perm', r.name, modeStr));
    grid on; box on
    text(obs,yl(2),sprintf('  obs=%.3f (r≈%.3f)\np_{two}=%.3g',obs,tanh(obs),pT), ...
        'VerticalAlignment','top','HorizontalAlignment','left');
end

sgtitle(sprintf('%s split-half shuffled-label null (animal-level Δz)', modeStr));
end


function rr = clamp_r(r)
rr = max(min(r,0.999999),-0.999999);
end

function [hA, hB] = split_half(trs)
trs = trs(:);
n = numel(trs);
if n < 2
    hA = []; hB = []; return;
end
perm = trs(randperm(n));
nA = floor(n/2);
hA = perm(1:nA);
hB = perm(nA+1:end);
if isempty(hA) || isempty(hB)
    hA = []; hB = [];
end
end

function v = mean_cat_PV(PV, trs, b, k)
% mean of PV{tr,b,k} across trials trs -> vector [nCells x 1]
vv = [];
for i = 1:numel(trs)
    tr = trs(i);
    x = PV{tr,b,k};
    if isempty(x), continue; end
    if isempty(vv)
        vv = double(x(:));
    else
        vv = vv + double(x(:));
    end
end
if isempty(vv)
    v = [];
else
    v = vv / max(1, numel(trs));
end
end

function M = mean_cat_V(V, trs, k, nC, B)
% mean of V{tr,k} across trials -> [nC x B]
acc = zeros(nC,B);
cnt = zeros(nC,B);

for i = 1:numel(trs)
    tr = trs(i);
    X = V{tr,k};
    if isempty(X), continue; end
    X = double(X);
    m = isfinite(X);
    acc(m) = acc(m) + X(m);
    cnt(m) = cnt(m) + 1;
end

if all(cnt(:)==0)
    M = [];
    return;
end

M = nan(nC,B);
m = cnt>0;
M(m) = acc(m) ./ cnt(m);
end

function [obsZ, permDeltaZ] = pv_metrics_fast(PV, PVmask, o)
% OPTION B IMPLEMENTATION (PV):
% For each (rr,b), build split-half templates A{k}, B{k} for each valid k.
% Then create a 1-to-1 "derangement" mapping f over ks (no k maps to itself).
% Metrics per k:
%   ININ   = corr(A{k},   B{k})
%   OUTOUT = corr(A{f(k)},B{k})
%   INOUT  = corr(A{k},   B{f(k)})
%
% obsZ: [3x1] mean Fisher-z across bins and resamples for [ININ OUTOUT INOUT]
% permDeltaZ: [NPermx1] null for Δ = mean(ININ)-mean(OUTOUT) using label-shuffle of B across ks.

[~,B,K] = size(PVmask);
NResamp  = o.NResamp;
NPerm    = o.NPerm;
MinT     = o.MinTrialsPerK;

% Precompute available trials per (b,k)
trsByBK = cell(B,K);
for b = 1:B
    for k = 1:K
        trsByBK{b,k} = find(PVmask(:,b,k));
    end
end

% Cache templates per (rr,b): templA/templB and ks list
CACHE = repmat(struct('ks',[],'A',[],'B',[],'ok',false), NResamp, B);

% Build cache (observed templates)
for rr = 1:NResamp
    for b = 1:B
        templA = cell(K,1);
        templB = cell(K,1);
        okK    = false(K,1);

        for k = 1:K
            trs = trsByBK{b,k};
            if numel(trs) < MinT, continue; end
            [hA,hB] = split_half(trs);
            if isempty(hA) || isempty(hB), continue; end

            vA = mean_cat_PV(PV, hA, b, k);
            vB = mean_cat_PV(PV, hB, b, k);
            if isempty(vA) || isempty(vB), continue; end

            templA{k} = vA;
            templB{k} = vB;
            okK(k)    = true;
        end

        ks = find(okK);
        if numel(ks) < 2
            CACHE(rr,b).ok = false;
        else
            CACHE(rr,b).ok = true;
            CACHE(rr,b).ks = ks;
            CACHE(rr,b).A  = templA;
            CACHE(rr,b).B  = templB;
        end
    end
end

% ---------- observed metrics ----------
z_res = nan(NResamp,3,B);

for rr = 1:NResamp
    for b = 1:B
        if ~CACHE(rr,b).ok, continue; end
        ks = CACHE(rr,b).ks;
        A  = CACHE(rr,b).A;
        Bc = CACHE(rr,b).B;

        % Option B: one-to-one derangement map
        fks = derangement_map(ks);

        [z_in, z_out, z_inout] = pv_map_metrics_z(A, Bc, ks, fks);

        if isempty(z_in), continue; end
        z_res(rr,1,b) = mean(z_in,    'omitnan');
        z_res(rr,2,b) = mean(z_out,   'omitnan');
        z_res(rr,3,b) = mean(z_inout, 'omitnan');
    end
end

m = squeeze(mean(z_res,3,'omitnan'));   % [NResamp x 3]
obsZ = mean(m,1,'omitnan')';
obsZ = obsZ(:);

% ---------- permutation null for Δ ----------
permDeltaZ = nan(NPerm,1);
if NPerm < 5 || all(~isfinite(obsZ)), return; end

for pidx = 1:NPerm
    deltas_rrb = nan(NResamp,B);

    for rr = 1:NResamp
        for b = 1:B
            if ~CACHE(rr,b).ok, continue; end
            ks = CACHE(rr,b).ks;
            A  = CACHE(rr,b).A;
            Bc = CACHE(rr,b).B;

            % shuffle labels of B across ks (same as your original null)
            ks_shuf = ks(randperm(numel(ks)));
            Bsh = cell(K,1);
            for ii = 1:numel(ks)
                Bsh{ks(ii)} = Bc{ks_shuf(ii)};
            end

            fks = derangement_map(ks);
            [z_in, z_out, ~] = pv_map_metrics_z(A, Bsh, ks, fks);
            if isempty(z_in) || isempty(z_out), continue; end

            deltas_rrb(rr,b) = mean(z_in,'omitnan') - mean(z_out,'omitnan');
        end
    end

    permDeltaZ(pidx) = mean(deltas_rrb(:),'omitnan');
end

end

function [z_in, z_out, z_inout] = pv_map_metrics_z(A, B, ks, fks)
% Compute Option-B PV metrics for one (rr,b) using a 1-to-1 mapping fks.
% ks and fks are same length; fks(i) != ks(i) for all i.

n = numel(ks);
z_in    = nan(n,1);
z_out   = nan(n,1);
z_inout = nan(n,1);

ii = 0;
for i = 1:n
    k  = ks(i);
    k2 = fks(i);

    if isempty(A{k}) || isempty(B{k}) || isempty(A{k2}) || isempty(B{k2})
        continue;
    end

    [r1,~] = safe_corr(A{k},  B{k},  3, 1e-12);  % ININ
    [r2,~] = safe_corr(A{k2}, B{k},  3, 1e-12);  % OUTOUT
    [r3,~] = safe_corr(A{k},  B{k2}, 3, 1e-12);  % INOUT

    if ~isfinite(r1) || ~isfinite(r2) || ~isfinite(r3), continue; end

    ii = ii + 1;
    z_in(ii)    = atanh(clamp_r(r1));
    z_out(ii)   = atanh(clamp_r(r2));
    z_inout(ii) = atanh(clamp_r(r3));
end

z_in    = z_in(1:ii);
z_out   = z_out(1:ii);
z_inout = z_inout(1:ii);
end
