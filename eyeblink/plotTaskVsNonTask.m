function [rMat,pMat] = plotTaskVsNonTask(ratNames, rateMode, preWindow)
% plotTaskVsNonTask   Compare Pre- vs Task-epoch rates and partial residuals
%   [rMat,pMat] = plotTaskVsNonTask(ratNames)
%   [rMat,pMat] = plotTaskVsNonTask(ratNames, rateMode)
%   [rMat,pMat] = plotTaskVsNonTask(ratNames, rateMode, preWindow)
%
% Columns:
%  1) normPre vs normTask
%  2) ratePre vs rateTask
%  3) residPre_MI vs residTask_MI
%  4) residPre_bpspike vs residTask_bpspike
%  5) residPre_bps vs residTask_bps

if nargin<2, rateMode = 'pre'; end
if nargin<3, preWindow = [-5 0]; else preWindow = preWindow(:)'; end

% set x-axis label based on rateMode
switch rateMode
    case 'pre',      xLab = 'Pre';
    case 'all',      xLab = 'All';
    case 'all_fast', xLab = 'All fast';
    otherwise        xLab = 'Pre';
end

taskWindow = [0 2];
nRats = numel(ratNames);
rMat = nan(nRats+1,5);
pMat = nan(nRats+1,5);

% pooled accumulators
ALL.normPre        = [];
ALL.normTask       = [];
ALL.ratePre        = [];
ALL.rateTask       = [];
ALL.MI_pre         = [];
ALL.MI_task        = [];
ALL.bpspike_pre    = [];
ALL.bpspike_task   = [];
ALL.bps_pre        = [];
ALL.bps_task       = [];

figure('Color','w','Position',[100 100 1800 300*(nRats+1)]);

% per-rat rows
for r=1:nRats
    rat = evalin('base', ratNames{r});
    days = autoDateList(rat);
    iAn  = find(strcmp(days, rat.An),1);
    days = days(iAn-2:iAn);

    % collect spikes & metrics
    spikesAll = {};
    MIpre = []; MIpost = [];
    bpsp_pre = []; bpsp_post = [];
    bps_pre = []; bps_post = [];
    csAll = [];
    allPos = []; allVel = [];
    trialInts = [];
    for d=1:3
        day = days{d};
        S = rat.Ca_peaks.(sprintf('CA_peaks_%s',day));
        M = rat.ratemask.(sprintf('ratemask_%s',day))==1;
        cs = rat.CS_times.(sprintf('CS_%s',day)); cs=cs(:);
        miS = rat.MI_noCSUS.(sprintf('MI_%s',day));
        miT = rat.MI_CSUS15.(sprintf('MI_%s',day));
        Bt_all      = rat.bitsper.(sprintf('MI_%s',day));
        Bt_t_all    = rat.bitsperCSUS.(sprintf('MI_%s',day));
        bpsp_pre_d  = Bt_all(1,:)';    bps_pre_d  = Bt_all(2,:)';
        bpsp_post_d = Bt_t_all(1,:)';  bps_post_d = Bt_t_all(2,:)';
        pos = rat.pos.(sprintf('pos_%s',day)); vel = ca_velocity(pos);
        allPos = [allPos; pos(:,1)]; allVel = [allVel, vel];
        for t=1:numel(cs)
            trialInts(end+1,:) = [cs(t)+preWindow(1), cs(t)+taskWindow(2)];
        end
        csAll = [csAll; cs];
        for c=1:size(S,1)
            if ~M(c), continue; end
            st = S(c,:); st = st(~isnan(st)&st>0);
            if numel(st)<3, continue; end
            spikesAll{end+1,1} = st;
            MIpre(end+1,1)      = miS(c);
            MIpost(end+1,1)     = miT(c);
            bpsp_pre(end+1,1)    = bpsp_pre_d(c);
            bpsp_post(end+1,1)   = bpsp_post_d(c);
            bps_pre(end+1,1)     = bps_pre_d(c);
            bps_post(end+1,1)    = bps_post_d(c);
        end
    end
    N = numel(spikesAll);

    % compute Pre/Task rates
    ratePre  = nan(N,1);
    rateTask = nan(N,1);
    for i=1:N
        st = spikesAll{i};
        cntP = sum(arrayfun(@(t) sum(st>=t+preWindow(1)&st<t+preWindow(2)), csAll));
        cntT = sum(arrayfun(@(t) sum(st>=t+taskWindow(1)&st<t+taskWindow(2)), csAll));
        ratePre(i)  = cntP/(numel(csAll)*diff(preWindow));
        rateTask(i) = cntT/(numel(csAll)*diff(taskWindow));
    end

    % session mean rate
    switch rateMode
        case 'pre'
            sessionDur = numel(csAll)*diff(preWindow);
            meanRate   = cellfun(@numel, spikesAll)/sessionDur;
        case 'all'
            sessionDur = totalNonTrialTime(allPos, trialInts);
            meanRate   = countSpikesOutside(spikesAll, allPos, trialInts)/sessionDur;
        case 'all_fast'
            [sessionDur, fastIdx] = totalNonTrialTime(allPos, trialInts, allVel,4);
            meanRate   = countSpikesAtIdx(spikesAll, fastIdx)/sessionDur;
    end
    normPre  = ratePre./meanRate;
    normTask = rateTask./meanRate;

    % accumulate pooled
    ALL.normPre      = [ALL.normPre; normPre];
    ALL.normTask     = [ALL.normTask; normTask];
    ALL.ratePre      = [ALL.ratePre; ratePre];
    ALL.rateTask     = [ALL.rateTask; rateTask];
    ALL.MI_pre       = [ALL.MI_pre; MIpre];
    ALL.MI_task      = [ALL.MI_task; MIpost];
    ALL.bpspike_pre  = [ALL.bpspike_pre; bpsp_pre];
    ALL.bpspike_task = [ALL.bpspike_task; bpsp_post];
    ALL.bps_pre      = [ALL.bps_pre; bps_pre];
    ALL.bps_task     = [ALL.bps_task; bps_post];

    % residualize
    [rmiP,rmiT]   = residualsAgainst(MIpre,   MIpost,   meanRate);
    [rbpP,rbpT]   = residualsAgainst(bpsp_pre,bpsp_post,meanRate);
    [rbpsP,rbpsT] = residualsAgainst(bps_pre, bps_post, meanRate);

    % plot this rat row
    for col=1:5
        ax = subplot(nRats+1,5,(r-1)*5+col); hold(ax,'on');
        switch col
            case 1, x=normPre;    y=normTask;    ttl='Normalized Rate';
            case 2, x=ratePre;    y=rateTask;    ttl='Raw Rate (Hz)';
            case 3, x=rmiP;       y=rmiT;        ttl='MI Residual';
            case 4, x=rbpP;       y=rbpT;        ttl='bps/spike Residual';
            case 5, x=rbpsP;      y=rbpsT;       ttl='bps Residual';
        end
        scatter(ax,x,y,15,'filled'); plot(ax,[min(x) max(x)],[min(y) max(y)],'k--'); axis tight;
        mask = isfinite(x)&isfinite(y);
        if col<3
            [rv,pv] = corr(x(mask), y(mask), 'Rows','complete');
        else
            [rv,pv] = partialcorr(x(mask), y(mask), meanRate(mask),'Rows','complete');
        end
        % only label first subplot
        if r==1 && col==1
            xlabel(ax,xLab);
            ylabel(ax,'Task');
        end
        % only title first row
        if r==1
            title(ax,ttl);
        end
        text(ax,0.05*max(x),0.9*max(y),sprintf('r=%.2f, p=%.3f',rv,pv),'FontSize',10);
        rMat(r,col)=rv; pMat(r,col)=pv;
        hold(ax,'off');
    end
end

% pooled row
row = nRats+1;
maskMI   = isfinite(ALL.MI_pre)&isfinite(ALL.MI_task)&isfinite(ALL.ratePre);
[rpmiP,rpmiT] = residualsAgainst(ALL.MI_pre(maskMI),ALL.MI_task(maskMI),ALL.ratePre(maskMI));
maskBpsp = isfinite(ALL.bpspike_pre)&isfinite(ALL.bpspike_task)&isfinite(ALL.ratePre);
[rpbpP,rpbpT] = residualsAgainst(ALL.bpspike_pre(maskBpsp),ALL.bpspike_task(maskBpsp),ALL.ratePre(maskBpsp));
maskBps  = isfinite(ALL.bps_pre)&isfinite(ALL.bps_task)&isfinite(ALL.ratePre);
[rpbpsP,rpbpsT] = residualsAgainst(ALL.bps_pre(maskBps),ALL.bps_task(maskBps),ALL.ratePre(maskBps));
for col=1:5
    ax = subplot(nRats+1,5,(row-1)*5+col); hold(ax,'on');
    switch col
        case 1, x=ALL.normPre;      y=ALL.normTask;    ttl='All Norm';
        case 2, x=ALL.ratePre;      y=ALL.rateTask;    ttl='All Rate';
        case 3, x=rpmiP;           y=rpmiT;           ttl='All MI Residual';
        case 4, x=rpbpP;           y=rpbpT;           ttl='All bps/spike Residual';
        case 5, x=rpbpsP;          y=rpbpsT;          ttl='All bps Residual';
    end
    scatter(ax,x,y,15,'r','filled'); plot(ax,[min(x) max(x)],[min(y) max(y)],'k--'); axis tight;
    mask = isfinite(x)&isfinite(y);
    if col<3
        [rv,pv] = corr(x(mask), y(mask), 'Rows','complete');
    else
        [rv,pv] = partialcorr(x(mask), y(mask), ALL.ratePre(mask),'Rows','complete');
    end
    text(ax,0.05*max(x),0.9*max(y),sprintf('r=%.2f, p=%.3f',rv,pv),'FontSize',10);
    xlabel(ax,xLab); ylabel(ax,'Task'); title(ax,ttl);
    rMat(row,col)=rv; pMat(row,col)=pv;
    hold(ax,'off');
end
end

% helper functions
function [rp,rt] = residualsAgainst(Apre,Apost,rate)
    mask = isfinite(Apre)&isfinite(Apost)&isfinite(rate);
    X = [ones(sum(mask),1), rate(mask)];
    bp = X \ Apre(mask);
    bt = X \ Apost(mask);
    rp = Apre(mask)  - X*bp;
    rt = Apost(mask) - X*bt;
end

function counts=countSpikesOutside(spikesAll,allPos,trials)
    N=numel(spikesAll); counts=nan(N,1);
    for i=1:N, st=spikesAll{i}; keep=true(size(st));
        for j=1:size(trials,1), keep=keep & ~(st>=trials(j,1)&st<trials(j,2)); end
        counts(i)=sum(keep);
    end
end

function counts=countSpikesAtIdx(spikesAll,idx)
    N=numel(spikesAll); counts=nan(N,1);
    for i=1:N, st=spikesAll{i}; [~,ix]=min(abs(st(:)-idx'),[],2); counts(i)=numel(ix); end
end

function [T,fastIdx]=totalNonTrialTime(allTimes,trials,velData,speedThresh)
    if nargin<4, good=true(size(allTimes)); fastIdx=[];
    else     good=velData(1,:)>=speedThresh; fastIdx=find(good); end
    times=allTimes(good); dt=median(diff(times));
    for i=1:size(trials,1), times(times>=trials(i,1)&times<trials(i,2))=[]; end
    T = numel(times)*dt;
end
