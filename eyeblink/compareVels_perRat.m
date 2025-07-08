function [meanIn,meanOut,pval,velIn,velOut] = compareVels_perRat(ratName,winSecs)
% compareVels_perRat  Compare velocity inside vs. outside CS-aligned windows
%
%   [meanIn,meanOut,pval,velIn,velOut] = velInVsOutTrials(ratName,[0 2])
%
%   * ratName   – variable name of animal in base workspace (char or string)
%   * winSecs   – 1×2 vector [t0 t1] relative to CS (default [0 2]  s)
%
%   Returns:
%     meanIn   – mean velocity  inside trial windows  (cm s⁻¹)
%     meanOut  – mean velocity outside trial windows  (cm s⁻¹)
%     pval     – two-sample t-test p (velIn vs velOut)
%     velIn    – column vector of all in-trial velocity samples
%     velOut   – column vector of all out-of-trial velocity samples
%
%   A box-plot is produced comparing the two distributions.  If the rat has
%   <3 “learned” sessions (final-2:final), the function returns empty arrays.

% -------------------------------------------------------------------------
if nargin<2,  winSecs = [0 2];  end
ratName = char(ratName);

animal = evalin('base',ratName);
pos    = @(day) animal.pos.(['pos_' day]);          % convenience handle
csT    = @(day) animal.CS_times.(['CS_' day]);      %     "        "

dateList = autoDateList(animal);
finIdx   = find(strcmp(dateList,animal.An));

if finIdx<3
    warning('%s has < 3 learned sessions',ratName)
    [meanIn,meanOut,pval,velIn,velOut] = deal([]);
    return
end

learnedDays = dateList(finIdx-2:finIdx);

velIn  = [];   % velocity samples in windows
velOut = [];   % velocity samples outside all windows

for d = 1:numel(learnedDays)
    day = learnedDays{d};
    v           = ca_velocity(pos(day));   % v(1,:) speed, v(2,:) time stamps
    vTime       = v(2,:);    vVal = v(1,:);
    trialMask   = false(size(vVal));       % mark indices that fall in ANY trial window

    csTimes = csT(day);
    for k = 1:numel(csTimes)
        if strcmp(ratName,'rat0816') && k==1,  continue, end
        cs  = csTimes(k);
        inWin = vTime >= cs+winSecs(1) & vTime <= cs+winSecs(2);
        velIn = [velIn; vVal(inWin)'];            %#ok<AGROW>
        trialMask = trialMask | inWin;
    end

    velOut = [velOut; vVal(~trialMask)'];         %#ok<AGROW>
end

% Remove NaNs (just in case)
velIn  = velIn(~isnan(velIn));
velOut = velOut(~isnan(velOut));

if isempty(velIn) || isempty(velOut)
    warning('%s: no valid velocity samples',ratName)
    [meanIn,meanOut,pval] = deal(NaN);
else
    meanIn  = nanmean(velIn)
    meanOut = nanmean(velOut)
    stdIn = nanmean(velIn)
    stdOut = nanmean(velOut)
    [h,pval,ci,stats]  = ttest2(velIn,velOut);          % two-sample t-test
    pval
    stats
end

% ----------------- plot --------------------------------------------------
figure('Color','w','Position',[300 400 350 450]);
boxplot([velOut; velIn], ...
        [repmat({'Out-of-Trial'},numel(velOut),1);
        repmat({'In-Trial'},numel(velIn),1)], ...
        'Colors',[0.2 0.6 0.8; 0.6 0.2 0.2]);
ylabel('Velocity  (cm/s)');
title(sprintf('%s  velocity  in vs out of trials\np = %.3g',ratName,pval));

end
