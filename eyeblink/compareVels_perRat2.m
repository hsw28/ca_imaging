function [meanIn,meanPre,pval,velIn,velPre] = compareVels_perRat2(ratName,winSecs)
% compareVels_perRat  Compare velocity inside vs. 2 s before CS instead of all non trial times
%
%   [meanIn,meanPre,pval,velIn,velPre] = compareVels_perRat(ratName,[0 2])
%
%   * velIn  – cm/s in the [t0→t1] s window after each CS
%   * velPre – cm/s in the [t0-2→t0] s window before each CS
%   * pval   – two‐sample t‐test velIn vs velPre
%
%   If fewer than 3 learned sessions, returns empties.

if nargin<2,  winSecs=[0 2];  end
ratName = char(ratName);

animal = evalin('base',ratName);
pos    = @(day) animal.pos.(['pos_' day]);
csT    = @(day) animal.CS_times.(['CS_' day]);

dateList = autoDateList(animal);
finIdx   = find(strcmp(dateList,animal.An));
if finIdx<3
  warning('%s has < 3 learned sessions',ratName)
  [meanIn,meanPre,pval,velIn,velPre] = deal([]);
  return
end

learnedDays = dateList(finIdx-2:finIdx);
velIn  = [];
velPre = [];

for d = 1:numel(learnedDays)
  day   = learnedDays{d};
  v     = ca_velocity(pos(day));     % [speed; time]
  vTime = v(2,:);  vVal = v(1,:);
  cs    = csT(day);

  for k = 1:numel(cs)
    if strcmp(ratName,'rat0816') && k==1,  continue,  end

    % in‐trial window
    inIdx = vTime >= cs(k)+winSecs(1)   & vTime <= cs(k)+winSecs(2);
    velIn  = [velIn;  vVal(inIdx).'];

    % pre‐trial window: exactly 2 s before the CS window
    preStart = cs(k) + winSecs(1) - 2;
    preEnd   = cs(k) + winSecs(1);
    preIdx   = vTime >= preStart & vTime <  preEnd;
    velPre  = [velPre; vVal(preIdx).'];
  end
end

% drop NaNs
velIn  = velIn (~isnan(velIn));
velPre = velPre(~isnan(velPre));

if numel(velIn)<1 || numel(velPre)<1
  warning('%s: no valid velocity samples',ratName)
  [meanIn,meanPre,pval] = deal(NaN);
else
  meanIn  = mean(velIn);
  meanPre = mean(velPre);
  stdIn   = std(velIn);
  stdPre  = std(velPre);

  [~,pval,ci,stats] = ttest2(velIn,velPre);  % if you prefer non‐parametric: signrank
end

% box‐plot
figure('Color','w','Position',[300 400 350 450]);
boxplot([velPre; velIn], ...
        [repmat({'Pre-CS'}, numel(velPre),1);
         repmat({'In-CS'},  numel(velIn),1)], ...
        'Colors',[0.2 0.6 0.8; 0.6 0.2 0.2]);
ylabel('Velocity (cm/s)');
title(sprintf('%s: pre-CS vs in-CS velocity (p = %.3g)',ratName,pval));

end
