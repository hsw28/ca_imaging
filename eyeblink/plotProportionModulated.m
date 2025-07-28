function change = plotProportionModulated()
% plotProportionModulated  Fraction of neurons with ↑FR in trace (vs null)
ratNames   = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
win        = [0 2];
minSpikes  = 0;
nPerm      = 500;
alpha      = 0.05;

nRats = numel(ratNames);
all_obsFold = [];
all_h       = [];
fracRat     = nan(nRats,1);
stdRat      = nan(nRats,1);

figure('Color','w','Position',[100 100 1200 600]);
change = NaN(nRats,3000);

for r = 1:nRats+1
  if r <= nRats
    rat   = evalin('base', ratNames{r});
    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An),1);
    days  = dates(idx-2:idx);

    cells = struct('st',{},'cs',{});
    for d = 1:3
      spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
      csTimes  = rat.CS_times.(sprintf('CS_%s',days{d}));
      ratemask = rat.ratemask.(sprintf('ratemask_%s',days{d}));

      [nCells,~] = size(spk);
      for c = 1:nCells
        if ratemask(c) == 0, continue, end                 % ← exclude cell
        times = spk(c,:);
        times = times(~isnan(times)&times>0);
        cells(end+1).st = times(:);
        cells(end).cs   = csTimes;
      end
    end
    subplot(2,3,r);
    titleTxt = ratNames{r};
  else
    cells = struct('st',{},'cs',{});
    for rr = 1:nRats
      rat   = evalin('base', ratNames{rr});
      dates = autoDateList(rat);
      idx   = find(strcmp(dates, rat.An),1);
      days  = dates(idx-2:idx);
      for d = 1:3
        spk      = rat.Ca_peaks.(sprintf('CA_peaks_%s',days{d}));
        csTimes  = rat.CS_times.(sprintf('CS_%s',days{d}));
        ratemask = rat.ratemask.(sprintf('ratemask_%s',days{d}));

        [nCells,~] = size(spk);
        for c = 1:nCells
          if ratemask(c) == 0, continue, end               % ← exclude cell
          times = spk(c,:);
          times = times(~isnan(times)&times>0);
          cells(end+1).st = times(:);
          cells(end).cs   = csTimes;
        end
      end
    end
    subplot(2,3,6);
    titleTxt = 'All rats';
  end

  % ------------- downstream analysis unchanged --------------------------
  nCells    = numel(cells);
  obsFold   = nan(nCells,1);
  h         = false(nCells,1);
  for i = 1:nCells
    st = cells(i).st;
    cs = cells(i).cs;
    validCS = cs(cs+win(2)<=max(st));
    nT = numel(validCS);
    if nT<1, continue; end

    FRt = nan(nT,1);
    for t = 1:nT
      t0 = validCS(t)+win(1);
      t1 = t0+diff(win);
      FRt(t) = sum(st>=t0 & st<t1)/diff(win);
    end

    maskCS = false(size(st));
    for t = 1:nT
      maskCS = maskCS | (st>=validCS(t)+win(1)&st<validCS(t)+win(2));
    end
    totalNon = (max(st)-min(st)) - nT*diff(win);
    FRr = (numel(st)-sum(maskCS)) / totalNon;

    obsFold(i) = mean(FRt)/FRr;
    if sum(st>=validCS(1)&st<validCS(end)+win(2)) < minSpikes
      if r<=nRats, change(r,i) = NaN; end
      continue
    elseif r<=nRats
      change(r,i) = obsFold(i);
    end

    nullF = nan(nPerm,1);
    tStarts = linspace(min(st),max(st)-diff(win),1000);
    for ip = 1:nPerm
      samp = randsample(tStarts,nT,true);
      FRs  = arrayfun(@(s) sum(st>=s&st<s+diff(win))/diff(win), samp);
      nullF(ip) = mean(FRs)/FRr;
    end
    pval = mean(nullF >= obsFold(i));
    h(i) = (pval<alpha);
  end

  all_obsFold = [all_obsFold; obsFold];
  all_h       = [all_h;       h];

  if r<=nRats, fracRat(r) = mean(h,'omitnan'); end

  maxFold = nanmax(obsFold);
  edges   = linspace(0, ceil(maxFold), 30);
  histogram(obsFold(~h), edges, 'Normalization','probability'); hold on;
  histogram(obsFold(h),  edges, 'Normalization','probability');
  xlabel('Fold‐change'); ylabel('Probability');
  title(titleTxt);
  legend('ns','sig','Location','Best');
end

fprintf('\n=== Fraction modulated by rat ===\n');
for r = 1:nRats
  fprintf('%s: %.3f\n', ratNames{r}, fracRat(r));
end
fracAll = mean(fracRat,'omitnan');
fprintf('All rats combined: %.3f\n', fracAll);
end
