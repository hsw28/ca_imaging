function [scorematrixA scorematrixB] = plotPCAscores(scores, envCSUSmatrix, princcopnumber)
  %plots principle component scores for all eyeblink trials
  %CS/US matrix should come from movingtimetraining or equivalent. will be a matrix where first row is either 1 or 2 to indicate env A or B, second row is 10 or 20 or 0 to indicate CS/US/neither

k=1;
cstimes = find(envCSUSmatrix(1,:)==10);
ustimes = find(envCSUSmatrix(1,:)==20);
allcs_start = [];
allus_start = [];

while k<=length(envCSUSmatrix)
  overtime = (find(cstimes>k));
  overtime = min(overtime);
    if length(overtime>0);
      overtime = overtime(1);
      overtime = cstimes(overtime);
      allcs_start(end+1) = overtime;
    end
  k=overtime+7;
end



scorematrixA = NaN(9, length(allcs_start));
scorematrixB = NaN(9, length(allcs_start));
env = NaN(1, length(allcs_start));
for i=1:length(allcs_start)
  curstart = allcs_start(i);
  curend = curstart+8;
  wantedscore = scores(curstart:curend, princcopnumber);
  if envCSUSmatrix(1,curstart) == 1
  scorematrixA(:,i) = wantedscore;
  else

  scorematrixB(:,i) = wantedscore;
  end
end

scorematrixA;
scorematrixB;


figure
plot(scorematrixA, 'Color', [0.3010 0.7450 0.9330], 'LineWidth', .2);
hold on
plot(scorematrixB, 'Color', [.9 .7 .7], 'LineWidth', .2);
plot(nanmean(scorematrixA'), 'Color', [0 0.4470 0.7410], 'LineWidth', 4)
plot(nanmean(scorematrixB'), 'Color', 'red', 'LineWidth', 4)
vline(5)
vline(5)
vline(5)
