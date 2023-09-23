function [scorematrixA scorematrixB CS_length] = plotPCAscores(scores, envCSUSmatrix, princcopnumber)
  %plots principle component scores for all eyeblink trials
  %CS/US matrix should come from movingtimetraining or equivalent. will be a matrix where first row is either 1 or 2 to indicate env A or B, second row is 10 or 20 or 0 to indicate CS/US/neither


cstimes = find(envCSUSmatrix(1,:)==10);
ustimes = find(envCSUSmatrix(1,:)==20);
allcs_start = [];
allcs_end = [];

k=1;
l=0;
while l<=length(envCSUSmatrix)
  undertime = (find(cstimes<l));
  undertime = max(undertime);
    if length(undertime>0)
      undertime = undertime(1);
      undertime = cstimes(undertime);
      allcs_end(end+1) = undertime;
      l= undertime+10;
      k = undertime+1;
    else
      l = l+10;
    end

    overtime = (find(cstimes>=k));
    overtime = min(overtime);
      if length(overtime>0)
        overtime = overtime(1);
        overtime = cstimes(overtime);
        allcs_start(end+1) = overtime;
      end

end


mean(allcs_end - allcs_start)
CS_length = round(mean(allcs_end - allcs_start))


allUS_end = [];
z=1;
q=1;
while z<=length(envCSUSmatrix) && q<3
  undertime = (find(ustimes<z));
  undertime = max(undertime);
    if length(undertime>0);
      undertime = undertime(1);
      undertime = ustimes(undertime);
      allUS_end(end+1) = undertime;
      if undertime+10 < length(envCSUSmatrix)
          z=undertime+10;
      else
        z = length(envCSUSmatrix);
        q = q+1;
      end
    else
      z = z+10;
    end
end


len = allUS_end(1)- allcs_start(1);
scorematrixA = NaN(len+1, length(allcs_start));
scorematrixB = NaN(len+1, length(allcs_start));
env = NaN(1, length(allcs_start));
for i=1:length(allcs_start)
  curstart = allcs_start(i);
  curend = allUS_end(i);
  wantedscore = scores(curstart:curend, princcopnumber);
  if envCSUSmatrix(1,curstart) == 1
  scorematrixA(1:length(wantedscore),i) = wantedscore;
  else
  scorematrixB(1:length(wantedscore),i) = wantedscore;
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
