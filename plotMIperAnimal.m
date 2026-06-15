function plotMIperAnimal(ratNames)
% ratNames: e.g., {'rat0314', 'rat0222'}

colors = [0.1 0.6 0.9]; % match previous plot palette
dayLabels = {'An-2', 'An-1', 'An'};

allMeans = [];

figure('Color','w', 'Position', [300 300 800 400]);

% --- Panel 1: Per animal MI
figure
hold on
%subplot(1,2,1); hold on;
for r = 1:length(ratNames)
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    An_day = rat.An;
    idx = find(strcmp(dateList, An_day));
    theseDays = dateList(idx-2:idx);

    for d = 1:3
  %   miVec = rat.MI_noCSUS15.(['MI_' theseDays{d}]);
  %   miVec = rat.MI_noCSUS30.(['MI_' theseDays{d}]);
  %  miVec = rat.MI_CSUS15.(['MI_' theseDays{d}]);
    %miVec = rat.MI_wCSUS.(['MI_' theseDays{d}]);
    % miVec = rat.bitsperCSUS.(['MI_' theseDays{d}]);
    %          miVec = miVec(1,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SPIKE
  %           miVec = miVec(2,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SEC

      miVec = rat.bitsper.(['MI_' theseDays{d}]);
            miVec = miVec(1,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SPIKE
        %    miVec = miVec(2,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SEC

    %  nanny = isnan(miVec2);
    %  miVec(nanny) = NaN;


        x = d + (r-1)*4;  % space between animals
        scatter(repmat(x, size(miVec)), miVec, 10, [0.1 0.6 0.9], 'filled', 'jitter','on', 'jitterAmount',0.2);
        boxchart(repmat(x, size(miVec)), miVec, 'BoxFaceColor', 'black', 'BoxWidth', 0.5, 'MarkerStyle', 'none');
    end
end
% Build custom labels
xticks = [];
xticklabels = {};

for r = 1:length(ratNames)
    xticks = [xticks, (1:3) + (r-1)*4];

    if r == 1
        xticklabels = [xticklabels, {'day n-2', 'day n-1', 'day n'}];
    else
        xticklabels = [xticklabels, {'', '', ''}];
    end
end

figure
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

% Add rat labels centered across their 3 days
for r = 1:length(ratNames)
    xCenter = mean((1:3) + (r-1)*4);
    text(xCenter, ylim(gca)*[0.95; 0], sprintf('rat%d', r), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

ylabel('Mutual Information (bits)');
title('MI per Animal per Day');
set(gca, 'Box', 'off');


% --- Panel 2: Mean MI per animal (across the 3 days)
figure
hold on
%subplot(1,2,2); hold on;

animalMeans = zeros(length(ratNames), 1);
animalSEMs = zeros(length(ratNames), 1);
animalSTDs = zeros(length(ratNames), 1);
animalMED = zeros(length(ratNames), 1);

for r = 1:length(ratNames)
    rat = evalin('base', ratNames{r});
    dateList = autoDateList(rat);
    An_day = rat.An;
    idx = find(strcmp(dateList, An_day));
    theseDays = dateList(idx-2:idx);

    allMI = [];
    for d = 1:3
        rate = rat.ratemask.(['ratemask_' theseDays{d}]);
           miVec = rat.MI_noCSUS15.(['MI_' theseDays{d}]);
        %  miVec = rat.MI_noCSUS30.(['MI_' theseDays{d}]);
%          miVec = rat.MI_CSUS15.(['MI_' theseDays{d}]);
        %  miVec = rat.MI_wCSUS.(['MI_' theseDays{d}]);

        %  miVec = rat.bitsperCSUS.(['MI_' theseDays{d}]);
          %      miVec = miVec(1,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SPIKE
        %        miVec = miVec(2,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SEC

          miVec = rat.bitsper.(['MI_' theseDays{d}]);
               miVec = miVec(1,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SPIKE
        %       miVec = miVec(2,:); %%% ONLY FOR BITS PER -- THIS IS BITS PER SEC

        miVec(rate==0) = NaN;


        allMI = [allMI; miVec(:)];
    end
    animalMeans(r) = nanmean(allMI);
    animalSEMs(r) = nanstd(allMI) / sqrt(numel(allMI));
    animalSTDs(r) = nanstd(allMI);
    animalMED(r)=nanmedian(allMI);
end

animalMeans = animalMeans
animalSTDs = animalSTDs
animalMED = animalMED


hBar = bar(1:length(ratNames), animalMeans, 'FaceColor', [0.5 0.5 0.5]);
errorbar(1:length(ratNames), animalMeans, animalSEMs, 'k', 'LineStyle', 'none', 'LineWidth', 1);

set(gca, 'XTick', 1:length(ratNames), 'XTickLabel', ratNames, 'XTickLabelRotation', 45);
ylabel('Mean MI (bits)');
title('Per-Animal Mean MI');
set(gca, 'Box', 'off');


end
