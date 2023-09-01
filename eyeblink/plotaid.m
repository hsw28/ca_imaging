function f = plotaid(scores1, scores2, scores3, scores1B, scores2B, scores3B)
%function f = plotaid(scores1, scores2, scores3, scores1B, scores2B, scores3B, overall_scores)

  %take three outputs from plotPCAscores which gets scores for individual trials
  %makes two graphs-- one with all trials, one with averages

subplot(1,2,1)
%plot3(overall_scores(1:3246,1), overall_scores(1:3246,2), overall_scores(1:3246,3), 'color', [.5 .5 .5])
%hold on
%plot3(overall_scores(3247:end,1), overall_scores(3247:end,2), overall_scores(3247:end,3), 'color', 'black')

plot3(scores1(1:6,1), scores2(1:6,1), scores3(1:6,1), 'color', 'magenta', 'LineWidth', .5)
plot3(scores1(6:end,1), scores2(6:end,1), scores3(6:end,1), 'color', 'cyan', 'LineWidth', .5)
plot3(scores1B(1:6,1), scores2B(1:6,1), scores3B(1:6,1), 'color', 'red', 'LineWidth', .5)
plot3(scores1B(6:end,1), scores2B(6:end,1), scores3B(6:end,1), 'color', 'green', 'LineWidth', .5)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('All trials, ENV A CS = pink, A-US = blue, B-CS = red, B-US = green')

hold on
for k=2:size(scores1,2)
  plot3(scores1(1:6,k), scores2(1:6,k), scores3(1:6,k), 'color', 'magenta', 'LineWidth', .5)
  plot3(scores1(6:end,k), scores2(6:end,k), scores3(6:end,k), 'color', 'cyan', 'LineWidth', .5)
  plot3(scores1B(1:6,k), scores2B(1:6,k), scores3B(1:6,k), 'color', 'red', 'LineWidth', .5)
  plot3(scores1B(6:end,k), scores2B(6:end,k), scores3B(6:end,k), 'color', 'green', 'LineWidth', .5)
end

subplot(1,2,2)
plot3(nanmean(scores1(1:6,:)'), nanmean(scores2(1:6,:)'), nanmean(scores3(1:6,:)'), 'color', 'magenta', 'LineWidth', 2)
hold on
plot3(nanmean(scores1B(1:6,:)'), nanmean(scores2B(1:6,:)'), nanmean(scores3B(1:6,:)'), 'color', 'red', 'LineWidth', 2)

plot3(nanmean(scores1(6:end,:)'), nanmean(scores2(6:end,:)'), nanmean(scores3(6:end,:)'), 'color', 'cyan', 'LineWidth', 2)
plot3(nanmean(scores1B(6:end,:)'), nanmean(scores2B(6:end,:)'), nanmean(scores3B(6:end,:)'), 'color', 'green', 'LineWidth', 2)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('Average of trials, ENV A CS = pink, A-US = blue, B-CS = red, B-US = green')
