function f = plotaid(scores, CSUS_vector)
%function f = plotaid(scores1, scores2, scores3, scores1B, scores2B, scores3B, overall_scores)


  %take three outputs from plotPCAscores w
  %makes two graphs-- one with all trials, one with averages

  %can also do it with isomap output as such: plotaid(Y.coords{6}', training21(wanted21))

set(0,'DefaultFigureVisible', 'off');
[scores1 scores1B] = plotPCAscores(scores, CSUS_vector, 1);
[scores2 scores2B] = plotPCAscores(scores, CSUS_vector, 2);
[scores3 scores3B CS_length] = plotPCAscores(scores, CSUS_vector, 3);

set(0,'DefaultFigureVisible', 'on');
figure


subplot(2,1,1)
%plot3(overall_scores(1:3246,1), overall_scores(1:3246,2), overall_scores(1:3246,3), 'color', [.5 .5 .5])
%hold on
%plot3(overall_scores(3247:end,1), overall_scores(3247:end,2), overall_scores(3247:end,3), 'color', 'black')

%plot3(scores1(1:CS_length,1), scores2(1:CS_length,1), scores3(1:CS_length,1), 'color', 'magenta', 'LineWidth', .5)
%plot3(scores1(CS_length:end,1), scores2(CS_length:end,1), scores3(CS_length:end,1), 'color', 'cyan', 'LineWidth', .5)
plot3(scores1B(1:CS_length,1), scores2B(1:CS_length,1), scores3B(1:CS_length,1), 'color', 'red', 'LineWidth', .5)
plot3(scores1B(CS_length:end,1), scores2B(CS_length:end,1), scores3B(CS_length:end,1), 'color', 'green', 'LineWidth', .5)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('All trials, ENV A CS = pink, A-US = blue, B-CS = red, B-US = green')

hold on
for k=2:size(scores1B,2)
%  plot3(scores1(1:CS_length,k), scores2(1:CS_length,k), scores3(1:CS_length,k), 'color', 'magenta', 'LineWidth', .5)
%  plot3(scores1(CS_length:end,k), scores2(CS_length:end,k), scores3(CS_length:end,k), 'color', 'cyan', 'LineWidth', .5)
  plot3(scores1B(1:CS_length,k), scores2B(1:CS_length,k), scores3B(1:CS_length,k), 'color', 'red', 'LineWidth', .5)
  plot3(scores1B(CS_length:end,k), scores2B(CS_length:end,k), scores3B(CS_length:end,k), 'color', 'green', 'LineWidth', .5)
end

subplot(2,1,2)
%plot3(nanmean(scores1(1:CS_length,:)'), nanmean(scores2(1:CS_length,:)'), nanmean(scores3(1:CS_length,:)'), 'color', 'magenta', 'LineWidth', 2)
%hold on
plot3(nanmean(scores1B(1:CS_length,:)'), nanmean(scores2B(1:CS_length,:)'), nanmean(scores3B(1:CS_length,:)'), 'color', 'red', 'LineWidth', 2)
hold on

%plot3(nanmean(scores1(CS_length:end,:)'), nanmean(scores2(CS_length:end,:)'), nanmean(scores3(CS_length:end,:)'), 'color', 'cyan', 'LineWidth', 2)
plot3(nanmean(scores1B(CS_length:end,:)'), nanmean(scores2B(CS_length:end,:)'), nanmean(scores3B(CS_length:end,:)'), 'color', 'green', 'LineWidth', 2)
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
xlabel('Dimension 1')
ylabel('Dimension 2')
zlabel('Dimension 3')
title('Average of trials, ENV A CS = pink, A-US = blue, B-CS = red, B-US = green')
