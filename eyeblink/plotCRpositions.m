function plotCRpositions(pos, CR_onsets)
% pos: [n x 3] matrix (columns: time, x, y)
% CR_onsets: vector of timestamps for CR onsets

time = pos(:,1);
x = pos(:,2);
y = pos(:,3);

% Interpolate X and Y at CR onset times
x_cr = interp1(time, x, CR_onsets, 'linear', 'extrap');
y_cr = interp1(time, y, CR_onsets, 'linear', 'extrap');

% Plot
%figure('Color', 'w');
plot(x, y, '-', 'Color', [0.7 0.7 0.7], 'MarkerSize', 2); hold on;
plot(x_cr, y_cr, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);

xlabel('X Position');
ylabel('Y Position');
title('Animal Position with CR Onset Locations');
axis equal;
set(gca, 'Box', 'off', 'FontSize', 12);
legend({'Position trace', 'CR onset'}, 'Location', 'best');
axis([50 150, 70 125])
