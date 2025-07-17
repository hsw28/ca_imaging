function plotCRtoUSmovement(pos, CR_onsets)
% pos: [n x 3] matrix of [time, x, y]
% CR_onsets: vector of timestamps (CR starts)

% Extract position vectors
t_pos = pos(:,1);
x = pos(:,2);
y = pos(:,3);

a = 0
b = .75

% Interpolate positions at CR and US
x_cr = interp1(t_pos, x, CR_onsets + a, 'linear', 'extrap');
y_cr = interp1(t_pos, y, CR_onsets + a, 'linear', 'extrap');
x_us = interp1(t_pos, x, CR_onsets + b, 'linear', 'extrap');
y_us = interp1(t_pos, y, CR_onsets + b, 'linear', 'extrap');

% Displacement vectors
dx = x_us - x_cr;
dy = y_us - y_cr;

% Plot full path
%figure('Color', 'w');
plot(x, y, '-', 'Color', [0.6 0.6 0.6], 'MarkerSize', 2);  % background path
hold on;

% Plot CR onset (circle) and movement arrow to US
for i = 1:length(CR_onsets)
    plot(x_cr(i), y_cr(i), 'ro', 'Color', [1 0 0], 'MarkerSize', 6, 'LineWidth', 1.5);       % CR start
end

% Draw arrows with quiver
quiver(x_cr, y_cr, dx, dy, 0, ...
       'Color', [1 0 0], ...
       'LineWidth', 2, ...
       'MaxHeadSize', 0.2, ...
       'AutoScale', 'on');


% Labels and aesthetics
xlabel('X Position');
ylabel('Y Position');
title('CR Onset and Movement to US (0.75s Later)');
%title('US Onset and Movement to 2s after CS');
axis equal;
set(gca, 'Box', 'off', 'FontSize', 12);
legend({'Full trajectory', 'CR onset', 'Movement to US'}, 'Location', 'best');

axis([50 150, 70 125])
