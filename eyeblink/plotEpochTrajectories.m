function plotEpochTrajectories(pos, CS_onsets)
% plotEpochTrajectories  Plot real movement trajectories in each trial epoch
%
% pos        [n×3] = [time, x, y]
% CS_onsets  [m×1] = timestamps of CS onset (one per trial)
%
% Epochs (relative to CS):
%   CS:    [0.00 0.25] s → blue
%   Trace: [0.25 0.75] s → green
%   US:    [0.75 0.85] s → red
%   Post:  [0.85 2.00] s → black

% unpack
t_pos = pos(:,1);
x_pos = pos(:,2);
y_pos = pos(:,3);

% define epochs
epochs = { [0    0.25], ...
           [0.25 0.75], ...
           [0.75 0.85], ...
           [0.85 2.00] };
colors = { [0 0 1], [0 0.7 0], [1 0 0], [0 0 0] };
names  = {'CS (0–250 ms)', 'Trace (250–750 ms)', ...
          'US (750–850 ms)', 'Post (850–2000 ms)'};

figure
% for each single-epoch figure
for e = 1:4
    subplot(1,5,e)
    hold on;
    % background full path in light gray
    plot(x_pos, y_pos, '-', 'Color',[0.6 0.6 0.6], 'LineWidth',1);

    for i = 1:numel(CS_onsets)
        t0 = CS_onsets(i);
        % extract that trial’s continuous segment
        t_seg = t_pos >= (t0 + epochs{e}(1)) & t_pos <= (t0 + epochs{e}(2));
        plot(x_pos(t_seg), y_pos(t_seg), '-', ...
             'Color', colors{e}, 'LineWidth',1.5);
    end
    title(['Trajectory during ', names{e}]);
    xlabel('X Position'); ylabel('Y Position');
    axis([55 145, 70 130])
    %axis equal;
    %box off;
end

% composite figure
subplot(1,5,5)
hold on;
plot(x_pos, y_pos, '-', 'Color',[0.6 0.6 0.6], 'LineWidth',1);

for e = 1:4
    for i = 1:numel(CS_onsets)
        t0 = CS_onsets(i);
        t_seg = t_pos >= (t0 + epochs{e}(1)) & t_pos <= (t0 + epochs{e}(2));
        plot(x_pos(t_seg), y_pos(t_seg), '-', ...
             'Color', colors{e}, 'LineWidth',1.5);
    end
end
    axis([55 145, 70 130])
legend(['Full path', names], 'Location','best');
title('All epochs overlaid');
xlabel('X Position'); ylabel('Y Position');
%axis equal;
box off;

% … after all subplots have been drawn …

% after drawing your figure…
set(gcf,'Renderer','painters');
print(gcf, 'fourperiods.svg', '-dsvg', '-r600');



end
