function f = velhisto(spiketimes, pos)
%makes a histogram of spik times and speed

spikestimes = spiketimes(:);
spikestimes = sort(spikestimes);


vel = ca_velocity(pos);
vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30.0005);
edges = [vel(2,1):2:vel(2,end)];

[N,edges] = histcounts(spikestimes,edges);

%interp1(vel(2,:),vel(1,:),edges)

%subplot(2,1,1)
figure
subplot(2,1,1)
plot(vel(2,:),(vel(1,:)*100./max(vel(1,:))), 'LineWidth', 2, 'Color', 'black')
hold on
histogram(spikestimes, 'BinEdges', edges, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', [0.3010 0.7450 0.9330])
xlabel('Time (Sec)')
ylabel('Speed / Calcium Events')
subplot(2,1,2)
plot(vel(2,:),(vel(1,:)*100./max(vel(1,:))), 'LineWidth', 2, 'Color', 'black')
hold on
histogram(spikestimes, 'BinEdges', edges, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', [0.3010 0.7450 0.9330])
xlabel('Time (Sec)')
ylabel('Speed / Calcium Events')


edges = [vel(2,1):.1:vel(2,end)];
[N,edges] = histcounts(spikestimes,edges);
interpv = interp1(vel(2,:),vel(1,:),edges);


f = Hcorr(interpv(1:end-1), N);
[Y,I] = max(f(1,:));
l0 = find(f(2,:)==0);
%f = f(1,l0);
f(:, I)'
