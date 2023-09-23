close all;
clear all;

rng(1);

N = 60;   % nuber of place cells
T = 2000;
dim = 2;

% Hyperparameters for the simulated foraging rat
s = 0.1;
r = 1;
noise = 0; %1.5; %1.5; % 0 for closed balls 

%%% First, calculate the trajectory
Tpool = round(10*T);
xpool = [];
while length(xpool) < Tpool
    pos_candidates = 2*(rand(1,dim)-0.5);
    if sum( pos_candidates.^2 ) < 1;       %for dim=3, pi/6 = 52.4% work out. But his is exponentially bad for higher dimensions!
        xpool = [xpool; pos_candidates];
        disp(length(xpool)/Tpool)
    end
end

trajectory = zeros(T, dim);
unusedflag = ones(Tpool,1, 'logical');

trajectory(1,:) = xpool(1,:);
unusedflag(1) = false;

[val, ind] = sort( (  sqrt( sum((trajectory(1,:)-xpool).^2, 2) ) - s*r  ).^2  );
trajectory(2,:) = xpool(ind(1),:);
unusedflag(ind(1)) = false;

for i = 3:T
    x_tminus = trajectory(i-2,:);
    x_t = trajectory(i-1,:);
    r = rand();
    allowed_list = xpool(unusedflag,:);
    [val, ind] = sort( (  sqrt( sum((x_t-allowed_list).^2, 2) ) - s*r  ).^2  );

    testpoint_inacceptable = true;
    idx_test = 0;
    while testpoint_inacceptable & (idx_test < length(allowed_list))
        idx_test = idx_test + 1;
        idx = ind(idx_test);
        x_cand = allowed_list(idx_test,:);

        not_moving_forward = ( sum( (x_cand-x_t).^2) > sum( (x_cand-x_tminus).^2) );
        step_too_large = ( sqrt(sum( (x_cand-x_t).^2)) > 2*s );
        testpoint_inacceptable = not_moving_forward | step_too_large;
    end

    if (idx_test == length(allowed_list))
        x_cand = allowed_list(ind(1),:);
    end

    trajectory(i,:) = x_cand;
    [M,idx_true] = min(sum((xpool-x_cand).^2,2));
    unusedflag(idx_true) = false;

    if sqrt(sum( (trajectory(i,:) - trajectory(i-1,:)).^2)) > 2*s*dim
        disp('oh boy...')
        disp(i)
        error('Rats are not quantum tunneling, man...')
    end

end



%% Then, get the place field activity

ff_pos = randi(2,N,dim)-1.5; % firing field positions: random edges of a centered hypercube.
ff_pos = unique(ff_pos,'rows');
while length(ff_pos) < N
    ff_pos = [ ff_pos; 2*(rand(1,dim)-0.5)];
end

ff_rate = 5;    %maximum firing rate
ff_width = 0.25;   %width of the firing rate in space.
bc_width = 0.25;

activity = zeros(N,T);
for neuron = 1:N
    activity(neuron,:) = ff_rate .* exp( - sum((trajectory - ff_pos(neuron,:)).^2,2) / (2*ff_width*ff_width) );
end

data_rat = activity' + noise*randn(size(activity')); %typical notation
data_rat(data_rat<0) = 0;


figure(1)

subplot(1,4,1)
plot(trajectory(:,1),trajectory(:,2),'-')
title("Trajectory of rat")
daspect([1 1 1])

subplot(1,4,2)
scatter(trajectory(:,1), trajectory(:,2), [], data_rat(:,1),'.')
title("Activity of place Field No.1 in arena")
daspect([1 1 1])

subplot(1,4,3)
imagesc(data_rat')
title("Activity of all place field")

subplot(1,4,4)
imagesc(cov(data_rat))
title("Covariance of the neuronal data")
daspect([1 1 1])

%% run mind on the data

mindparameters.dt = 1;
mindparameters.pca.n = 0.95;
mindparameters.dim_criterion = .95;
mindparameters.ndir = 2;
mindparameters.min_leaf_pts = 40;
mindparameters.ntrees = 100;
mindparameters.verbose = true;
mindparameters.lm.n = 2000;

mindparameters.rwd.type = 'discrete';
mindparameters.rwd.sym = 'avg';
mindparameters.rwd.all_geo = true;
mindparameters.rwd.d = 2;
mindparameters.rwd.var_scale = 0.1;

mindparameters.embed.type = 'rwe';
mindparameters.embed.d = nan;
mindparameters.embed.mode = 'mds';
mindparameters.embed.local = false;

mindparameters.learnmapping = true;
mindparameters.mapping.k = [1:10, 15:5:50];
mindparameters.mapping.lambda = [0, 10.^(-8:.5:0)];
mindparameters.mapping.mode = 'lle';
mindparameters.mapping.nfolds_lle = 10;

mindparameters.prune_lm_by_time = false;

dembed = [3]; % embedding dimensions


data = struct();
data.t = reshape(1:length(data_rat),length(data_rat),1);
data.f = data_rat;

result = struct();
result.forestdat = mindAsFunction(data, mindparameters);
result.mindparameters = mindparameters;

embedparameters = mindparameters; 
embedparameters.embed.d = dembed;
[~, result.allembed] = embedAsFunction(result.forestdat, embedparameters);
fprintf('finished running embedAsFunction\n');

y = result.allembed(1).y;

trajectory = [trajectory; trajectory];
pos_x = trajectory(result.forestdat.lm.idx,1);
pos_y = trajectory(result.forestdat.lm.idx,2);

%%
figure(1)
mi = min(y(:));
ma = max(y(:));

subplot(1,2,1)
scatter3(y(:,1), y(:,2), y(:,3),[], pos_x,'.')
daspect([1 1 1])
xlim([mi, ma])
ylim([mi,ma])
zlim([mi,ma])
title("Position X on manifold")

subplot(1,2,2)
scatter3(y(:,1), y(:,2), y(:,3),[], pos_y,'.')
daspect([1 1 1])
xlim([mi, ma])
ylim([mi,ma])
zlim([mi,ma])
title("Position Y on manifold")


figure(2)

subplot(1,3,1)
scatter(pos_x, pos_y, [], y(:,1), '.')
daspect([1 1 1])
title('dim 1 on trajectory')

subplot(1,3,2)
scatter(pos_x, pos_y, [], y(:,2), '.')
daspect([1 1 1])
title('dim 2 on trajectory')

subplot(1,3,3)
scatter(pos_x, pos_y, [], y(:,3), '.')
daspect([1 1 1])
title('dim 3 on trajectory')