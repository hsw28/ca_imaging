

%construct the trainingsdataset and train MIND
[X Y Z] = lorenz(28, 10, 8/3, [0 1 1.05]);

dI = 10;
dataDFF = zeros(length(X(1:dI:end)),3);
dataDFF(:,1) = X(1:dI:end);
dataDFF(:,2) = Y(1:dI:end);
dataDFF(:,3) = Z(1:dI:end);

mindparameters.dt = 1; % step distance to the past
mindparameters.pca.n = 0.99; % keep data that explains 95% variance
mindparameters.dim_criterion = 0.95;
mindparameters.ndir = 3; % number of hyperplanes tested
mindparameters.min_leaf_pts = 1000; % 100; % minimum number of leaves
mindparameters.ntrees = 100;
mindparameters.verbose = true;
mindparameters.lm.n = round(size(dataDFF,1)/2); % number of landmarks

mindparameters.rwd.type = 'discrete';
mindparameters.rwd.sym = 'avg';
mindparameters.rwd.all_geo = true;
mindparameters.rwd.d = 2;
mindparameters.rwd.var_scale = 0.1;

mindparameters.embed.type = 'rwe';
mindparameters.embed.d = nan;
mindparameters.embed.mode = 'mds';
mindparameters.embed.local = false;
mindparameters.embed.opts = statset('MaxIter',400);    %%%% Iterations for MDS

mindparameters.learnmapping = true;
if size(dataDFF,1)>100
    mindparameters.mapping.k = [5:20, 12:5:50]; % range to get optimal values
else
    mindparameters.mapping.k = [5:20, 5:3:30]; % range to get optimal values
end
mindparameters.mapping.lambda = [10.^(-8:.5:0)]; % range to get optimal values
mindparameters.mapping.mode = 'lle';
mindparameters.mapping.nfolds_lle = 10;               

mindparameters.prune_lm_by_time = false;

dembed = [3];
mindparameters.dembed = dembed;


%% run MIND
data = struct();
times = reshape(1:size(dataDFF,1), size(dataDFF,1),1)./15;
data.t = times;
data.f = dataDFF;

dat = struct();
dat.forestdat = mindAsFunction(data, mindparameters);
dat.mindparameters = mindparameters;

embedparameters = mindparameters;
embedparameters.embed.d = dembed;
[~, dat.allembed] = embedAsFunction(dat.forestdat, embedparameters);
fprintf('finished running embedAsFunction\n');

%% evaluate on test dataset

[X1 Y1 Z1] = lorenz(28, 10, 8/3, [16,11,43]);

dataDFF_full = zeros(length(X1),3);
dataDFF_full(:,1) = X1;
dataDFF_full(:,2) = Y1;
dataDFF_full(:,3) = Z1;

dim = 1;
pca_coords = dat.forestdat.pca.model.transform(dataDFF_full,  mindparameters.pca.n);
embd = dat.allembed(dim).f2m.map.transform(pca_coords);
a = dat.allembed(dim).m2f.map.transform( embd );
b = dat.forestdat.pca.model.inverse_transform(a);
reconstructed = b;

figure(1)  % cross-validated reconstruction-error plotted along the attractor
err = sqrt(sum((reconstructed-dataDFF_full).^2,2));
scatter3(X1,Y1,Z1,[],err,'.')
title('reconstruction error')


%%
figure(2)   % intrinsic coordinates plotted along the attractor
rb = jet(160);
ax1 = subplot(1,3,1);
scatter3(X1(1:10:end),Y1(1:10:end),Z1(1:10:end),[],embd(1:10:end,1),'.')
daspect([1 1 1])
ax2 = subplot(1,3,2);
scatter3(X1(1:10:end),Y1(1:10:end),Z1(1:10:end),[],embd(1:10:end,2),'.')
daspect([1 1 1])
ax3 = subplot(1,3,3);
scatter3(X1(1:10:end),Y1(1:10:end),Z1(1:10:end),[],embd(1:10:end,3),'.')
colormap(rb);
caxis([0,0.1])
daspect([1 1 1])
hlink = linkprop([ax1,ax2,ax3],{'CameraPosition','CameraUpVector'}); 
rotate3d on



