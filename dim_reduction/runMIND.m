function f = runMIND(CAdata, time)
%% HOW TO run mind on the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(time,2) > size(time,1)
  time = time';
end



if size(time,1) ~= size(CAdata,1)
  CAdata = CAdata';
end


data = struct();

data.t = time; % [n x 1]
data.f = CAdata; %[n x cell #]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%PARAMETERS
%%FROM HPC PAPER
%In brief, we first learned a generative model of transition probabilities from population activity s(t) = [s1(t), …, sN(t)]
%of N neurons at time 0 < t < T, to the activity s(t + Δt) using the previously developed random forest method13 with a
%few modifications. First, when splitting the neural state space into regions using a set of hyperplanes organized in a
%decision tree, we assessed 20 random hyperplane orientations at every node of the tree and selected the orientation that
%best split the data. This improved performance with the large numbers of neurons typically encountered in calcium imaging.
%Second, we set the minimum number of leaves in each random tree to 500. Third, to define transitions, we considered all
%states Δt = 67 ms apart (one frame at a 15-Hz frame rate). Fourth, we fit manifolds to all data points, not only a subset
%of landmarks. All other hyperparameters were chosen as previously described

dembed = [6]; % embedding dimensions

mindparameters.dt = 1; % step distance to the past
mindparameters.pca.n = 0.95; %how many principal components to retain (or fraction of variance)
mindparameters.dim_criterion = .95; %I think this is the came as pca.n but for hybrid pca
mindparameters.ndir = 20; %number of hyperplanes tested
mindparameters.min_leaf_pts = 500; %minimum number of leaves
mindparameters.ntrees = 100; %number of trees --  no risk of overfitting by increasing the number of trees
mindparameters.verbose = true;
mindparameters.lm.lmf = 1; %landmark fraction
mindparameters.lm.n = length(data.t); %number of landmarks?? can be given as, for ex, round(size(dataDFF,1)/2). we dont want to use landmarks so I think this should be equal to number of data points?



mindparameters.rwd.type = 'discrete'; %type of random walk ('continuous' or 'discrete')
mindparameters.rwd.sym = 'avg'; % how to symmetrize global distances ('avg' or 'min')
mindparameters.rwd.all_geo = true; % if true, all distances will be geodesic distances. otherwise,
% distances between connected points will be ideal local distances and
% distances between non-connected points will be filled in with
% geodesic distances.
mindparameters.rwd.d = 2; % dimensionality of space in which random walk is performed
% (only used for continuous random walk)
mindparameters.rwd.var_scale = 0.1; % variance of diffusion kernel, expressed as fraction of maximum
% possible variance (only used for continuous random walk; shouldn't
% matter much)

mindparameters.embed.type = 'rwe';  %no idea -- told to ignore
mindparameters.embed.d = dembed;       %what MIND embeds into. In the scripts it is set twice. mindparameters.embed.d = nan; just initializes it
mindparameters.embed.mode = 'mds';  %options are mds for non-classic multidimensional scaling
%cmds for classic, or heirarchical
mindparameters.embed.local = false; %local embedding for MDS?
mindparameters.embed.opts = statset('MaxIter',400);    % Iterations for MDS, optional

mindparameters.learnmapping = true;   %learn mapping from global pca space to manifold
%mindparameters.mapping.k = [5:20, 15:5:50]; % choices of k for k nearest neighbors
%another example has it like so:
if size(data,1)>100
        mindparameters.mapping.k = [5:20, 12:5:50]; % range to get optimal values
     else
        mindparameters.mapping.k = [5:20, 5:3:30]; % range to get optimal values
     end
mindparameters.mapping.lambda = [0, 10.^(-8:.5:0)]; % choices of regularization parameter
mindparameters.mapping.mode = 'lle'; %for local linear embedding, other options are w-n kernel, tlle which is lle with time blocked cv?, or gp
mindparameters.mapping.nfolds_lle = 10; %choices of regularization parameter
mindparameters.prune_lm_by_time = false; % prune landmarks so they are separated in time, if requested
%mindparameters.prune_lm_by_time = true;


%%%%%%%%%%

result = struct();
result.forestdat = mindAsFunction(data, mindparameters);
result.mindparameters = mindparameters;
embedparameters = mindparameters;
embedparameters.embed.d = dembed;
[~, result.allembed] = embedAsFunction(result.forestdat, embedparameters);
fprintf('finished running embedAsFunction\n');

f= result.allembed;
y = result.allembed(1).y;

figure

scatter3(y(:,1), y(:,2), y(:,3), 5, [1:1:length(y)]);
