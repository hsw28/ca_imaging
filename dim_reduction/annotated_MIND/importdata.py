import scipy.io
import numpy as np
from scipy.io import savemat
# to create MIND embeddings
from mind import MIND
# For plotting purposes
import matplotlib.pyplot as plt

##FROM HPC PAPER
#In brief, we first learned a generative model of transition probabilities from population activity s(t) = [s1(t), …, sN(t)]
#of N neurons at time 0 < t < T, to the activity s(t + Δt) using the previously developed random forest method13 with a
#few modifications. First, when splitting the neural state space into regions using a set of hyperplanes organized in a
#decision tree, we assessed 20 random hyperplane orientations at every node of the tree and selected the orientation that
#best split the data. This improved performance with the large numbers of neurons typically encountered in calcium imaging.
#Second, we set the minimum number of leaves in each random tree to 500. Third, to define transitions, we considered all
#states Δt = 67 ms apart (one frame at a 15-Hz frame rate). Fourth, we fit manifolds to all data points, not only a subset
#of landmarks. All other hyperparameters were chosen as previously described


### my covariance matrix has a high Condition number (close to singular -- rows are dependant)

# Load the .mat file
mat_data = scipy.io.loadmat('/Users/Hannah/Programming/data_eyeblink/test3.mat')

# Access the MATLAB matrix from the loaded data
matlab_matrix = mat_data['test3']

# Convert the MATLAB matrix to a NumPy array
numpy_array = np.array(matlab_matrix)

# Now, numpy_array contains the MATLAB matrix as a NumPy array
X = numpy_array
Xn=np.transpose(X)
## WANT FIRST VALUE TO BE NUMBER OF TIME POINTS, SECOND TO BE NUMBER OF CELLS == NxD
X = Xn[:-1,:]
Xp = Xn[1:,:]



# parameters
#param_forest_nb_trees: number of tree in the MIND random forest
#param_tree_nb_leaf: min number of states in a node below which the node is considered a leaf and won't be split
#param_tree_nb_v: number of random candidates {v} (hyperplane vectors) for node partition
#param_tree_nb_c: number of random candidates {c} (intercepts) for node partition
#param_emb_dim: dimensionality of the embedding space.
#param_emb_fraction: fraction of all states (X) used to create the embedding (default = 1.)
#param_emb_p_threshold: min probability under which probability are considered equal to 0.
#param_emb_fraction_vote: min fraction of non-zero probaiblities across trees for p(x2|x2) to not be considered null (see explanation above). default value = 0.5.
#param_emb_mds_only: boolean indicating whether restricting the calculation of the coordinate to multidimensional scaling (default). if False, additional gradient-based optimization is performed.

param_forest_nb_trees = 10 #is 500 in hpc paper
param_tree_nb_leaf = 400 #is 500 in hpc paper
param_tree_nb_v = 5 #20???
param_tree_nb_c = 7

param_emb_dim = 5
param_emb_fraction = 1.
param_emb_p_threshold = 0. #i have to make this super small to not get inf distances
param_emb_fraction_vote = 0.
param_emb_mds_only = True #was false
#///////////////////////////////////////////


# create a MIND object with nb_trees trees
#  and tree parition parameters nb_leaf, nb_v, and nb_c.
mind = MIND(nb_trees = param_forest_nb_trees,
            nb_leaf = param_tree_nb_leaf,
            nb_v = param_tree_nb_v,
            nb_c = param_tree_nb_c)


# fit MIND object to our synthetic data.
mind.fit(X, Xp,
         emb_dim = param_emb_dim, emb_fraction = param_emb_fraction,
         p_threshold = param_emb_p_threshold, fraction_vote = param_emb_fraction_vote,
         mds_only = False)

savemat('BRAIN_data.mat', {'my_variable': mind.emb_Y})
