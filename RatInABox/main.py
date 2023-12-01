import numpy as np
import scipy.io
import argparse
import itertools
from envA_rectangle import simulate_envA
from envB_oval import simulate_envB
from CombinedPlaceTebcNeurons import CombinedPlaceTebcNeurons
from trial_marker import determine_cs_us
from learningTransfer import assess_learning_transfer
from actualVexpected import compare_actual_expected_firing

"""
Simulation Script for Neuronal Firing Rate Analysis
Note: advise using a conda environment:
    conda create -n ratinabox python=3.9
    conda activate ratinabox
    conda install numpy
    conda install scipy
    conda install matplotlib
    export PYTHONPATH="${PYTHONPATH}:/Users/Hannah/Programming/RatInABox"
    pip install shapely


Usage:
    python main.py [--balance_values BALANCE_VALUES] [--balance_dist BALANCE_DIST] [--balance_std BALANCE_STD]
                   [--responsive_values RESPONSIVE_VALUES] [--responsive_type RESPONSIVE_TYPE]

Arguments:
    --balance_values  : Comma-separated list of balance values or means for Gaussian distribution.
                        Example: --balance_values 0.3,0.5,0.7
                        If not provided, a default value of 0.5 is used.
    --balance_dist    : Specifies the type of distribution for the balance factor.
                        Options are 'fixed' and 'gaussian'.
                        Default is 'fixed'.
    --balance_std     : Standard deviation for the Gaussian distribution of the balance factor.
                        Only used if --balance_dist is set to 'gaussian'.
                        Default value is 0.1.
    --responsive_values: Comma-separated list of responsive rates or probabilities for distributions.
                         Example: --responsive_values 0.4,0.6,0.8
                         If not provided, a default value of 0.5 is used.
    --responsive_type : Type of distribution for the responsive rate.
                        Options are 'fixed', 'binomial', 'normal', 'poisson'.
                        Default is 'fixed'.

Examples:
    python main.py --balance_values 0.3,0.5,0.7 --balance_dist gaussian --balance_std 0.1
                   --responsive_values 0.4,0.6,0.8 --responsive_type binomial

    python main.py --balance_values 0.5 --balance_dist fixed --responsive_values 0.5 --responsive_type fixed

Description:
    The script conducts simulations to evaluate how different configurations of balance factors and responsive rates affect neuronal firing patterns. Balance can be set as a fixed value or as a mean for a Gaussian distribution. The responsive rate determines the proportion of neurons responsive to tEBC signals and can be set as a fixed value or sampled from specified distributions.

    The script loads position data from a MATLAB file, performs simulations in two environments, and assesses learning transfer and spatial coding accuracy. The script supports a grid search over multiple balance and responsive rate values, allowing a comprehensive analysis of various parameter combinations. Results are printed to the console.

Requirements:
    - Ensure all necessary modules and custom classes are correctly imported and configured.
    - Replace 'path_to_your_matlab_file.mat' with the actual path to your MATLAB file.
    - Adjust environment settings and neuron parameters as needed in the script.
"""



def parse_list(arg_value):
    if isinstance(arg_value, list):
        return [float(item) for item in arg_value]
    else:
        return [float(item) for item in arg_value.split(',')]


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Simulation Script for Neuronal Firing Rate Analysis')
parser.add_argument('--balance_values', type=parse_list, help='List of balance values or means for Gaussian distribution')
parser.add_argument('--balance_dist', choices=['fixed', 'gaussian'], default='fixed', help='Distribution type for balance')
parser.add_argument('--balance_std', type=float, default=0.1, help='Standard deviation for Gaussian balance distribution')
parser.add_argument('--responsive_values', type=parse_list, help='List of responsive rates or probabilities for distributions')
parser.add_argument('--responsive_type', choices=['fixed', 'binomial', 'normal', 'poisson'], default='fixed', help='Type of distribution for responsive rate')
args = parser.parse_args()

def get_distribution_values(dist_type, params, size):
    if dist_type == 'fixed':
        return np.full(size, params[0])
    elif dist_type == 'gaussian':
        mean, std = params
        return np.clip(stats.norm(mean, std).rvs(size=size), 0, 1)
    elif dist_type == 'binomial':
        p = params[0]
        return np.random.binomial(1, p, size=size)
    elif dist_type == 'normal':
        mean, std = params
        return np.clip(stats.norm(mean, std).rvs(size=size), 0, 1)
    elif dist_type == 'poisson':
        lam = params[0]
        return np.clip(stats.poisson(lam).rvs(size=size), 0, 1)

# Load MATLAB file and extract position data
matlab_file_path = '/Users/Hannah/Programming/data_eyeblink/rat314/ratinabox_data/pos314.mat'  # Replace with your MATLAB file path
data = scipy.io.loadmat(matlab_file_path)
position_data_envA = data['envA314_522']  # Adjust variable name as needed
position_data_envB = data['envB314_524']  # Adjust variable name as needed

# Set parameters
# Set parameters
num_neurons = 80
balance_values = parse_list(args.balance_values) if args.balance_values else [0.5]
responsive_values = parse_list(args.responsive_values) if args.responsive_values else [0.5]

# Perform grid search over balance and responsive rates
for balance_value, responsive_val in itertools.product(balance_values, responsive_values):
    balance_distribution = get_distribution_values(args.balance_dist, [balance_value, args.balance_std], num_neurons)
    responsive_distribution = get_distribution_values(args.responsive_type, [responsive_val], num_neurons)

    # Simulate in Environment A and Environment B
    response_envA = simulate_envA(position_data_envA, balance_distribution, responsive_distribution)
    response_envB = simulate_envB(position_data_envB, balance_distribution, responsive_distribution)

    # Assess learning transfer and other metrics
    similarity_score = assess_learning_transfer(response_envA, response_envB, balance_value, responsive_val)
    print(f"Balance: {balance_value}, Responsive Rate: {responsive_val}, Learning Transfer: {similarity_score}")
