import numpy as np
import pandas as pd
from ratinabox.Environment import Environment
from ratinabox.Agent import Agent
from your_custom_neuron_class import CombinedPlaceTebcNeurons  # Import your custom class


#modeling environment A (rectangle)
#using equation from https://www.biorxiv.org/content/10.1101/2023.10.08.561112v1.full :
'''
Place and grid cell rate maps were generated from a real exploration trajectory using
the open source Python software RatInABox. The respective activity rates are then used
to train a logistic regressor to predict the real activity of each individual neurons.
To evaluate each model performance, we computed a F1 score for each neuron using
a place input model, which penalizes both incorrect classifications of active and inactive periods.
'''

#allows me to upload my own trajectory <-- I HAVE TO SCALE THIS

# Function to determine CS/US presentation
def determine_stimulus(trial_marker):
    cs_present = 1 <= trial_marker <= 5
    us_present = 6 <= trial_marker <= 10
    return cs_present, us_present

# Convert inches to meters for EnvA
width_in_meters_A = 31 * 0.0254
height_in_meters_A = 20 * 0.0254

# Create EnvA with the converted size
envA = Environment(size=(width_in_meters_A, height_in_meters_A), boundary_conditions='solid')

# Create an agent in EnvA
agentA = Agent(environment=envA)

# Number of neurons
N = 900  # Adjust as per your data

# Create combined neurons
combined_neurons_A = CombinedPlaceTebcNeurons(environment=envA, N=N)

# Load your trajectory data for EnvA
trajectory_data_A = pd.read_csv('your_trajectory_envA.csv')

# Main simulation loop for EnvA
for index, row in trajectory_data_A.iterrows():
    current_time = row['timestamp']
    agentA.position = np.array([row['x'], row['y']])
    trial_marker = row['trial_marker']

    cs_present, us_present = determine_stimulus(trial_marker)
    combined_neurons_A.update_firing_rate(agentA.position, cs_present, us_present, current_time)

    # Record firing rates and other relevant data
    # ...
