import numpy as np
from ratinabox.Environment import Environment
from ratinabox.Agent import Agent
from CombinedPlaceTebcNeurons import CombinedPlaceTebcNeurons

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


def simulate_envA(position_data, balance_distribution, responsive_distribution):


    # Define environment parameters for a rectangular environment
    env_params = {
        'boundary': [[0, 0], [0, 20 * 0.0254], [31 * 0.0254, 20 * 0.0254], [31 * 0.0254, 0]],  # Adjust the coordinates as needed
        'boundary_conditions': 'solid'
    }
    env = Environment(params=env_params)

     # Create an agent in the environment
    agent = Agent(env)  # Pass the environment object 'env' to the Agent

    # Number of neurons
    N = 900  # Adjust as needed

    # Create CombinedPlaceTebcNeurons with the environment
    combined_neurons = CombinedPlaceTebcNeurons(agent, N, balance_distribution, responsive_distribution)


    # Initialize an array to store firing rates
    firing_rates = np.zeros((N, position_data.shape[1]))

    # Iterate over the trajectory data
    for index in range(position_data.shape[1]):  # Iterate over columns
        # Update agent's position
        agent.position = np.array([position_data[1, index], position_data[2, index]])  # x and y

        # Update neuron's state
        combined_neurons.update_state(agent.position, position_data[0, index])  # timestamp

        # Store firing rates
        firing_rates[:, index] = combined_neurons.get_firing_rates()

    # Return the firing rates for further analysis
    return firing_rates
