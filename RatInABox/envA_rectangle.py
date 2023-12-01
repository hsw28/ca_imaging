import numpy as np
from ratinabox.Environment import Environment
from ratinabox.Agent import Agent
from CombinedPlaceTebcNeurons import CombinedPlaceTebcNeurons
from trial_marker import determine_cs_us

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
    N = 80  # Adjust as needed
    firing_rates = np.zeros((N, position_data.shape[1]))

    # Check the format of responsive_distribution and adjust if necessary
    if isinstance(responsive_distribution, (float, int)):
        # If it's a single value, repeat it N times
        responsive_distribution = np.full(N, responsive_distribution)
    elif isinstance(responsive_distribution, (list, np.ndarray)):
        # Convert to numpy array if it's not already
        responsive_distribution = np.array(responsive_distribution)
        if responsive_distribution.size != N:
            raise ValueError(f"Length of responsive_distribution must be equal to the number of neurons (N={N})")
    else:
        raise TypeError("responsive_distribution must be a float, int, list, or numpy.ndarray")

    # Create CombinedPlaceTebcNeurons with the environment
    combined_neurons = CombinedPlaceTebcNeurons(agent, N, balance_distribution, responsive_distribution)

    # Initialize last CS and US times
    last_CS_time = None
    last_US_time = None

    # Iterate over the trajectory data
    for index in range(position_data.shape[1]):
        # Current timestamp
        current_time = position_data[0, index]

        # Determine if CS or US is present
        trial_marker = position_data[3, index]  # Assuming the 4th row has trial markers
        cs_present, us_present = determine_cs_us(trial_marker)

        # Update last CS/US time if necessary
        if cs_present and (last_CS_time is None or current_time > last_CS_time):
            last_CS_time = current_time
        if us_present and (last_US_time is None or current_time > last_US_time):
            last_US_time = current_time

        # Calculate time since CS and US
        time_since_CS = current_time - last_CS_time if last_CS_time is not None else -1
        time_since_US = current_time - last_US_time if last_US_time is not None else -1

        # Update agent's position and neuron states
        agent.position = np.array([position_data[1, index], position_data[2, index]])
        combined_neurons.update_state(agent.position, time_since_CS, time_since_US)


        # Store firing rates
        firing_rates[:, index] = combined_neurons.get_firing_rates()

    # Return the firing rates for further analysis
    return firing_rates
