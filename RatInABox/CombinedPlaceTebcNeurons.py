import numpy as np
from ratinabox.Neurons import Neurons, PlaceCells
from tebc_response import response_profiles

'''
Python class template for CombinedPlaceTebcNeurons that integrates both place cell and tEBC
cell functionalities. This class is designed to be used with the RatInABox framework.
- It includes a balance parameter to adjust the contribution of place cell activity versus
tEBC cell activity for each neuron.
- also includes tebc_responsive_rate that specifies the percentage of neurons that are responsive to tEBC signals.

# Example usage
num_neurons = 100
balance = 0.5  # Example balance factor
tebc_responsive_rate = 0.6  # Example: 60% of neurons are tEBC-responsive
combined_neurons = CombinedPlaceTebcNeurons(num_neurons, place_cells, balance, tebc_responsive_rate)

'''


class CombinedPlaceTebcNeurons(Neurons):
    def __init__(self, agent, N, balance_distribution, responsive_distribution):
        self.agent = agent
        self.num_neurons = N

        # Define parameters for PlaceCells and initialize them
        place_cells_params = {
            "n": N,
            "name": "PlaceCells",
            "description": "gaussian",
            "widths": 0.20,
            "place_cell_centres": None,
            "wall_geometry": "geodesic",
            "min_fr": 0,
            "max_fr": 1,
        }
        self.place_cells = PlaceCells(agent, place_cells_params)

        # Set balance and responsiveness distributions
        self.balance_distribution = balance_distribution
        self.responsive_distribution = responsive_distribution

        # Assign TEBC responsiveness and cell types
        self.tebc_responsive_neurons, self.cell_types = self.assign_tebc_responsiveness_and_types()
        self.firing_rates = np.zeros(N)

    def assign_tebc_responsiveness_and_types(self):
        # Check if responsive_distribution is a single value or an array
        if isinstance(self.responsive_distribution, (float, int)):
            responsive_probs = np.full(self.num_neurons, self.responsive_distribution)
        else:
            responsive_probs = np.array(self.responsive_distribution)
            if responsive_probs.ndim != 1 or len(responsive_probs) != self.num_neurons:
                raise ValueError("responsive_distribution must be a 1D array of length num_neurons")
        responsive_probs = np.clip(responsive_probs, 0, 1)
        responsive_neurons = np.random.rand(self.num_neurons) < responsive_probs

        cell_type_probs = [0.051, 0.032, 0.373, 0.155, 0.199, 0.050, 0.093, 0.047]
        cell_types = np.random.choice(range(1, 9), size=self.num_neurons, p=cell_type_probs)
        return responsive_neurons, cell_types

    def update_state(self, agent_position, time_since_CS, time_since_US):
        self.agent.position = agent_position
        self.place_cells.update()

        for i in range(self.num_neurons):
            # Check if the history for each neuron is populated
            if len(self.place_cells.history['firingrate']) > i and len(self.place_cells.history['firingrate'][i]) > 0:
                place_response = self.place_cells.history['firingrate'][i][-1]
            else:
                place_response = 0  # default value if history not populated

            tebc_response = 0
            if self.tebc_responsive_neurons[i]:
                cell_type = self.cell_types[i]
                response_func = response_profiles[cell_type]['response_func']
                tebc_response = response_func(time_since_CS)

            self.firing_rates[i] = (1 - self.balance_distribution[i]) * place_response + self.balance_distribution[i] * tebc_response

    def calculate_firing_rate(self, agent_position, time_since_CS, time_since_US):
        firing_rates = np.zeros(self.num_neurons)
        for i in range(self.num_neurons):
            place_response = self.firing_rates[i]  # Directly use the updated firing rates
            tebc_response = 0
            if self.tebc_responsive_neurons[i]:
                cell_type = self.cell_types[i]
                response_func = response_profiles[cell_type]['response_func']
                tebc_response = response_func(time_since_CS, time_since_US)
            firing_rates[i] = (1 - self.balance_distribution[i]) * place_response + self.balance_distribution[i] * tebc_response
        return firing_rates

    def get_firing_rates(self):
        # Return the current firing rates of all neurons
        return self.firing_rates
