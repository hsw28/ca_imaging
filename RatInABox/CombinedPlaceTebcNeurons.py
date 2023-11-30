import numpy as np
from ratinabox.Neurons import Neurons, PlaceCells
from tebc_response import default_tEBC_response
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


class CombinedPlaceTebcNeurons:
    def __init__(self, num_neurons, place_cells, balance, tebc_responsive_rate):
        self.num_neurons = num_neurons
        self.place_cells = place_cells
        self.balance = balance
        self.tebc_responsive_rate = tebc_responsive_rate
        self.tebc_responsive_neurons = self.assign_tebc_responsiveness()

    def assign_tebc_responsiveness(self):
        # Randomly assign a subset of neurons to be responsive to tEBC
        num_responsive = int(self.num_neurons * self.tebc_responsive_rate)
        responsive_neurons = np.random.choice([True, False], size=self.num_neurons,
                                              p=[self.tebc_responsive_rate, 1 - self.tebc_responsive_rate])
        return responsive_neurons

    def calculate_firing_rate(self, agent_position, time_since_CS, time_since_US):
        firing_rates = np.zeros(self.num_neurons)
        for i in range(self.num_neurons):
            place_response = self.place_cells.calculate_firing_rate(agent_position, i)
            tebc_response = 0
            if self.tebc_responsive_neurons[i]:
                tebc_response = default_tEBC_response(self.cell_types[i], time_since_CS, time_since_US)
            firing_rates[i] = (1 - self.balance) * place_response + self.balance * tebc_response
        return firing_rates
