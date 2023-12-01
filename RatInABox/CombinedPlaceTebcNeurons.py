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
    def __init__(self, agent, num_neurons, balance_distribution, responsive_distribution):
        super().__init__(agent, num_neurons)
        self.place_cells = PlaceCells(agent.environment, num_neurons)
        self.balance_distribution = balance_distribution
        self.responsive_distribution = responsive_distribution
        self.tebc_responsive_neurons = self.assign_tebc_responsiveness()

    def assign_tebc_responsiveness_and_types(self):
        # Assign a subset of neurons to be responsive to tEBC and their cell types
        responsive_neurons = np.random.choice([True, False], size=self.num_neurons,
                                              p=[self.tebc_responsive_rate, 1 - self.tebc_responsive_rate])
        cell_type_probs = [0.051, 0.032, 0.374, 0.155, 0.199, 0.050, 0.093, 0.047]
        cell_types = np.random.choice(range(1, 9), size=self.num_neurons, p=cell_type_probs)
        return responsive_neurons, cell_types


    def update_state(self, agent_position, time_since_CS, time_since_US):
        # Update the state of the place cells based on the agent's position
        self.place_cells.update_state(agent_position)

        # Iterate over each neuron to update its state
        for i in range(self.num_neurons):
            # Calculate the place cell response
            place_response = self.place_cells.calculate_firing_rate(agent_position, i)

            # Initialize tEBC response
            tebc_response = 0

            # Check if the neuron is tEBC-responsive and update its state
            if self.tebc_responsive_neurons[i]:
                # Determine the cell type for this neuron
                cell_type = ...  # Assign the cell type based on your criteria

                # Get the response function for the cell type
                response_func = response_profiles[cell_type]['response_func']

                # Calculate the tEBC response based on time since CS and US
                tebc_response = response_func(time_since_CS)

            # Combine place and tEBC responses based on the balance factor
            combined_response = (1 - self.balance_distribution[i]) * place_response + self.balance_distribution[i] * tebc_response

            # Update the firing rate of the neuron
            self.firing_rates[i] = combined_response



    def calculate_firing_rate(self, agent_position, time_since_CS, time_since_US):
        firing_rates = np.zeros(self.num_neurons)
        for i in range(self.num_neurons):
            place_response = self.place_cells.calculate_firing_rate(agent_position, i)
            tebc_response = 0
            if self.tebc_responsive_neurons[i]:
                cell_type = ...  # Determine the cell type here
                response_func = response_profiles[cell_type]['response_func']
                tebc_response = response_func(time_since_CS, time_since_US)
            firing_rates[i] = (1 - self.balance_distribution[i]) * place_response + self.balance_distribution[i] * tebc_response
        return firing_rates
