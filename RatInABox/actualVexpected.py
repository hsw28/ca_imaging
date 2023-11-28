import numpy as np

def calculate_spatial_accuracy(actual_firing_rates, expected_firing_rates):
    # Calculate a metric (e.g., correlation) between actual and expected firing rates
    accuracy = np.corrcoef(actual_firing_rates, expected_firing_rates)[0, 1]
    return accuracy

# Example data structures to store firing rates
actual_firing_rates_envA = {}  # {balance_level: [firing_rates]}
actual_firing_rates_envB = {}  # Similar structure for EnvB
expected_firing_rates = {}     # {neuron_id: [expected_firing_rates_based_on_position]}

# After running simulations and collecting data
accuracy_scores_envA = {}
accuracy_scores_envB = {}
for balance_level in balance_levels:
    actual_rates_A = actual_firing_rates_envA[balance_level]
    actual_rates_B = actual_firing_rates_envB[balance_level]
    for neuron_id in range(number_of_neurons):
        expected_rates = expected_firing_rates[neuron_id]
        accuracy_scores_envA[neuron_id, balance_level] = calculate_spatial_accuracy(actual_rates_A[neuron_id], expected_rates)
        accuracy_scores_envB[neuron_id, balance_level] = calculate_spatial_accuracy(actual_rates_B[neuron_id], expected_rates)

# Analyze accuracy scores
# Higher scores indicate better spatial coding accuracy
for (neuron_id, balance_level), score in accuracy_scores_envA.items():
    print(f"EnvA - Neuron: {neuron_id}, Balance Level: {balance_level}, Accuracy Score: {score}")
for (neuron_id, balance_level), score in accuracy_scores_envB.items():
    print(f"EnvB - Neuron: {neuron_id}, Balance Level: {balance_level}, Accuracy Score: {score}")
