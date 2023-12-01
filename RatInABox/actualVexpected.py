import numpy as np

def calculate_spatial_accuracy(actual_firing_rates, expected_firing_rates):
    # Calculate a metric (e.g., correlation) between actual and expected firing rates
    accuracy = np.corrcoef(actual_firing_rates, expected_firing_rates)[0, 1]
    return accuracy

def compare_actual_expected_firing(actual_firing_rates_envA, actual_firing_rates_envB, expected_firing_rates, balance_levels, number_of_neurons):
    # Calculate accuracy scores for each balance level and neuron
    accuracy_scores_envA = {}
    accuracy_scores_envB = {}

    for balance_level in balance_levels:
        actual_rates_A = actual_firing_rates_envA[balance_level]
        actual_rates_B = actual_firing_rates_envB[balance_level]
        for neuron_id in range(number_of_neurons):
            expected_rates = expected_firing_rates[neuron_id]
            accuracy_scores_envA[neuron_id, balance_level] = calculate_spatial_accuracy(actual_rates_A[neuron_id], expected_rates)
            accuracy_scores_envB[neuron_id, balance_level] = calculate_spatial_accuracy(actual_rates_B[neuron_id], expected_rates)

    return accuracy_scores_envA, accuracy_scores_envB
