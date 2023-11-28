# Import necessary modules and custom classes
from envA_rectangle import simulate_envA
from envB_oval import simulate_envB
from CombinedPlaceTebcNeurons import CombinedPlaceTebcNeurons
from trial_marker import determine_cs_us
from learningTransfer import assess_learning_transfer
from actualVexpected import compare_actual_expected_firing

# Parameters
balance_levels = [0.2, 0.4, 0.6, 0.8, 1.0]  # Example balance levels

# Results containers
learning_transfer_results = {}
spatial_coding_accuracy_results = {}

# Main simulation loop
for balance in balance_levels:
    # Update balance parameter in CombinedPlaceTebcNeurons
    CombinedPlaceTebcNeurons.set_balance(balance)

    # Simulate in Environment A
    firing_rates_envA, position_data_envA = simulate_envA()

    # Simulate in Environment B
    firing_rates_envB, position_data_envB = simulate_envB()

    # Assess learning transfer
    learning_transfer_score = assess_learning_transfer(firing_rates_envA, firing_rates_envB)
    learning_transfer_results[balance] = learning_transfer_score

    # Compare actual vs. expected firing rates for spatial coding accuracy
    accuracy_envA = compare_actual_expected_firing(firing_rates_envA, position_data_envA)
    accuracy_envB = compare_actual_expected_firing(firing_rates_envB, position_data_envB)
    spatial_coding_accuracy_results[balance] = (accuracy_envA, accuracy_envB)

    # Additional analysis (e.g., decoding accuracy) can be added here

# Output results
for balance, transfer_score in learning_transfer_results.items():
    print(f"Balance Level: {balance}, Learning Transfer Score: {transfer_score}")

for balance, (accuracy_A, accuracy_B) in spatial_coding_accuracy_results.items():
    print(f"Balance Level: {balance}, Spatial Accuracy EnvA: {accuracy_A}, EnvB: {accuracy_B}")

# Further analysis to determine optimal balance
# ...

# Iterate and refine based on findings
# ...
