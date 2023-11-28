import numpy as np

def calculate_similarity(response_envA, response_envB):
    # Calculate a similarity metric (e.g., correlation) between responses in EnvA and EnvB
    similarity = np.corrcoef(response_envA, response_envB)[0, 1]
    return similarity

# Example data structure to store responses
responses_envA = {}  # {balance_level: [firing_rates_during_CS_and_US]}
responses_envB = {}  # Similar structure for EnvB

# After running simulations and collecting data
similarity_scores = {}
for balance_level in balance_levels:
    response_envA = responses_envA[balance_level]
    response_envB = responses_envB[balance_level]
    similarity_scores[balance_level] = calculate_similarity(response_envA, response_envB)

# Analyze similarity scores
# Higher scores indicate better learning transfer
for balance_level, score in similarity_scores.items():
    print(f"Balance Level: {balance_level}, Similarity Score: {score}")
