import numpy as np

def calculate_similarity(response_envA, response_envB):
    # Calculate a similarity metric (e.g., correlation) between responses in EnvA and EnvB
    similarity = np.corrcoef(response_envA, response_envB)[0, 1]
    return similarity

def assess_learning_transfer(responses_envA, responses_envB, balance_levels):
    # Calculate similarity scores for each balance level
    similarity_scores = {}
    for balance_level in balance_levels:  # Use a different name for the loop variable
        response_envA = responses_envA[balance_level]
        response_envB = responses_envB[balance_level]
        similarity_scores[balance_level] = calculate_similarity(response_envA, response_envB)

    return similarity_scores
