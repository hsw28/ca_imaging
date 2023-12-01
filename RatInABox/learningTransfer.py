import numpy as np

def calculate_similarity(response_envA, response_envB):
    # Ensure that both arrays are of the same length
    min_length = min(len(response_envA), len(response_envB))
    response_envA = response_envA[:min_length]
    response_envB = response_envB[:min_length]

    # Calculate a similarity metric (e.g., correlation) between responses in EnvA and EnvB
    similarity = np.corrcoef(response_envA, response_envB)[0, 1]
    return similarity

def assess_learning_transfer(response_envA, response_envB, balance_value, responsive_value):
    # Calculate the similarity score for the given balance and responsive values
    similarity_score = calculate_similarity(response_envA, response_envB)
    return similarity_score


#i dont think this is actually want i want to measure
