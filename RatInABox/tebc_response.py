import numpy as np

#cell type responses based on manuscript Sequence of Single Neuron Changes in CA1 Hippocampus of Rabbits During Acquisition of Trace Eyeblink Conditioned Responses


def gaussian_response(time, peak_per_trial, start, end, peak_time, baseline):
    # Standard deviation is set to span 1/6th of the response window, assuming a normal distribution
    sd = (end - start) / 6
    if start <= time <= end:
        return peak_per_trial * np.exp(-((time - peak_time)**2) / (2 * sd**2)) + baseline
    return baseline

def fluctuating_response(time, start, end, baseline, fluctuation_strength, fluctuation_frequency):
    if start <= time <= end:
        fluctuation = fluctuation_strength * np.sin(time / fluctuation_frequency)
        return baseline + fluctuation
    return baseline

def bimodal_response(time, peaks_per_trial, start, end, peak_times, baseline):
    sd1 = (peak_times[0] - start) / 3
    sd2 = (end - peak_times[1]) / 3
    if start <= time <= end:
        return peaks_per_trial[0] * np.exp(-((time - peak_times[0])**2) / (2 * sd1**2)) + \
               peaks_per_trial[1] * np.exp(-((time - peak_times[1])**2) / (2 * sd2**2)) + baseline
    return baseline

def uniform_response(time, start, end, baseline, response_level):
    if start <= time <= end:
        return response_level + baseline
    return baseline

def low_level_response(time, start, end, baseline):
    if start <= time <= end:
        return baseline
    return baseline

# Define the response profiles for each cell type based on the histogram data
response_profiles = {
    1: {'response_func': lambda t: gaussian_response(t, peak_per_trial=0.005, start=-50, end=300, peak_time=50, baseline=0.005)},
    2: {'response_func': lambda t: gaussian_response(t, peak_per_trial=0.0025, start=-50, end=150, peak_time=50, baseline=0.005)},
    3: {'response_func': lambda t: fluctuating_response(t, start=-100, end=900, baseline=0.03, fluctuation_strength=0.001, fluctuation_frequency=50)},
    4: {'response_func': lambda t: bimodal_response(t, peaks_per_trial=[0.001, 0.00075], start=-50, end=450, peak_times=[-25, 200], baseline=0.025)},
    5: {'response_func': lambda t: gaussian_response(t, peak_per_trial=0.00125, start=0, end=150, peak_time=50, baseline=0.005)},
    6: {'response_func': lambda t: uniform_response(t, start=-100, end=900, baseline=0.018, response_level=0.002)},
    7: {'response_func': lambda t: gaussian_response(t, peak_per_trial=0.001, start=-50, end=350, peak_time=100, baseline=0.01)},
    8: {'response_func': lambda t: low_level_response(t, start=-100, end=900, baseline=0.015)},
}

# Example usage:
cell_type = 1
time_since_CS = 100  # Time in ms
# Assuming the time_since_US doesn't affect the response for simplicity
# If it does, you will need to modify the functions to incorporate the effect of US
firing_rate = response_profiles[cell_type]['response_func'](time_since_CS)
print(firing_rate)
