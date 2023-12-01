def determine_cs_us(trial_marker):
    """
    Determines if the conditioned stimulus (CS) or unconditioned stimulus (US) should be presented.

    Args:
    trial_marker (int): The marker indicating the trial phase from the position file.
                        0 for intertrial, 1-5 for CS, 6-10 for US.

    Returns:
    tuple: (cs_present, us_present) indicating the presence of CS and US.
    """
    cs_present = 1 <= trial_marker <= 5
    us_present = 6 <= trial_marker <= 10

    return cs_present, us_present
