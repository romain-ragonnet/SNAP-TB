import numpy as np


def generate_an_activation_profile(age, params):
    """
    This method generates the outcome of the LTBI episode. Return a time to TB activation (in days. Activation will effectively
    occur if the date of death of the individual is later than the date of activation.
    The method and parameterisation are based on Model 6 from Ragonnet et al 2017.
    """
    stages = ['_child', '_teen', '_adult']
    if age < 5.:
        stage_index = 0
    elif age < 15.:
        stage_index = 1
    else:
        stage_index = 2
    # high risk or low risk
    high_risk = np.random.binomial(n=1, p=(1. - params['g' + stages[stage_index]]), size=1)
    if high_risk == 1:
        # use age-specific epsilon to generate an activation time
        time_to_activation = np.random.exponential(scale=1./params['epsi' + stages[stage_index]], size=1)
    else:
        # we use the reactivation rates that vary with age
        if stage_index == 0:
            time_to_activation = np.random.exponential(scale=1./params['nu' + stages[0]], size=1)
            if time_to_activation/365.25 + age > 5.:  # we use the reactivation rates for the [5-15] age-category
                time_to_activation = (5. - age)*365.25 +\
                                     np.random.exponential(scale=1./params['nu' + stages[1]], size=1)
                if time_to_activation / 365.25 + age > 15.:  # we use the reactivation rates for the [5-15] age-category
                    time_to_activation = (15. - age) * 365.25 + \
                                         np.random.exponential(scale=1. / params['nu' + stages[2]], size=1)
        elif stage_index == 1:
            time_to_activation = np.random.exponential(scale=1. / params['nu' + stages[1]], size=1)
            if time_to_activation / 365.25 + age > 15.:  # we use the reactivation rates for the [5-15] age-category
                time_to_activation = (15. - age) * 365.25 + \
                                     np.random.exponential(scale=1. / params['nu' + stages[2]], size=1)
        else:
            time_to_activation = np.random.exponential(scale=1. / params['nu' + stages[2]], size=1)

    return time_to_activation


if __name__ == "__main__":
    n_active = 0
    for _ in range(100000):
        age = np.random.uniform(low=0., high=70.)
        time_to_act = generate_an_activation_profile(age)
        print "Age: " + str(age) + " / Time in years: " + str(round(time_to_act/365.25))

        time_to_live = 365.25*(70 - age)
        if time_to_act < time_to_live:
            n_active += 1

    print n_active
