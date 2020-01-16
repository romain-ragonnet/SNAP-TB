import numpy as np
# from math import exp
# from scipy.optimize import fsolve
# import dill

#
# def generate_time_to_activation(age, params, index=None):
#     """
#     Randomly generates a time form infection to activation that depends on individual's age, hiv status and diabetes.
#     We use the inverse transform sampling method, based on the cdf of the times to activation.
#     index provided the list index of the stored time that will be used
#     """
#     if age < 5.:
#         string = '_child'
#     elif age < 15.:
#         string = '_teen'
#     else:
#         string = '_adult'
#
#     if not params['generate_activation_times']:  # read a time from the loaded dictionary
#         # ind = np.random.randint(low=0, high=100000)
#         return params['activation_times_dic'][string][index]
#     else:
#         activation_params = {'epsi': None, 'kappa': None, 'nu': None, 'mu': None}
#         activation_params['mu'] = 1 / (params['life_expectancy'] * 365.25)
#
#         for param in activation_params.keys():
#             if param != 'mu':
#                 activation_params[param] = params[param + string]
#
#         # inverse transform sampling
#         activation_params['u'] = np.random.rand()
#         if activation_params['u'] >= get_overall_proportion(activation_params):
#             time = None
#         else:
#             time = fsolve(proportion_of_active_tb_over_time, 100., args=activation_params)
#
#         return time
#
# def get_overall_proportion(a_p):
#     """
#     Calculates the total risk of activation as provided in Ragonnet et al., Epidemics 2017, Table S5, Model 4.
#     a_p stands for activation_params coming from 'generate_time_to_activation'.
#     """
#     num = a_p['epsi'] * (a_p['nu'] + a_p['mu']) + a_p['kappa'] * a_p['nu']
#     deno = (a_p['kappa'] + a_p['epsi'] + a_p['mu']) * (a_p['nu'] + a_p['mu'])
#     return num/deno
#
# def proportion_of_active_tb_over_time(t, a_p):
#     """
#     The activation function as provided in Ragonnet et al., Epidemics 2017, Table S4, Model 4.
#     a_p stands for activation_params coming from 'generate_time_to_activation'.
#     """
#     frac_1 = a_p['epsi'] / (a_p['kappa'] + a_p['epsi'] + a_p['mu'])
#     frac_2 = a_p['kappa'] * a_p['nu'] / ((a_p['nu'] - a_p['kappa'] - a_p['epsi']) * (a_p['kappa'] + a_p['epsi'] + a_p['mu']))
#     frac_3 = a_p['kappa'] * a_p['nu'] / ((a_p['nu'] + a_p['mu']) * (a_p['kappa'] + a_p['epsi'] - a_p['nu']))
#     rate_1 = -(a_p['kappa'] + a_p['epsi'] + a_p['mu'])
#     rate_2 = -(a_p['nu'] + a_p['mu'])
#
#     i_t = (frac_1 + frac_2) * (1. - exp(rate_1 * t)) + frac_3 * (1. - exp(rate_2 * t))
#     return i_t - a_p['u']
#
# def generate_and_store_times_to_activation(n_times, params):
#     """
#     Generate and save times to activation in a pickle object that could be reused later.
#     :param n_times: number of infection histories generated for each age category
#     :param params: a dictionary of parameters containing the latency parameters
#     """
#     print str(n_times) + " new activation times are being generated."
#     params['generate_activation_times'] = True
#     example_age = {'_child':3., '_teen': 10., '_adult': 30.}
#     activation_times = {'_child': [], '_teen': [], '_adult': []}
#
#     for age_cat, age in example_age.iteritems():
#         for i in range(n_times):
#             activation_times[age_cat].append(generate_time_to_activation(age, params))
#
#     file_name = "pickled_activation_times.pickle"
#     file_stream = open(file_name, "wb")
#     dill.dump(activation_times, file_stream)
#     file_stream.close()
#     print "Activation times have been saved successfully."

def generate_an_activation_profile(age, params):
    """
    This method generates the outcome of the LTBI episode. Return a time to TB activation (in days. Activation will effectively
    occur if the date of death of the individual is later than the date of activation.
    The method and parameterisation are based on Model 6 from Ragonnet et al 2017.
    """
    # params = {
    #     'epsi_child': 1.9e-2,
    #     'epsi_teen': 1.4e-2,
    #     'epsi_adult': 5.6e-3,
    #     'nu_child': 3.4e-11,
    #     'nu_teen': 6.4e-6,
    #     'nu_adult': 3.3e-6,
    #     'g_child': 0.65,
    #     'g_teen': 0.81,
    #     'g_adult': 0.95
    # }
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
        # use epsilon to generate an activation time
        time_to_activation = np.random.exponential(scale=1./params['epsi' + stages[stage_index]], size=1)
    else:
        # we use the reactivation rates that vary with age
        time_to_activation = float("inf")
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
