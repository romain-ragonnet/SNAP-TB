from numpy import random, nonzero, exp
import tb_activation

class individual:
    """
    This class defines the individuals
    """
    def __init__(self, id, household_id, dOB):
        # individual non TB-specific characteristics
        self.id = id  # integer
        self.household_id = household_id     # integer
        self.dOB = dOB   # integer corresponding to the date of birth

        # variables below are not used for the moment
        self.diabetes = False   # boolean
        self.hiv = False    # boolean

        # individual TB-specific characteristics
        self.vaccinated = False
        self.vaccine_immunity_duration = 15.  # in years
        self.ltbi = False
        self.tb_strain = None  # 'ds' or 'mdr'
        self.active_tb = False
        self.tb_organ = None
        self.tb_detected = False
        self.tb_treated = False
        self.die_with_tb = False
        self.contacts_while_tb = {'household': set([]), 'school': set([]), 'workplace': set([]), 'community': set([])}

        # individual group characteristics (school and work related)
        self.group_ids = []  # list of ids of the groups to which the individual is currently belonging
        self.is_ever_gonna_work = False

        # programmed events
        self.programmed = {'death': None, 'leave_home': None, 'go_to_school': None, 'leave_school': None,
                           'leave_work': None}

    def set_dOB(self, age, current_time, time_step=None):
        """
        Given a specific age (in years) at the time of initialisation, defines the date of birth of an individual.
        The dOB attribute is therefore assigned with a negative integer, corresponding to the number of days elapsed
        between the birth and the model initialisation.
        """
        self.dOB = current_time - round(age*365.25)
        if age == 0.:
            self.dOB = current_time - round(random.uniform(low=0., high=time_step))

    def set_death_date(self, life_duration):
        """
        Defines the death date. This concerns natural death so this date might change in case of active TB
        """
        proposed_dOD = self.dOB + round(life_duration*365.25)
        if proposed_dOD > 0.:
            self.programmed['death'] = proposed_dOD
        else:
            self.programmed['death'] = random.randint(0, 10.*365.25, 1)[0]

    def set_date_leaving_home(self, params, minimum_date=None):
        """
        Define the date at which individuals will be leaving the household
        :param minimum_date: if specified, lower bound for the leaving date
        """
        self.programmed['leave_home'] = round(self.dOB + 365.25 *\
                                                   random.uniform(params['minimal_age_leave_hh'],
                                                                  params['maximal_age_leave_hh']))
        if minimum_date is not None:
            self.programmed['leave_home'] = max(self.programmed['leave_home'], minimum_date)

    def set_school_and_work_details(self, params):
        """
        Define the three following dates:
         - go_to_school_date
         - out_of_school_date
         - out_of_work_date
        """
        age_go_to_school = params['school_age'] * 365.25 + random.uniform(-1., 1.) * 365.25  # in days
        age_out_of_school = params['active_age_low'] * 365.25 + random.uniform(-3., 3.) * 365.25  # in days
        age_out_of_work = random.uniform(55., 70.) * 365.25  # in days

        self.programmed['go_to_school'] = self.dOB + round(age_go_to_school)
        self.programmed['leave_school'] = self.dOB + round(age_out_of_school)
        self.programmed['leave_work'] = self.dOB + round(age_out_of_work)

        active = random.binomial(1, params['perc_active']/100.)
        if active == 1:
            self.is_ever_gonna_work = True

    def get_age_in_years(self, time):
        return (time - self.dOB)/365.25

    def get_relative_susceptibility(self, time, params):
        """
        Return the relative susceptibility to infection of the individual. This quantity depends on the following factors:
         vaccination status, age ...
        Returns the relative susceptibility. Baseline is for a non-vaccinated individual.
        """
        rr = 1.
        if self.vaccinated:
            age = self.get_age_in_years(time)
            if age <= params['bcg_start_waning_year']:
                efficacy = params['bcg_maximal_efficacy']
            elif age <= params['bcg_end_waning_year']:
                efficacy = params['bcg_maximal_efficacy'] - \
                           (age - params['bcg_start_waning_year']) * \
                                                            params['bcg_maximal_efficacy'] / \
                                                            (params['bcg_end_waning_year'] - params['bcg_start_waning_year'])
            else:
                efficacy = 0.
            rr *= 1. - efficacy

        if self.ltbi:
            rr *= params['latent_protection_multiplier']
        return rr

    def get_relative_infectiousness(self, params, time):
        """
            Return the relative infectiousness of the individual. This quantity depends on the following factors:
            detection status, treatment status. Smear status
            Returns the relative infectiousness. Baseline is for an undetected Smear-positive TB case.
        """
        age = self.get_age_in_years(time)

        # sigmoidal scale-up
        rr = 1. / (1. + exp(-(age - params['infectiousness_switching_age'])))

        # alternate profile for infectiousness
        if params['linear_scaleup_infectiousness']:
            if age <= 10.:
                rr = 0.
            elif age >= 15:
                rr = 1.
            else:
                rr = 0.2 * age - 2.

        # organ-manifestation
        if self.tb_organ == '_smearneg':
            rr *= params['rel_infectiousness_smearneg']
        elif self.tb_organ == '_extrapulmonary':
            return 0.

        # detection status
        if 'detection' in self.programmed.keys():
            if time >= self.programmed['detection']:
                rr *= params['rel_infectiousness_after_detect']

        return rr

    def infect_individual(self, time, params, strain):
        """
        The individual gets infected with LTBI at time "time".
        """
        self.ltbi = True
        self.tb_strain = strain
        self.determine_activation(time, params)

    def determine_activation(self, time, params):
        """
        Determine whether and when the infected individual will activate TB.
        """
        time_to_activation = tb_activation.generate_an_activation_profile(self.get_age_in_years(time), params)

        #  if time_to_activation is not None:
        if round(time + time_to_activation) < self.programmed['death']:
            # the individual will activate TB
            self.programmed['activation'] = round(time + time_to_activation)

    def test_individual_for_ltbi(self, params):
        """
        Apply screening to an individual. This individual may or may not be infected.
        return: a boolean variable indicating the result of the test (True for positif test)
        """
        if self.ltbi:  # infected individual.
            test_result = bool(random.binomial(n=1, p=params['ltbi_test_sensitivity']))
        else:  # not infected
            if self.vaccinated:  # bcg affects specificity
                test_result = bool(random.binomial(n=1, p=1. - params['ltbi_test_specificity_if_bcg']))
            else:
                test_result = bool(random.binomial(n=1, p=1. - params['ltbi_test_specificity_no_bcg']))
        return test_result

    def get_preventive_treatment(self, params, time=0, delayed=False):
        """
        The individual receives preventive treatment. If the individual is currently infected, infection vanishes with
        probability pt_efficacy. The efficacy parameter represents a combined rate of adherence and treatment efficacy.
        Return the date of prevented activation, which is the date of TB activation in the event that the individual was
        meant to progress to active TB.
        """
        date_prevented_activation = None
        if self.ltbi:
            success = random.binomial(n=1, p=params['pt_efficacy'])
            if success == 1:
                self.ltbi = False
                if 'activation' in self.programmed.keys():
                    if not delayed or (delayed and (time + params['pt_delay_due_to_tst']) <= self.programmed['activation']):
                        date_prevented_activation = self.programmed['activation']
                        del self.programmed['activation']

        return date_prevented_activation

    def activate_tb(self):
        """
        Make the individual activate TB
        """
        self.active_tb = True
        self.ltbi = False  # convention

    def define_tb_outcome(self, time, params, tx_success_prop):
        """
        This method determines the outcome of the individual's active TB episode, accounting for both natural history
         and clinical management.
        :return:
        A dictionary containing the information to be recorded from the model perspective in case the programmed events
        dictionary needs to be updated. The keys of the returned dictionary may be "death", "recovery" and/or "detection".
        The values are the dates of the associated events.
        """
        # Smear-positive or Smear-negative
        organ_probas = [params['perc_smearpos']/100.,  # smear_pos
                        (100. - params['perc_smearpos'] - params['perc_extrapulmonary'])/100.,  # smear_neg
                        params['perc_extrapulmonary']/100.]  # extrapulmonary

        draw = random.multinomial(1, organ_probas)
        index = int(nonzero(draw)[0])
        self.tb_organ = ['_smearpos', '_smearneg', '_extrapulmonary'][index]

        # Natural history of TB
        if self.tb_organ == '_smearpos':
            organ_for_natural_history = '_smearpos'
        else:
            organ_for_natural_history = '_closed_tb'

        if params['new_tbnh_parameters']:
            t_to_sp_cure = round(365.25 * random.exponential(scale=1. / params['rate_sp_cure' +
                                                                               organ_for_natural_history]))
            t_to_tb_death = round(365.25 * random.exponential(scale=1. / params['rate_tb_mortality' +
                                                                                organ_for_natural_history]))
            if t_to_sp_cure <= t_to_tb_death:
                sp_cure = 1
                t_s = t_to_sp_cure
                t_m = float('inf')
            else:
                sp_cure = 0
                t_s = float('inf')
                t_m = t_to_tb_death
        else:
            sp_cure = random.binomial(n=1, p=params['perc_sp_cure' + organ_for_natural_history]/100.)
            if sp_cure == 1:
                t_s = round(365.25 * random.exponential(scale=1. / params['lambda_timeto_sp_or_death' +
                                                                          organ_for_natural_history]))
                t_m = float('inf')  # infinite
            else:
                t_s = float('inf')  # infinite
                t_m = round(365.25 * random.exponential(scale=1. / (params['lambda_timeto_sp_or_death' +
                                                                           organ_for_natural_history] *
                                                                    params['tb_mortality_multiplier'])))

        # Random generation of programmatic durations
        [t_d, t_t] = random.exponential(scale=[365.25/params['lambda_timeto_detection' + organ_for_natural_history],
                                               params['time_to_treatment']])
        t_d = round(t_d)
        t_t = round(t_t)

        # In case of spontaneous cure occurring before detection
        if sp_cure == 1 and t_d >= t_s:
            # In case spontaneous cure occurs before natural death
            if time + t_s < self.programmed['death']:
                self.programmed['recovery'] = time + t_s
                return {'recovery': time + t_s, 'time_active': t_s}
            else:  # natural death occurs before sp cure. Death will be registered as TB death
                self.die_with_tb = True
                return {'time_active': self.programmed['death'] - time}
        # In case of tb death before detection
        elif sp_cure == 0 and t_d >= t_m:
            self.die_with_tb = True
            if time + t_m < self.programmed['death']:
                return {'death': time + t_m, 'time_active': t_m}
            else:
                return {'time_active': self.programmed['death'] - time}
        # In case of detection effectively happening
        elif time + t_d < self.programmed['death']:
            self.programmed['detection'] = time + t_d
            to_be_returned = {'detection': self.programmed['detection'],
                              'time_active': self.programmed['death'] - time}
            # In case of sp cure occurring after detection and before natural death
            if sp_cure == 1 and time + t_s < self.programmed['death']:
                self.programmed['recovery'] = time + t_s
                to_be_returned['recovery'] = self.programmed['recovery']
                to_be_returned['time_active'] = t_s
            # In case of treatment effectively happening
            if time + t_d + t_t < self.programmed['death']:
                strain_multiplier = 1.
                if self.tb_strain == 'mdr':
                    strain_multiplier = params['perc_dst_coverage'] / 100.
                    strain_multiplier *= params['relative_treatment_success_rate_mdr']
                tx_cure = random.binomial(n=1, p=tx_success_prop * strain_multiplier)
                if tx_cure == 1:
                    if t_d + t_t < t_s: # will not overwrite the sp_cure date if it happens before treatment
                        self.programmed['recovery'] = time + t_d + t_t
                        to_be_returned['recovery'] = self.programmed['recovery']
                        to_be_returned['time_active'] = t_d + t_t
                elif tx_cure == 0 and self.tb_strain == 'ds':  # there is a risk of DR amplification
                    ampli = random.binomial(n=1, p=params['perc_risk_amplification'] / 100.)
                    if ampli == 1:
                        self.tb_strain = 'mdr'  # may be improved in the future as the amplification should occur later
                        to_be_returned['dr_amplification'] = time + t_d + t_t
            return to_be_returned
        # Otherwise, natural death occurs before detection
        else:
            return {'time_active': self.programmed['death'] - time}

    # def record_contacts(self, contact_dict):
    #     """
    #     update the attribute self.contacts_while_tb according to the content of contact_dict
    #     """
    #     for location in self.contacts_while_tb.keys():
    #         for contact_id in contact_dict[location].keys():
    #             if contact_id not in self.contacts_while_tb[location]:
    #                 self.contacts_while_tb[location].append(contact_id)

    def detect_tb(self):
        self.tb_detected = True

    def recover(self):
        """
        Make the individual recover from TB
        """
        self.active_tb = False
        self.ltbi = False
        self.tb_strain = None
        self.tb_organ = None
        self.tb_detected = False
        self.tb_treated = False
        self.die_with_tb = False
        self.contacts_while_tb = {'household': set([]), 'school': set([]), 'workplace': set([]), 'community': set([])}

        if 'activation' in self.programmed.keys():
            del(self.programmed['activation'])
        if 'recovery' in self.programmed.keys():
            del (self.programmed['recovery'])

    def reset_params(self):
        """
        Reset most params to default values. The dOB and dOD are initialised in another process.
        """
        self.diabetes = False
        self.hiv = False
        self.ltbi = False
        self.active_tb = False
        self.vaccinated = False
        self.group_ids = []

        self.contacts_while_tb = {'household': set([]), 'school': ([]), 'workplace': ([]), 'community': ([])}

        self.tb_organ = None

        self.tb_detected = False
        self.tb_treated = False
        self.die_with_tb = False

        self.is_ever_gonna_work = False

        self.programmed = {}

    def assign_vaccination_status(self, coverage):
        """
        When an individual is born, vaccination occurs with probability 'coverage'
        """
        vacc = random.binomial(1, coverage)
        if vacc == 1:
            self.vaccinated = True

    def is_in_subgroup(self, subgroup, time=None):
        """
        Binary test to inform whether the individual belongs to a specific subgroup.
        subgroup is one of "children" (0-15 yo), "young_children" (0-5 yo)), "adult" (15yo+)
        return: boolean variable
        """
        test = False

        if subgroup == 'children':
            age = self.get_age_in_years(time)
            if age <= 15.:
                test = True
        elif subgroup == 'young_children':
            age = self.get_age_in_years(time)
            if age <= 5.:
                test = True
        elif subgroup == 'adult':
            age = self.get_age_in_years(time)
            if age > 15.:
                test = True
        else:
            print "subgroup is not admissible"

        return test
