from sys import exit, stdout
import agent
import household
import toolkit
import numpy as np
import copy
from math import ceil, floor
import os
from itertools import repeat
from datetime import datetime
from importData import get_agecategory
import dill
import time
from calibration_targets import calib_targets

def age_preference_function(age_difference, sigma):
    """
    Given the age difference between two individuals, computes the relative probability of contact.
    Reference: no age difference.
    location_type is either "school" or "workplace"
    """
    relative_risk = np.exp(-age_difference ** 2 / (2 * sigma ** 2))
    return relative_risk

class Model:
    """
    Defines the whole population and interactions between individuals
    Also account for the household structure
    """
    def __init__(self, data, i_seed, scenario, i_run, initialised=True):
        self.i_seed = i_seed
        self.scenario = scenario
        self.i_run = i_run
        self.timer = time.time()
        self.params = {}
        self.age_pyramid = data.age_pyramid
        self.fertility_replacement = True  # default model behaviour. May switch during simulation
        self.constant_birth_rate = False  # may be switch on after burn-in if age-pyramid required.
        self.stopped_simulation = False  # in case simulation has been forced to stop
        self.stopped_time = None  # when stopped
        self.status_file_created = False
        self.sigmoidal_birth_rate_function = None
        self.scale_up_functions = data.scale_up_functions
        self.scale_up_functions_current_time = {}  # store the output of scale_up_functions at each time-step
        self.remaining_calibration_targets = {}  # keys are years and values are dictionaries with targets

        self.contact_rates_matrices = data.contact_rates_matrices
        self.prem_contact_rate_functions = data.prem_contact_rate_functions
        self.sd_agepref_work = data.sd_agepref_work
        self.individuals = {} # dictionary of all individuals keyed by their unique ID
        self.households = {} # dictionary of all households keyed by their unique ID
        self.last_ind_id = 0
        self.dates_of_birth = {}
        self.ind_by_agegroup = {}
        self.initialise_agegroups()
        self.age_pyramid_date = 0
        self.birth_numbers = 0  # reset at each step

        self.birth_numbers_function = data.birth_numbers_function

        self.last_ind_id = -1
        self.population = 0
        self.n_households = 0
        self.greatest_ind_id = 0
        self.eligible_hh_for_birth = {}  # keys: hh ids, values: hh sizes
        self.individuals_want_to_move = {}  # dictionary of ind ids as keys and desired moving times at values
        self.empty_households = []  # list of empty households

        self.time = int(0)
        self.last_year_completed = 0
        self.time_reset_records = float('inf')  # time at which recording attributes are reset
        self.records_have_been_reset = False
        self.has_been_stored = False

        self.mean_age = 0.  # average age of the population
        self.prop_under_5 = 0.

        self.contact_matrices = {}
        self.n_contacts = {}
        self.initialise_contact_matrices()

        self.pool_of_life_durations = None

        self.groups = {}  # keyed by group ids. valued with the ids of individuals that belong to each group.
        self.groups_by_type = {'schools': [], 'workplaces': []}  #  list the group ids by group type
        self.group_types = {}  #  reverse version of groups_by_type. Keys are group ids and values are group types

        self.programmed_events = {'death': {}, 'leave_home': {}, 'go_to_school': {}, 'leave_school': {}, 'leave_work': {}}

        # storage units for model outputs
        self.timeseries_log = {
            'times': []
        }  # to store time_series type output variables
        self.checkpoint_outcomes = {'ages': {}, 'household_sizes': {},
                                    'school_sizes': {},
                                    'workplace_sizes': {},
                                    'prevalence_by_age': {}}  # contain output measures that are not time-series

        self.activation_stats = {'n_infections': 0, 'n_activations': 0}
        self.ltbi_age_stats = {'ltbi_ages': [], 'ending_tb_ages': []}
        self.ltbi_age_stats_have_been_recorded = False
        self.prevalence_by_age_recorded = False

    #  ________________________________________________________________
    #               Model initialisation
    def initialise_contact_matrices(self):
        """
        This method initialises the contact matrices storage
        """
        for key in ['contact', 'transmission', 'transmission_end_tb']:
            self.contact_matrices[key] = {}
            self.n_contacts[key] = {}
            for location in ['school', 'workplace', 'household', 'community']:
                self.contact_matrices[key][location] = np.zeros((101, 101))  # null matrix 100x100
                self.n_contacts[key][location] = 0

    def print_run_characteristics(self):
        for organ in ['_smearpos', '_closed_tb']:
            print "The average time to detection for \"" + organ + "\" is " + \
                  str(round(1. / self.params['lambda_timeto_detection' + organ],2)) + " years"

    def evaluate_all_scale_up_functions(self):
        for scale_up_key, func in self.scale_up_functions.iteritems():
            remaining_years = (float(self.age_pyramid_date) - self.time) / 365.25
            time_in_year_date = self.params['current_year'] - remaining_years
            if self.params['time_variant_programmatic']:
                self.scale_up_functions_current_time[scale_up_key] = func(time_in_year_date)
            else:
                self.scale_up_functions_current_time[scale_up_key] = func(2050.)

        self.process_cdr()
        self.process_organ_proportions()

    def initialise_model(self, data):
        self.collect_params(data)
        self.population = self.params['population']
        self.evaluate_all_scale_up_functions()

        if self.scenario != 'init':
            self.collect_scenario_specific_params(data)
        self.initialise_timeseries_storage()
        self.generate_households()
        self.populate_households()
        self.build_schools_and_workplaces()

        self.tb_has_started = False
        self.run()

    """
            Methods related to model initialisation (parameter processing + storage initialisation)
    """
    def initialise_timeseries_storage(self):
        """
        read the list of timeseries to record and initialise their storage accordingly
        """
        for ts_name in self.params['timeseries_to_record']:
            self.timeseries_log[ts_name] = []

    def reset_params(self, data):
        """
        Re-initialise the parameters dictionary.
        """
        self.params = {}
        self.collect_params(data)

    def collect_params(self, data):
        """
        collect the params stored in the data attribute of the model_runner object
        and populates the params attribute of the model object
        """
        # collect common parameters
        for key, value in data.console.iteritems():
            self.params[key] = value
        for key, value in data.common_parameters.iteritems():
            self.params[key] = value
        self.prem_contact_rate_functions = data.prem_contact_rate_functions
        self.params['age_pyramid'] = data.age_pyramid
        self.params['activation_times_dic'] = data.activation_times_dic
        self.process_n_iterations()
        self.process_calibration_targets()
        # self.process_cdr()
        self.process_age_pyramid_date()
        self.process_reset_records_date()
        self.process_rr_transmission_by_location()
        self.pool_of_life_durations = data.pool_of_life_durations

    def collect_scenario_specific_params(self, data):
        # collect or rewrite scenario-specific parameters
        if self.scenario != 'init':
            for key, value in data.scenarios[self.scenario].iteritems():
                self.params[key] = value
        if 'time_step' in data.scenarios[self.scenario].keys() or 'n_years' in data.scenarios[self.scenario].keys():
            self.process_n_iterations()

    def process_n_iterations(self):
        """
        Calculate the number of iterations needed
        """
        if self.initialised:
            n_years = self.params['n_years']
        else:
            n_years = self.params['duration_burning_demo'] + self.params['duration_burning_tb']

        self.params['n_iterations'] = int(ceil(n_years * 365.25 / self.params['time_step']))
        self.params['checkpoints'] = [(cp*365.25 - cp*365.25 % self.params['time_step']) for cp in self.params['years_checkpoints']]

    def process_age_pyramid_date(self):
        if self.params['running_mode'] != 'manual':
            self.age_pyramid_date = self.params['duration_burning_demo'] + self.params['duration_burning_tb'] +\
                                    self.params['n_years']
        else:
            self.age_pyramid_date = self.params['duration_burning_demo'] + self.params['duration_burning_tb']
        self.age_pyramid_date = int(self.age_pyramid_date * 365.25)

    def process_reset_records_date(self):
        if self.params['reset_records_time'] == 'after_demo_burning':
            self.time_reset_records = int(365.25*self.params['duration_burning_demo'])
        elif self.params['reset_records_time'] == 'after_tb_burning':
            self.time_reset_records = int(365.25*(self.params['duration_burning_demo'] +\
                                                  self.params['duration_burning_tb']))

    def process_calibration_targets(self):
        self.remaining_calibration_targets = {}
        if self.params['country'] in calib_targets.keys():
            for target in calib_targets[self.params['country']]:
                if target['year'] not in self.remaining_calibration_targets.keys():
                    self.remaining_calibration_targets[target['year']] = []
                self.remaining_calibration_targets[target['year']].append(target)
  
    def process_rr_transmission_by_location(self):
        self.params['rr_transmission_by_location'] = {}
        for location in self.contact_matrices['contact'].keys():
            prop_physical = self.params['prop_physical_' + location + '_contact']
            self.params['rr_transmission_by_location'][location] =\
                prop_physical + (1 - prop_physical) * self.params['rr_transmission_nonphysical_contact']

        ref_rr = self.params['rr_transmission_by_location']['household']
        for key, val in self.params['rr_transmission_by_location'].iteritems():
            self.params['rr_transmission_by_location'][key] = val/ref_rr

    def initialise_agegroups(self):
        self.ind_by_agegroup = {}
        for cat_ind in range(16):
            age_cat = 'X_' + str(cat_ind+1)
            self.ind_by_agegroup[age_cat] = []

    def update_ind_by_agegroup(self):
        self.initialise_agegroups()
        for key, value in self.dates_of_birth.iteritems():
            age = (self.time - value) / 365.25
            self.ind_by_agegroup[get_agecategory(age,'prem')].append(key)

    """
             Methods related to population initialisation (households and individuals + generation)
    """
    def generate_households(self):
        """
        Generate households
        """
        nb_hh = int(round(self.params['population'] / self.params['household_size']))
        for i in range(nb_hh - 1):
            h_id = i
            self.households[h_id] = household.household(h_id)
            self.n_households += 1
        # Create an extra household which is a retirement home
        self.households[h_id + 1] = household.household(h_id + 1, retirement_home=True)
        self.n_households += 1

    def populate_households(self):
        """
        Populate existing households
        """
        # Allocate the adults first
        parenting_households = [] # store the hh ids where kids can be allocated during the initialisation phase
        for h in self.households.values():
            if h.retirement_home:
                continue
            parents_age = np.random.uniform(self.params['minimal_age_leave_hh'], 100.)  # btwn 18 and 100.
            h.repopulate_date = self.time - 365.25*(parents_age - self.params['minimal_age_leave_hh'])
            if self.time - h.repopulate_date < 365.25*20.:
                self.eligible_hh_for_birth[h.id] = 2
            if parents_age < 60.:  # we generate a couple and we allow kids to join
                for j in [1, 2]:
                    self.add_new_individual_in_hh(h_id=h.id, age=parents_age)
                parenting_households.append(h.id)
            else:  # we generate a couple
                self.add_new_individual_in_hh(h_id=h.id, age=parents_age)
                self.add_new_individual_in_hh(h_id=h.id, age=parents_age)

        # Allocate the remaining individuals as kids
        n_kids = self.population - (self.last_ind_id + 1)
        for i in range(n_kids):
            hh_id = parenting_households[i % len(parenting_households)]
            kids_age = np.random.uniform(0., 40.)
            self.add_new_individual_in_hh(h_id=hh_id, age=kids_age)

    def add_new_individual_in_hh(self, h_id, age, ind_id=None):
        if ind_id is None:
            ind_id = self.last_ind_id + 1
            self.last_ind_id += 1

        self.individuals[ind_id] = agent.individual(id=ind_id, household_id=h_id, dOB=0.)

        self.set_birth_and_death(ind_id, age=age)
        if age == 0.:
            age_cat = 'X_1'
        else:
            age_cat = get_agecategory(age, 'prem')
        self.ind_by_agegroup[age_cat].append(ind_id)

        self.households[h_id].individual_ids.append(ind_id)
        self.households[h_id].size += 1
        self.households[h_id].last_baby_time = self.time

    def build_schools_and_workplaces(self):
        """
        Create the different groups (schools and workplaces).
        :return:
        """
        # Build schools and assign the different households to the schools
        n_schools = int(ceil(self.params['n_schools'] * self.population / 1.e5))
        self.groups_by_type['schools'] = range(1, n_schools + 1)
        for hh_id in self.households.keys():
            if self.households[hh_id].retirement_home:
                continue
            self.households[hh_id].school_id = np.random.choice(self.groups_by_type['schools'], 1)[0]
        for school_id in self.groups_by_type['schools']:
            self.groups[school_id] = []
            self.group_types[school_id] = 'school'

        # Build workplaces
        prop_in_working_age = sum([self.age_pyramid['X_' + str(i)] for i in range(5, 13)])
        n_working_individuals = self.population * prop_in_working_age * self.params['perc_active']/100.
        n_workplaces = int(ceil(n_working_individuals / (self.params['n_colleagues'] + 1.)))
        self.groups_by_type['workplaces'] = range(n_schools + 1, n_schools + n_workplaces + 1)
        for workplace_id in self.groups_by_type['workplaces']:
            self.groups[workplace_id] = []
            self.group_types[workplace_id] = 'workplace'

    def adjust_nb_of_schools_and_workplaces(self):
        """
        When population size changes, we may need to readjust the number of schools and workplaces
        """
        for group_type in ['school', 'workplace']:
            if group_type == 'school':
                ideal_nb_of_groups = int(ceil(self.params['n_schools'] * self.population / 1.e5))
            else:
                prop_in_working_age = sum([self.age_pyramid['X_' + str(i)] for i in range(5, 13)])
                n_working_individuals = self.population * prop_in_working_age * self.params['perc_active'] / 100.
                ideal_nb_of_groups = int(ceil(n_working_individuals / (self.params['n_colleagues'] + 1.)))
            if len(self.groups_by_type[group_type + 's']) > ideal_nb_of_groups:
                self.remove_groups(group_type=group_type,
                                   n_to_remove=(len(self.groups_by_type[group_type + 's']) - ideal_nb_of_groups))
            elif len(self.groups_by_type[group_type + 's']) < ideal_nb_of_groups:
                self.create_groups(group_type=group_type,
                                   n_to_add=(ideal_nb_of_groups - len(self.groups_by_type[group_type + 's'])))

    def create_groups(self, group_type, n_to_add):
        for _ in range(n_to_add):
            # create the group and add to relevant model attributes
            new_group_id = max(self.groups.keys()) + 1
            self.groups[new_group_id] = []  # create an empty group
            self.groups_by_type[group_type + 's'].append(new_group_id)
            self.group_types[new_group_id] = group_type

            if group_type == 'school':  # need to assign the new school to some households
                # nb_assigned_hh = np.random.poisson(lam=self.params['n_households_per_' + group_type])
                nb_assigned_hh = np.random.poisson(lam=float(len(self.households))/float(len(self.groups_by_type['schools'])))

                reassigned_hh = np.random.choice(self.households.keys(), nb_assigned_hh, replace=False)
                for hh_id in reassigned_hh:
                    self.households[hh_id].school_id = new_group_id

    def populate_new_group(self, new_group_id, group_type, reassigned_hh=None):
        """
        This method will assign some individuals to the new group. Individuals are randomly picked among those who are
         already part of another group of the same type.
         If the group is a school, we just pick the children from the households that have been associated with the new school.
         reassigned_hh contains the list of households that have just been assigned to the school.
        """
        if group_type == 'school':  # we take the kids from the relevant households
            for hh_id in reassigned_hh:
                for ind_id in self.households[hh_id].individual_ids:
                    # attend a group?
                    if len(self.individuals[ind_id].group_ids) > 0:
                        prev_group_id = self.individuals[ind_id].group_ids[0]
                        if self.group_types[prev_group_id] == 'school':  # is the group a school?
                            self.make_individual_change_group(ind_id, prev_group_id, new_group_id)
        else:  # this is a workplace. We pick workers randomly
            nb_workers_to_move = np.random.poisson(lam=self.params['n_colleagues'] + 1.)
            picked_workplaces = np.random.choice(self.groups_by_type['workplaces'],
                                                 min(nb_workers_to_move, len(self.groups_by_type['workplaces'])),
                                                 replace=False)
            for prev_workplace_id in picked_workplaces:
                if len(self.groups[prev_workplace_id]) > 0:
                    ind_id = np.random.choice(self.groups[prev_workplace_id], 1)[0]
                    self.make_individual_change_group(ind_id, prev_workplace_id, new_group_id)

    def make_individual_change_group(self, ind_id, prev_group_id, new_group_id):
        """
        Individual ind_id is moving from group prev_group_id to new_group_id. We update all relevant attributes
        """
        # from the individual's perspective
        self.individuals[ind_id].group_ids = [new_group_id]

        # from the model perspective
        self.groups[prev_group_id] = [indd_id  for indd_id in self.groups[prev_group_id] if indd_id != ind_id]
        self.groups[new_group_id].append(ind_id)

    def remove_groups(self, group_type, n_to_remove):
        closing_group_ids = np.random.choice(self.groups_by_type[group_type + 's'], n_to_remove, replace=False)
        for prev_group_id in closing_group_ids:
            # remove group_id from groups_by_type and group_types
            self.groups_by_type[group_type + 's'] = [g_id for g_id in self.groups_by_type[group_type + 's'] if
                                                     g_id != prev_group_id]
            del(self.group_types[prev_group_id])
            ind_ids_to_move = self.groups[prev_group_id]
            if group_type == 'school':  # we need to assign relevant households to new schools
                self.reassign_schools_after_school_closure(prev_group_id)
            for ind_id in ind_ids_to_move:
                if group_type == 'workplace':  # new workplace chosen at random
                    new_group_id = np.random.choice(self.groups_by_type[group_type + 's'], 1)[0]
                else:  # this is a school. A new school has been assigned to the individual's household
                    new_group_id = self.households[self.individuals[ind_id].household_id].school_id
                self.make_individual_change_group(ind_id, prev_group_id, new_group_id)

            del(self.groups[prev_group_id])

    def reassign_schools_after_school_closure(self, closing_school_id):
        """
        closing_school_id was just removed. We need to assign another school to the households that were associated to
        the closing school.
        Note that the attribute self.groups_by_type['schools'] has already been updated so closing_school_id is no longer
        in this list.
        """
        for h in self.households.values():
            if h.school_id == closing_school_id:
                h.school_id = np.random.choice(self.groups_by_type['schools'], 1)[0]

    def set_initial_tb_states(self):

        if self.params['init_n_tb_cases'] > 0:
            diseased_indices = np.random.choice(self.individuals.keys(), self.params['init_n_tb_cases'], replace=False)
            for ind in diseased_indices:
                tb_strain = ['ds', 'mdr'][int(np.random.binomial(n=1, p=self.params['init_mdr_perc'] / 100.))]
                self.individuals[ind].tb_strain = tb_strain
                self.make_individual_activate_tb(ind, init=True)
                self.individuals[ind].programmed['activation'] = self.time

            # TB activations have triggered decrements in ltbi prevalence. We do not want this in the initialisation phase.
            # self.ltbi_prevalence += self.params['init_n_tb_cases']

    def set_initial_infection_states(self):
        """
        Assign an infection status (infected or susceptible) to each individual of the initial population
        """
        n_ltbi = round(self.population * self.params['init_ltbi_prev']/100.)
        if n_ltbi > 0:
            infected_indices = np.random.choice(self.individuals.keys(), int(n_ltbi), replace=False)
            for ind in infected_indices:
                if not self.individuals[ind].active_tb:
                    tb_strain = ['ds', 'mdr'][int(np.random.binomial(n=1, p=self.params['init_mdr_perc']/100.))]
                    self.infect_an_individual(ind, strain=tb_strain)
                    self.ltbi_prevalence += 1

    def get_other_household_members(self, ind_id):
        """
        Provides the list of ids of the individuals living with "ind_id".
        return: a list of ids
        """
        hh = self.households[self.individuals[ind_id].household_id]
        return [ind for ind in hh.individual_ids if ind != ind_id]

    def list_eligible_households_for_birth(self):
        """
        This method lists all the households that are eligible to receive a newborn individual. The object attribute
         self.eligible_hh_for_birth is updated
        :return: nothing
        """
        self.eligible_hh_for_birth = {}
        for h in self.households.values():
            if not h.retirement_home and \
                            (self.time - h.repopulate_date) < 365.25 * self.params['duration_hh_eligible_for_birth'] and \
                    h.size > 0 and (self.time - h.last_baby_time) > h.minimum_time_to_next_baby:
                self.eligible_hh_for_birth[h.id] = h.size

    def pick_eligible_household_for_birth(self):
        """
        This method will randomly pick an eligible household among the ones listed in self.eligible_hh_for_birth.
        The hosuehold of small size will be favored
        :return: a household id
        """
        if len(self.eligible_hh_for_birth) > 0:
            hh_ids = self.eligible_hh_for_birth.keys()
            hh_sizes = self.eligible_hh_for_birth.values()

            probas = [1./float(x) for x in hh_sizes]
            probas = [x / sum(probas) for x in probas]

            draw = np.random.multinomial(1, probas)
            index = int(np.nonzero(draw)[0])
            hh_id = hh_ids[index]
        else:
            hh_size = 0
            while hh_size == 0:
                hh_id = np.random.choice(self.households.keys(), 1)[0]
                if not self.households[hh_id].retirement_home:
                    hh_size = self.households[hh_id].size
        return hh_id

    def take_elderly_to_retirement_home(self):
        """
        Each year, some of the individuals who live alone and who are >70 years old are going to the retirement home.
        The proba of going to retirement home is calculated so that 'perc_single_over90yo_in_retirement_home' % of single
        individuals >90 yo are in a retirement home.
        """
        proba_go_to_retirement_home = 1 - pow(1. - self.params['perc_single_over90yo_in_retirement_home']/100., 1./20.)
        for h in self.households.values():
            if not h.retirement_home and h.size == 1:
                ind_id = h.individual_ids[0]
                age = self.individuals[ind_id].get_age_in_years(self.time)
                if age >= 70.:
                    if bool(np.random.binomial(1, proba_go_to_retirement_home)):
                        # update individual's attributes
                        self.individuals[ind_id].household_id = self.n_households - 1
                        # update previous hh attributes
                        self.households[h.id].size = 0
                        self.households[h.id].individual_ids = []
                        self.empty_households.append(h.id)

                        # update retirement home attributes
                        self.households[self.n_households - 1].size += 1
                        self.households[self.n_households - 1].individual_ids.append(ind_id)

    def lighten_want_to_move_home(self):
        """
        If the queue self.individuals_want_to_move becomes too long, we make some individuals stay home but
        allow births to happen in the original home
        """
        # If an individual has been waiting for a new home for more that 1 year, we allow birth to happen in his/her
        # household
        updates_ind_ids = []
        for ind_id, time in self.individuals_want_to_move.iteritems():
            if self.time - time > 365.25:
                h = self.households[self.individuals[ind_id].household_id]
                if h.size > 0:
                    # Allow birth to happen
                    h.repopulate_date = self.time
                    self.eligible_hh_for_birth[h.id] = h.size
                updates_ind_ids.append(ind_id)

        # clean self.individuals_want_to_move dictionary
        for ind_id in updates_ind_ids:
            del(self.individuals_want_to_move[ind_id])

    def run(self):
        """
        run the initialised model for n_iterations time-steps (weeks)
        """
        self.status_file_created = False
        # If the model is already initialised, we need to update the number of iterations
        if self.initialised:
            self.process_n_iterations()

        if self.params['force_tb_init']:
            self.tb_has_started = False  # the tb initialisation process will happen in any case

        for i in range(self.params['n_iterations']):

            if self.time >= self.time_reset_records and not self.records_have_been_reset:
                self.reset_recording_attributes()
                self.records_have_been_reset = True
            self.move_forward()
            if self.params['stop_if_condition'] and not self.stopped_simulation:
                stop = self.shall_we_stop()
                if stop:
                    self.stop_running_model()

        if not self.ltbi_age_stats_have_been_recorded:
            self.record_ltbi_ages()

        # post-simulation treatment (does not apply to initialisation runs)
        if self.initialised and not self.stopped_simulation:
            # Contact heatmaps
            if self.params['plot_contact_heatmap']:
                self.record_all_contacts()  # contacts occurring over the last time-step are recorded

            # Moving average for some timeseries
            self.clean_timeseries()
            self.average_timeseries()

            if self.params['running_mode'] in ['find_a_calibrated_model', 'run_ks_based_calibration',
                                               'run_lhs_calibration']:
                self.store_calibrated_model()
            # average time active TB
            # print "Average time a TB case is active: " + str(round(self.time_active['total_time_active'] / self.time_active['total_n_cases'], 2)) + " days"

            if self.activation_stats['n_infections'] > 0:
                perc_activation = 100.* self.activation_stats['n_activations'] / self.activation_stats['n_infections']
                print "The percentage of activations among all infections is " + str(round(perc_activation, 2)) + "%"

        # change status file name when simulation completed
        if self.status_file_created:
            dir_path = os.path.join('outputs', self.params['project_name'])
            file_path = os.path.join(dir_path, 'progress_seed' + str(self.i_seed) + '_' + self.scenario + '_run' + str(self.i_run) + '.txt')
            if self.params['running_mode'] == 'manual' and self.params['load_calibrated_models']:
                os.remove(file_path)
            else:
                new_file_path = os.path.join(dir_path, 'complete_seed' + str(self.i_seed) + '_' + self.scenario + '_run' + str(self.i_run) + '.txt')
                os.rename(file_path, new_file_path)

    def move_forward(self):
        if self.params['run_universal_methods']:
            self.run_universal_methods()  # what needs to be done for every single individual at every step
        self.store_variables()
        self.reset_attributes()
        self.time += self.params['time_step']
        self.evaluate_all_scale_up_functions()

        if not self.stopped_simulation:
            # the below commands are not run if the simulation has been forced to stop

            # demo only
            if self.time >= (self.age_pyramid_date - 90.*365.25):  # we may want to target the age-pyramid
                if self.params['birth_process'] == 'agepyramid_driven':
                    self.fertility_replacement = False

            if self.time >= self.age_pyramid_date:
                if self.params['birth_process'] == 'agepyramid_driven':
                    if self.sigmoidal_birth_rate_function is None: # create the sigmoidal birth-rate function
                        # self.params['fixed_birth_rate_based_on'] = 'latest_birth_rate'  # prov
                        latest_birth_rate = self.timeseries_log['birth_rate'][-1]
                        if self.params['fixed_birth_rate_based_on'] == 'mortality_rate':
                            mu = 1 / self.params['life_expectancy']  # per year
                            future_birth_rate = mu * 1000.
                        elif self.params['fixed_birth_rate_based_on'] == 'latest_birth_rate':
                            # average birth rates over the last 10 years
                            n_iter_to_consider = int(floor(365.25 * 10. / self.params['time_step']))
                            future_birth_rate = np.mean(self.timeseries_log['birth_rate'][-n_iter_to_consider:])
                        else:
                            future_birth_rate = self.params['birth_rate']
                        self.sigmoidal_birth_rate_function = toolkit.make_sigmoidal_curve(
                            y_low=latest_birth_rate, y_high=future_birth_rate, x_start=self.time,
                            x_inflect=self.time+365.25*5., multiplier=1.)
                    self.constant_birth_rate = True
                    self.params['birth_rate'] = self.sigmoidal_birth_rate_function(self.time)
                    if not self.ltbi_age_stats_have_been_recorded:
                        self.record_ltbi_ages()
                        self.ltbi_age_stats_have_been_recorded = True
                else:
                    self.fertility_replacement = True

            if self.get_current_date() >= self.params['prevalence_by_age_record_time'] and not self.prevalence_by_age_recorded:
                self.record_tb_prevalence_by_age()
                print "Prevalence by age recorded"
                self.prevalence_by_age_recorded = True

            if self.time > (float(self.last_year_completed + 1) * 365.25):  # one year has been completed
                self.run_yearly_methods()

            # TB below
            if self.time >= self.params['duration_burning_demo'] * 365.25:
                # turn transmission on if requested
                self.transmission = self.params['transmission']
                if not self.tb_has_started: # TB is starting right now then
                    if self.params['init_n_tb_cases'] > 0 or self.params['init_ltbi_prev'] > 0:
                        if not self.initialised and self.params['duration_burning_tb'] == 0:
                            pass
                        else:
                            self.set_initial_tb_states()
                            self.set_initial_infection_states()
                            self.tb_has_started = True

            if self.tb_has_started:
                if self.time >= (self.params['duration_burning_demo'] + self.params['duration_burning_tb'] + self.params['intervention_start_delay'] ) * 365.25:
                    self.apply_interventions()

                self.trigger_programmed_activations()
                if self.transmission:
                    self.spread_infections()
                self.trigger_programmed_detections()
                self.trigger_programmed_recoveries()
                self.trigger_programmed_tb_deaths()

            if self.time > (float(self.last_year_completed + 1) * 365.25):  # one year has been completed
                self.last_year_completed += 1
                # major simulation stop?
                file_path = os.path.join('outputs', self.params['project_name'], 'keep_running.txt')
                if not os.path.exists(file_path):
                    exit('Process exit from model.py: "keep_running" file was deleted')
                if self.params['print_time']:
                    print self.scenario + ' run ' + str(self.i_run) + ': year ' + str(self.last_year_completed) + \
                          ' (' + str(int(self.get_current_date())) + ') completed. Absolute tb_prevalence: ' + \
                          str(int(self.tb_prevalence)) + ' (' + \
                          str(int(round(1.e5 * self.tb_prevalence / self.population))) + ' /100k)'
                    if os.name != 'nt':  # 'nt' for windows system
                        self.write_status_file()

            # demo only again
            self.trigger_programmed_deaths()

            if not self.fertility_replacement:
                self.trigger_births()

            self.update_want_to_move_list()
            self.trigger_individuals_move_home()
        # self.check_workplace_ages() # debugging
        # self.check_school_ages() # debugging

    def stop_running_model(self):
        print "!!!!!!!!!!!!!!!!!!!!      Stopping model run " + str(self.i_run) + " of scenario " + self.scenario
        print "!!!!!!!!!!!!!!!!!!!!      at time: " + str(self.time)
        # record contacts now if requested
        if self.params['plot_contact_heatmap']:
            self.record_all_contacts()  # contacts occurring over the last time-step are recorded
        self.stopped_simulation = True
        self.stopped_time = self.time

    def shall_we_stop(self):
        """
        Will make simulation stop if some conditions are verified
        """
        stop = False
        if self.timeseries_log['tb_prevalence'][-1] >= self.params['prevalence_max']:
            stop = True
            print "Model run will be forced to stop because tb prevalence is too high."
        if self.tb_prevalence == 0 and len(self.programmed_events['activation'].values()) == 0 and self.time > 365.25*(
            self.params['duration_burning_demo'] + self.params['duration_burning_tb'] + 1.) and self.params['transmission']:
            stop = True
            print "Model run will be forced to stop because there is no more TB."
        return stop

    def check_individuals(self):
        pass

    def run_universal_methods(self):
        """
        Measures and/or actions that have to be performed for every single individual in the system.
        Grouping these methods in a same loop allows us to save computing time.
        """
        # initialisation of the different variables
        self.mean_age = 0.
        self.prop_under_5 = 0.

        # loop through all individuals
        for ind_id, ind in self.individuals.iteritems():
            age = ind.get_age_in_years(self.time)
            self.mean_age += age
            if age < 5.:
                self.prop_under_5 += 1.

        # post-loop treatment of the variables
        self.mean_age /= self.population
        self.prop_under_5 /= self.population

    def run_yearly_methods(self):
        """
        Run methods that need to be run once a year
        """
        self.trigger_programmed_go_to_school()
        self.trigger_programmed_leave_school()
        self.trigger_programmed_leave_work()

        if not self.fertility_replacement:  # if fertility_replacement, population size is constant
            self.build_or_destroy_households()
        self.list_eligible_households_for_birth()

        self.take_elderly_to_retirement_home()
        self.lighten_want_to_move_home()

        self.update_ind_by_agegroup()

    def store_variables(self):
        """
        Populate dictionaries that store the model outputs
        """
        self.timeseries_log['times'].append(self.time)

        # Population
        if 'population' in self.params['timeseries_to_record']:
            self.timeseries_log['population'].append(self.population)

        # Birth-rates
        if 'birth_rate' in self.params['timeseries_to_record']:
            self.timeseries_log['birth_rate'].append(
                365.25 * 1000. * self.birth_numbers / (self.population * self.params['time_step']))

        # LTBI prevalence
        if 'ltbi_prevalence' in self.params['timeseries_to_record']:
            self.timeseries_log['ltbi_prevalence'].append(100.*self.ltbi_prevalence/self.population)

        # TB incidence
        if 'tb_incidence' in self.params['timeseries_to_record']:
            self.timeseries_log['tb_incidence'].append(
                365.25 * 1.e5 * self.tb_incidence / (self.population * self.params['time_step']))

        # TB prevalence
        if 'tb_prevalence' in self.params['timeseries_to_record']:
            self.timeseries_log['tb_prevalence'].append(1.e5*self.tb_prevalence/self.population)

        if 'tb_prevalence_ds' in self.params['timeseries_to_record']:
            self.timeseries_log['tb_prevalence_ds'].append(1.e5 * self.tb_prevalence_ds / self.population)

        if 'tb_prevalence_mdr' in self.params['timeseries_to_record']:
            self.timeseries_log['tb_prevalence_mdr'].append(1.e5 * self.tb_prevalence_mdr / self.population)

        if 'prop_mdr_prevalence' in self.params['timeseries_to_record']:
            prop_mdr = 0. if self.tb_prevalence == 0. else self.tb_prevalence_mdr / self.tb_prevalence
            self.timeseries_log['prop_mdr_prevalence'].append(prop_mdr)

        # TB deaths
        if 'tb_deaths' in self.params['timeseries_to_record']:
            self.timeseries_log['tb_deaths'].append(
                365.25 * 1.e5 * self.tb_deaths / (self.population * self.params['time_step']))

        # preventive treatments provided
        if 'n_pt_provided' in self.params['timeseries_to_record']:
            self.timeseries_log['n_pt_provided'].append(self.n_pt_provided)
        if 'n_useful_pt_provided' in self.params['timeseries_to_record']:
            self.timeseries_log['n_useful_pt_provided'].append(self.n_useful_pt_provided)

        if self.params['run_universal_methods']:
            # mean age
            if self.params['plot_ts_mean_age']:
                self.timeseries_log['mean_age'].append(self.mean_age)
            # proportion of people under 5yo
            if self.params['plot_ts_prop_under_5']:
                self.timeseries_log['prop_under_5'].append(self.prop_under_5)

        if self.params['running_mode'] in ['find_a_calibrated_model', 'run_ks_based_calibration', 'run_lhs_calibration']:
            self.check_calibration_targets()

        # checkpoints data
        if self.time in self.params['checkpoints']:
            self.generate_checkpoint_outputs()

    def reset_attributes(self):
        """
        Reset attributes that have to be re-initialised at each time-step.
        """
        self.birth_numbers = 0
        self.tb_incidence = 0.
        self.tb_deaths = 0.
        self.n_pt_provided = 0.
        self.n_useful_pt_provided = 0.

    def set_birth_and_death(self, ind_id, age):
        """
        update the dOB and dOD attributes of the individual ind_id, given a current age
        """
        self.individuals[ind_id].set_dOB(age, self.time, self.params['time_step'])
        self.individuals[ind_id].set_death_date(np.random.choice(self.pool_of_life_durations, 1)[0])
        self.add_event_to_programmed_events('death', ind_id)
        if self.time > 0:
            self.individuals[ind_id].set_date_leaving_home(params=self.params)
            self.add_event_to_programmed_events('leave_home', ind_id)

        self.individuals[ind_id].assign_vaccination_status(self.scale_up_functions_current_time['bcg_coverage_prop'])
        self.individuals[ind_id].set_school_and_work_details(self.params)
        self.update_school_and_work_programs(ind_id)
        self.dates_of_birth[ind_id] = self.individuals[ind_id].dOB

    def add_event_to_programmed_events(self, event_type, ind_id):
        """
        :param event_type: One of "death", "leave_home",
        :param ind_id:
        :return:
        """
        event_type_individual = event_type
        if event_type == 'tb_death':
            event_type_individual = 'death'
        if self.individuals[ind_id].programmed[event_type_individual] in self.programmed_events[event_type].keys():
            self.programmed_events[event_type][self.individuals[ind_id].programmed[event_type_individual]].append(ind_id)
        else:
            self.programmed_events[event_type][self.individuals[ind_id].programmed[event_type_individual]] = [ind_id]

    def build_or_destroy_households(self):
        """
        This method will adjust the current number of households in order to match the average household size
        """
        nb_households_modified = False
        ideal_nb_hh = int(round(self.population/self.params['household_size']))
        current_n_hh = copy.copy(self.n_households)

        if current_n_hh > ideal_nb_hh and len(self.empty_households) > 0:  # we need to remove households so we try to remove some empty ones.
            nb_households_modified = True
            n_hh_to_remove = current_n_hh - ideal_nb_hh
            hh_ids_to_remove = [hh_id for i, hh_id in enumerate(self.empty_households) if i+1 <= n_hh_to_remove]
            for hh_id in hh_ids_to_remove:
                if hh_id in self.eligible_hh_for_birth.keys():
                    del self.eligible_hh_for_birth[hh_id]
                self.empty_households = [hh for hh in self.empty_households if hh != hh_id]
                del self.households[hh_id]
                self.n_households -= 1
        elif current_n_hh < ideal_nb_hh:  # we need to build new households
            nb_households_modified = True
            nb_hh_to_build = ideal_nb_hh - current_n_hh
            for _ in range(nb_hh_to_build):
                new_hh_id = max(self.households.keys()) + 1
                self.households[new_hh_id] = household.household(id=new_hh_id)
                self.empty_households.append(new_hh_id)
                self.households[new_hh_id].school_id = np.random.choice(self.groups_by_type['schools'], 1)[0]
                self.n_households += 1

        if nb_households_modified:
            self.adjust_nb_of_schools_and_workplaces()

    def update_want_to_move_list(self):
        """
        List the individuals who want to leave their current household and update the self.individuals_want_to_move
        dictionary accordingly.
         This method does not make individuals move home. It only add them to a queue.
        """
        for leave_home_time in [d for d in self.programmed_events['leave_home'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['leave_home'][leave_home_time]:
                self.individuals_want_to_move[ind_id] = self.time
            del self.programmed_events['leave_home'][leave_home_time]

    def trigger_individuals_move_home(self):
        """
        Make individuals move to another household if available. Individuals can only move by couple.
        Individuals that have been waited for longer are prioritised.
        """

        while len(self.individuals_want_to_move) > 1 and len(self.empty_households) > 0:
            ind_id_1 = min(self.individuals_want_to_move, key=self.individuals_want_to_move.get)
            del(self.individuals_want_to_move[ind_id_1])
            ind_id_2 = min(self.individuals_want_to_move, key=self.individuals_want_to_move.get)
            del (self.individuals_want_to_move[ind_id_2])
            self.make_couple_move_to_household(ind_ids=[ind_id_1, ind_id_2], hh_id=self.empty_households[0])
            del(self.empty_households[0])

    def make_couple_move_to_household(self, ind_ids, hh_id):
        """
        The two individuals listed in ind_ids are moving to the household hh_id. All relevant attributes are updated
        accordingly.
        :param ind_ids: a list of two individual ids
        :param hh_id: a household id
        """
        for ind_id in ind_ids:
            # update previous household attributes
            prev_hh_id = self.individuals[ind_id].household_id
            self.households[prev_hh_id].size -= 1
            self.households[prev_hh_id].individual_ids = [i for i in self.households[prev_hh_id].individual_ids
                                                          if i != ind_id]
            if self.households[prev_hh_id].size == 0:
                self.empty_households.append(prev_hh_id)

            # update individual's attributes
            self.individuals[ind_id].household_id = hh_id

        for ind_id in ind_ids:
            # update new household attributes
            self.households[hh_id].individual_ids.append(ind_id)

        # update new household's attributes
        self.households[hh_id].size += 2
        self.households[hh_id].repopulate_date = self.time

    def apply_interventions(self):
        """
        This method runs the different interventions at each steps
        """
        if self.params['mass_pt_program']:
            self.screen_and_treat_population_for_ltbi()

    def screen_and_treat_population_for_ltbi(self):
        """
        Run the mass preventive therapy program
        """
        screened_ind_ids = self.pick_screened_individuals()
        for ind_id in screened_ind_ids:
            # screening
            test_result = self.individuals[ind_id].test_individual_for_ltbi(self.params)

            # treatment
            if test_result:
                self.provide_preventive_treatment(ind_id, delayed=False)

    def pick_screened_individuals(self):
        """
        Creates a list of individual ids that are provided with screening for LTBI through the mass program.
        This program may be targeted at a specific subgroup.
        return: a list of ids
        """
        n_screened_individuals = round((self.params['mass_pt_screening_rate'] / 100.) * self.population \
                                       * self.params['time_step'] / 365.25)
        if self.params['subgroup_for_mass_pt'] == 'all':  # screening applies to all
            screened_individuals = np.random.choice(self.individuals.keys(), n_screened_individuals, replace=False)
        else:   # screening applies to specific groups
            ind_ids = range(int(self.population))
            np.random.shuffle(ind_ids)
            n_recorded = 0
            n_tried = 0
            screened_individuals = []
            while n_recorded < n_screened_individuals:
                ind_id = ind_ids[n_tried]
                if self.individuals[ind_id].is_in_subgroup(subgroup=self.params['subgroup_for_mass_pt'], time=self.time):
                    screened_individuals.append(ind_id)
                    n_recorded += 1
                n_tried += 1
                if n_tried == self.population + 1:
                    print "Could not find enough individuals to screen!"
                    exit('Program stopped')

        return screened_individuals

    def update_school_and_work_programs(self, ind_id):
        """
        Update the programmed_go_to_school programmed_leave_school and programmed_leave_work dictionaries.
        """
        keys_to_loop = ['go_to_school', 'leave_school']
        if self.individuals[ind_id].is_ever_gonna_work:
            keys_to_loop.append('leave_work')

        for key in keys_to_loop:
            if self.individuals[ind_id].programmed[key] in self.programmed_events[key].keys():
                self.programmed_events[key][self.individuals[ind_id].programmed[key]].append(ind_id)
            else:
                self.programmed_events[key][self.individuals[ind_id].programmed[key]] = [ind_id]

    def trigger_programmed_go_to_school(self):
        for go_to_school_time in [d for d in self.programmed_events['go_to_school'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['go_to_school'][go_to_school_time]:
                self.make_individual_go_to_school(ind_id)
            del self.programmed_events['go_to_school'][go_to_school_time]

    def make_individual_go_to_school(self, ind_id):
        # identify the group id
        group_id = self.households[self.individuals[ind_id].household_id].school_id
        self.make_individual_enter_group(ind_id, group_id)

    def trigger_programmed_leave_school(self):
        for leave_school_time in [d for d in self.programmed_events['leave_school'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['leave_school'][leave_school_time]:
                self.make_individual_leave_school(ind_id)
            del self.programmed_events['leave_school'][leave_school_time]

    def make_individual_leave_school(self, ind_id):
        for g_id in self.individuals[ind_id].group_ids:
            if self.group_types[g_id] == 'school':
                group_id = g_id
        # self.households[self.individuals[ind_id].household_id].school_id
        self.make_individual_leave_group(ind_id, group_id)
        if self.individuals[ind_id].is_ever_gonna_work:
            self.make_individual_go_to_work(ind_id)

    def make_individual_go_to_work(self, ind_id):
        group_id = np.random.choice(self.groups_by_type['workplaces'], 1)[0]
        self.make_individual_enter_group(ind_id, group_id)

    def trigger_programmed_leave_work(self):
        for leave_work_time in [d for d in self.programmed_events['leave_work'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['leave_work'][leave_work_time]:
                self.make_individual_leave_work(ind_id)
            del self.programmed_events['leave_work'][leave_work_time]

    def make_individual_leave_work(self, ind_id):
        group_id = self.individuals[ind_id].group_ids[0]
        self.make_individual_leave_group(ind_id, group_id)

    def make_individual_enter_group(self, ind_id, group_id):
        """
        Make individual "ind_id" enter the group "group_id"
        """
        self.individuals[ind_id].group_ids.append(group_id)
        self.groups[group_id].append(ind_id)

    def make_individual_leave_group(self, ind_id, group_id):
        """
        Make individual "ind_id" leave the group "group_id"
        """
        self.individuals[ind_id].group_ids = [gr_id for gr_id in self.individuals[ind_id].group_ids if gr_id != group_id]
        self.groups[group_id] = [ids for ids in self.groups[group_id] if ids != ind_id]

    def trigger_programmed_activations(self):
        """
        The individuals listed in the programmed_events['activation'] dictionary for the times elapsed since the last iteration time
        have to activate TB.
        """
        for activation_time in [d for d in self.programmed_events['activation'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['activation'][activation_time]:
                self.make_individual_activate_tb(ind_id)
            del self.programmed_events['activation'][activation_time]

    def make_individual_activate_tb(self, ind_id, init=False):
        """
        Triggers activation in the individual ind_id by changing the relevant attribute and updating the "active_cases"
        dictionary.
        """
        if self.individuals[ind_id].ltbi:  # may be ltbi-neg if pt has been used with delay or if initialisation of TB states.
            self.ltbi_prevalence -= 1
        self.active_cases.append(ind_id)
        self.individuals[ind_id].activate_tb()
        self.tb_prevalence += 1
        if self.individuals[ind_id].tb_strain == 'ds':
            self.tb_prevalence_ds += 1
        else:
            self.tb_prevalence_mdr += 1

        self.tb_incidence += 1

        tb_outcome = self.individuals[ind_id].define_tb_outcome(time=self.time, params=self.params,
                                                                tx_success_prop=self.scale_up_functions_current_time['treatment_success_prop'])
        if 'time_active' in tb_outcome.keys():
            self.time_active['total_n_cases'] += 1
            self.time_active['total_time_active'] += tb_outcome['time_active']

        if 'dr_amplification' in tb_outcome.keys():
            self.tb_prevalence_ds -= 1  # should be improved in the future as dr_amplification should occur later at treatment
            self.tb_prevalence_mdr += 1

        self.update_programmed_events(tb_outcome, ind_id)

        if self.params['plot_all_tb_ages'] and not init:
            self.all_tb_ages.append(self.individuals[ind_id].get_age_in_years(self.time))

    def record_all_contacts(self):
        if self.stopped_simulation:
            return

        print "Start recording contact patterns"

        n_recorded_index = int(ceil(self.population * self.params['perc_sampled_for_contacts'] / 100.))
        recorded_ind_ids = np.random.choice(self.individuals.keys(), n_recorded_index, replace=False)

        # initialisation
        for location in self.contact_matrices['contact'].keys():
            self.contact_matrices['contact'][location] = np.zeros((101, 101))  # null matrix 100x100
            self.n_contacts['contact'][location] = 0

        for i, ind_id in enumerate(recorded_ind_ids):
            contact_dict = self.get_contacts_during_last_step(ind_id, infectious_only=False)
            self.update_contact_matrices(ind_id, contact_dict)

            # count contacts
            for location in self.n_contacts['contact'].keys():
                self.n_contacts['contact'][location] += sum(contact_dict[location].values())

            if os.name == 'nt':  # only on local windows machine
                perc_completion = int(round(100*(i+1)/n_recorded_index, 0))
                dynamic_string = str(perc_completion) + "% of contacts recorded."
                stdout.write("\r\x1b[K" + dynamic_string.__str__())
                stdout.flush()
        # # rescale the contact matrices
        # for contact_type in self.contact_matrices.keys():
        #     self.contact_matrices[contact_type] /= n_contacts

    def get_contacts_during_last_step(self, ind_id, infectious_only=True):
        """
        Generate contacts between ind_id and other individuals of the whole population.
        If infectious_only is True, only the infectiousness period of ind_id intercepted with the last time step
        is considered. If not, the full time step is considered
        return: a dicitonary {id1: nb_contacts1, id2: nb_contacts2}
        """
        contact_dict = {}
        for location in self.contact_matrices['contact'].keys():
            contact_dict[location] = {}
        if infectious_only:
            start_recording = max(self.individuals[ind_id].programmed['activation'], self.time - self.params['time_step'])
            if 'recovery' not in self.individuals[ind_id].programmed.keys():
                end_recording = min(self.individuals[ind_id].programmed['death'], self.time)
            else:
                end_recording = min(self.individuals[ind_id].programmed['recovery'], self.individuals[ind_id].programmed['death'],
                                         self.time)
            recording_duration = max(end_recording - start_recording, 0)
        else:
            recording_duration = 1.  # we want to display daily contacts

        if recording_duration > 0:
            contact_dict['household'] = self.get_household_contacts(ind_id, recording_duration, infectious_only)
            contact_dict['community'] = self.get_community_contacts(ind_id, recording_duration, infectious_only)

            # does ind_id attend a group?
            if len(self.individuals[ind_id].group_ids) > 0:
                group_id = self.individuals[ind_id].group_ids[0]
                group_type = self.group_types[group_id]
                contact_dict[group_type] = self.get_network_contacts(ind_id, recording_duration, group_id=group_id,
                                                                     group_type=group_type, infectious_only=infectious_only)
        return contact_dict

    def get_household_contacts(self, ind_id, recording_duration, infectious_only=True):
        """
        Generate contacts between ind_id and people living in the same household. recording_duration is the nb of days
        during which contacts are recorded.
        return: a dictionary {id1: nb_contacts1, id2: nb_contacts2}
        """
        household_contacts = {}
        if self.households[self.individuals[ind_id].household_id].size > 1:
            contact_ids = self.get_other_household_members(ind_id)
            for c_id in contact_ids:
                # we assume that everyone living with ind_id contact ind_id once a day
                nb_contacts = np.random.binomial(recording_duration, self.params['perc_hh_contacted']/100.)
                if nb_contacts > 0:
                    household_contacts[c_id] = nb_contacts
                    if self.params['contact_tracing_pt_program'] and infectious_only:
                        self.individuals[ind_id].contacts_while_tb['household'].add(c_id)
        return household_contacts

    def get_network_contacts(self, ind_id, recording_duration, group_id, group_type, infectious_only=True):
        """
        Generate contacts between ind_id and people from the same group (school or workplace).
        recording_duration is the nb of days during which contacts are recorded.
        group_type is either 'school' or 'workplace'
        return: a dictionary {id1: nb_contacts1, id2: nb_contacts2}
        """
        network_contacts = {}
        age_ind_id = self.individuals[ind_id].get_age_in_years(self.time)
        if len(self.groups[group_id]) > 1:
            if group_type == 'workplace' and self.params['country'] != "None":
                sigma = self.sd_agepref_work[self.params['country']]
                total_contact_rate = self.prem_contact_rate_functions['work'](age_ind_id)
            else:
                sigma = self.params['sd_agepreference_school']
                total_contact_rate = self.prem_contact_rate_functions['school'](age_ind_id)

            # calculate f_sigma(x_i,j) for each individual
            f_sigmas = {}
            for potential_contact_id in self.groups[group_id]:
                 if potential_contact_id != ind_id:
                    age_contact_id = self.individuals[potential_contact_id].get_age_in_years(self.time)
                    f_sigmas[potential_contact_id] = age_preference_function(age_contact_id - age_ind_id, sigma)

            sum_f_sigmas = sum(f_sigmas.values())
            raw_proba_of_contact = total_contact_rate / sum_f_sigmas
            raw_proba_of_contact = (min(raw_proba_of_contact, 1.))

            for potential_contact_id in self.groups[group_id]:
                if potential_contact_id != ind_id:
                    nb_contacts = np.random.binomial(recording_duration, f_sigmas[potential_contact_id] *\
                                                     raw_proba_of_contact)
                    if nb_contacts > 0:
                        network_contacts[potential_contact_id] = nb_contacts
                        if self.params['contact_tracing_pt_program'] and infectious_only:
                            self.individuals[ind_id].contacts_while_tb[group_type].add(potential_contact_id)
        return network_contacts

    def get_community_contacts(self, ind_id, recording_duration, infectious_only=True):
        """
        Generate contacts between ind_id and random people from the community.
        return: a dictionary {id1: nb_contacts1, id2: nb_contacts2}
        """
        community_contacts = {}
        age_ind_id = self.individuals[ind_id].get_age_in_years(self.time)
        Prem_col_index = min(int(floor(age_ind_id / 5.)), 15)
        contact_rates = self.contact_rates_matrices['other_locations'][Prem_col_index, :] * recording_duration
        nb_of_contacts_to_draw = np.random.poisson(lam=contact_rates)  # nb of contacts per age category

        for i in range(16):
            if nb_of_contacts_to_draw[i] > 0.:
                age_cat = "X_" + str(i + 1)
                if len(self.ind_by_agegroup[age_cat]) <= nb_of_contacts_to_draw[i]:
                    contact_ids = self.ind_by_agegroup[age_cat]
                else:
                    contact_ids = np.random.choice(self.ind_by_agegroup[age_cat], size=nb_of_contacts_to_draw[i],
                                                   replace=False)
                for contact_id in contact_ids:
                    if contact_id != ind_id:
                        community_contacts[contact_id] = 1  # we assume a single contact
                        if self.params['contact_tracing_pt_program'] and infectious_only:
                            self.individuals[ind_id].contacts_while_tb['community'].add(contact_id)

        return community_contacts

    def get_daily_community_contact_rate(self, ind_id, age):
        """
        Returns the daily contact rate with individuals from outside of the household and outside of the network. This
        parameter depends on the age of the index individual.
        """
        raw_contact_rate = self.contact_rate_function(age)
        hh_size = self.households[self.individuals[ind_id].household_id].size
        adjusted_contact_rate = raw_contact_rate - float(hh_size - 1)
        return max(0., adjusted_contact_rate)

    def apply_transmission(self, contact_dict, relative_infectiousness, index_id):
        """
        At this stage, contacts of ind_id have been defined (information contained in contact_dict) but we still don't
        know whether they are associated with transmission. This method will trigger possible transmission events.
        """
        index_age = floor(self.individuals[index_id].get_age_in_years(self.time))
        for location in contact_dict.keys():
            for contacted_id, nb_contacts in contact_dict[location].iteritems():
                transmission_proba = self.params['proba_infection_per_contact'] *\
                                     self.individuals[contacted_id].get_relative_susceptibility(self.time, self.params) *\
                                     relative_infectiousness

                # relative contact fitness according to location
                transmission_proba *= self.params['rr_transmission_by_location'][location]

                success_proba = 1. - (1. - transmission_proba)**nb_contacts
                transmission = np.random.binomial(1, success_proba)

                if transmission == 1:  # transmission does occur
                    contact_age = floor(self.individuals[contacted_id].get_age_in_years(self.time))
                    self.n_contacts['transmission'][location] += 1
                    if index_age <= 100. and contact_age <= 100.:
                        self.contact_matrices['transmission'][location][int(index_age), int(contact_age)] += 1
                    # diseased (or future diseased) individuals are not affected with reinfection
                    if not self.individuals[contacted_id].active_tb and 'activation' not in \
                            self.individuals[contacted_id].programmed.keys():
                        if not self.individuals[contacted_id].ltbi:  # This is a newly infected individual
                            self.ltbi_prevalence += 1
                        self.infect_an_individual(contacted_id, strain=self.individuals[index_id].tb_strain)
                        if 'activation' in self.individuals[contacted_id].programmed.keys():  # responsible for a new TB case
                            self.n_contacts['transmission_end_tb'][location] += 1
                            if index_age <= 100. and contact_age <= 100.:
                                self.contact_matrices['transmission_end_tb'][location][int(index_age), int(contact_age)] += 1
                        if self.params['ideal_pt_program'] and self.time >= \
                                        365.25*(self.params['duration_burning_demo'] +
                                                    self.params['duration_burning_tb'] + self.params['intervention_start_delay']):  # PT is provided to all infectees, as soon as they get infected
                            self.provide_preventive_treatment(contacted_id, delayed=False)

    def update_contact_matrices(self, ind_id, contact_dict):
        """
        contact_dict contains the recorded contacts of individual 'ind_id' during the last time step.
        The contact matrix is updated accordingly.
        """
        age_ind_id = floor(self.individuals[ind_id].get_age_in_years(self.time))
        if age_ind_id <= 100.:
            for location in contact_dict.keys():
                for contacted_id, nb_contacts in contact_dict[location].iteritems():
                    age_contacted = floor(self.individuals[contacted_id].get_age_in_years(self.time))
                    if age_contacted <= 100.:
                        self.contact_matrices['contact'][location][int(age_ind_id), int(age_contacted)] += nb_contacts

    def trigger_programmed_deaths(self):
        """
        The individuals listed in the programmed_deaths dictionary for the times elapsed since the last iteration time
        have to die.
        """
        for death_time in [d for d in self.programmed_events['death'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['death'][death_time]:
                self.make_individual_die(ind_id)
            del self.programmed_events['death'][death_time]

    def trigger_programmed_tb_deaths(self):
        """
        The individuals listed in the programmed_deaths dictionary for the times elapsed since the last iteration time
        have to die.
        """
        for death_time in [d for d in self.programmed_events['tb_death'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['tb_death'][death_time]:
                self.make_individual_die(ind_id)
            del self.programmed_events['tb_death'][death_time]

    def make_individual_die(self, ind_id):
        """
        Reset the individual characteristics and generate a new dOD. This simulates the death of
        individual ind_id and generates a newborn in the same household. The newborn individual keeps
        the same ind_id
        """
        if self.individuals[ind_id].active_tb:
            self.tb_prevalence -= 1
            if self.individuals[ind_id].tb_strain == 'ds':
                self.tb_prevalence_ds -= 1
            else:
                self.tb_prevalence_mdr -= 1
            self.active_cases = [ind for ind in self.active_cases if ind != ind_id]
            self.tb_deaths += 1

        if self.individuals[ind_id].ltbi:
            self.ltbi_prevalence -= 1

        del self.dates_of_birth[ind_id]
        self.clean_programmed_dictionaries(ind_id)
        self.remove_from_groups(ind_id)

        age_cat = get_agecategory(self.individuals[ind_id].get_age_in_years(self.time), 'prem')
        if ind_id not in self.ind_by_agegroup[age_cat]: # the individual was still recorded in the previous age_category
            age_cat = get_agecategory(self.individuals[ind_id].get_age_in_years(self.time - 365.25), 'prem')
        self.ind_by_agegroup[age_cat] = [ind for ind in self.ind_by_agegroup[age_cat] if ind != ind_id]

        previous_hh_id = self.individuals[ind_id].household_id
        self.households[previous_hh_id].size -= 1

        self.households[previous_hh_id].individual_ids = [i for i in self.households[previous_hh_id].individual_ids \
                                                          if i != ind_id]
        self.population -= 1

        # The previous household may become empty and eligible for a new couple to move in
        if self.households[previous_hh_id].size == 0 and not self.households[previous_hh_id].retirement_home:
            if previous_hh_id in self.eligible_hh_for_birth.keys():
                self.households[previous_hh_id].repopulate_date = -1.e8
                del(self.eligible_hh_for_birth[previous_hh_id])
            self.empty_households.append(previous_hh_id)

        del self.individuals[ind_id]

        if self.fertility_replacement:
            self.make_individual_bear()

    def clean_programmed_dictionaries(self, ind_id):
        """
        Individual ind_id is about to die. We need to clean up a few dicitonaries.
        """

        keys_to_loop = self.programmed_events.keys()
        for key in keys_to_loop:
            if key in self.individuals[ind_id].programmed.keys():
                if self.individuals[ind_id].programmed[key] in self.programmed_events[key].keys():
                    self.programmed_events[key][self.individuals[ind_id].programmed[key]] = \
                        [ids for ids in self.programmed_events[key][self.individuals[ind_id].programmed[key]]
                         if ids != ind_id]

        if ind_id in self.individuals_want_to_move.keys():
            del(self.individuals_want_to_move[ind_id])

    def remove_from_groups(self, ind_id):
        for group_id in self.individuals[ind_id].group_ids:
            self.groups[group_id] = [ind for ind in self.groups[group_id] if ind != ind_id]

    def trigger_births(self):

        if self.constant_birth_rate:
            average_nb_births_per_step = self.params['birth_rate']*self.population*self.params['time_step'] / (365.25 * 1000.)
        else:  # age-pyramid driven
            time_to_pyramid = self.age_pyramid_date - self.time
            time_to_pyramid = abs(time_to_pyramid)   # not too clean but needed when one step goes over age_pyramid_date
            average_nb_births_per_step = self.birth_numbers_function(time_to_pyramid)

        nb_births = np.random.poisson(average_nb_births_per_step)
        for _ in repeat(None, nb_births):  # supposed to be faster than a classic for loop
            self.make_individual_bear()

    def make_individual_bear(self, ind_id=None):
        self.population += 1
        self.birth_numbers += 1
        hh_id = self.pick_eligible_household_for_birth()
        self.add_new_individual_in_hh(h_id=hh_id, age=0., ind_id=ind_id)
        if hh_id in self.empty_households:
            self.empty_households = [h for h in self.empty_households if h != hh_id]

        if hh_id in self.eligible_hh_for_birth.keys():
            del self.eligible_hh_for_birth[hh_id]

    def update_programmed_events(self, event_dict, ind_id=None):
        """
        update the programmed events dictionary.
        :param event_dict keys are the type of event ("death",...) and the values are the dates of the events
        """
        if 'death' in event_dict.keys():
            self.program_tb_death(ind_id, event_dict['death'])
        if 'recovery' in event_dict.keys():
            self.add_recovery_to_programmed_recoveries(ind_id, event_dict['recovery'])
        if 'detection' in event_dict.keys():
            self.add_detection_to_programmed_detections(ind_id, event_dict['detection'])

    def add_detection_to_programmed_detections(self, ind_id, detection_date):
        if detection_date in self.programmed_events['detection'].keys():
            self.programmed_events['detection'][detection_date].append(ind_id)
        else:
            self.programmed_events['detection'][detection_date] = [ind_id]

    def trigger_programmed_detections(self):
        """
        The individuals listed in the programmed_recoveries dictionary for the times elapsed since the last iteration time
        have to recover from TB.
        """
        for detection_time in [d for d in self.programmed_events['detection'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['detection'][detection_time]:
                self.detect_individual(ind_id)
            del self.programmed_events['detection'][detection_time]

    def detect_individual(self, ind_id):
        self.individuals[ind_id].detect_tb()
        if self.params['contact_tracing_pt_program'] and self.time >= (self.params['duration_burning_demo'] +
                                                                           self.params['duration_burning_tb'] +
                                                                           self.params['intervention_start_delay'] ) * 365.25:
            self.trigger_detection_linked_pt(ind_id)

    def trigger_detection_linked_pt(self, ind_id):
        """
        Trigger the implementation of (testing and) preventive therapy in contacts of ind_id.
        param ind_id: The index active case
        """
        if self.params['contact_type_for_contact_tracing_pt'] == 'all':
            relevant_contact_types = ['household', 'school', 'workplace']
        elif self.params['contact_type_for_contact_tracing_pt'] == 'network':
            relevant_contact_types = ['school', 'workplace']
        else:  # school or workplace or hosuehold
            relevant_contact_types = [self.params['contact_type_for_contact_tracing_pt']]

        all_identified_contacts = []
        for contact_type in relevant_contact_types:
            if self.params['agegroup_for_contact_tracing_pt'] == 'all':
                relevant_contacts = list(self.individuals[ind_id].contacts_while_tb[contact_type])
                relevant_contacts = [ind for ind in relevant_contacts if ind in self.individuals.keys()]
            else:
                relevant_contacts = []
                for contact_id in self.individuals[ind_id].contacts_while_tb[contact_type]:
                    if contact_id in self.individuals.keys():
                        if self.individuals[contact_id].is_in_subgroup(
                                subgroup=self.params['agegroup_for_contact_tracing_pt'], time=self.time):
                            relevant_contacts.append(contact_id)

            if len(relevant_contacts) > 0:
                nb_identified_contacts = int(round(len(relevant_contacts)*self.params['perc_coverage_tracing_' + contact_type]/100.))

                identified_contacts = np.random.choice(relevant_contacts, size=nb_identified_contacts, replace=False)
                for contact_id in identified_contacts:
                    if contact_id not in all_identified_contacts:
                        all_identified_contacts.append(contact_id)

        for contact_id in all_identified_contacts:
            shall_we_test = False
            if self.params['agegroup_tst_before_pt_in_contact'] == 'all':
                shall_we_test = True
            elif self.params['agegroup_tst_before_pt_in_contact'] == 'none':
                shall_we_test = False
            else:
                if self.individuals[contact_id].is_in_subgroup(subgroup=self.params['agegroup_tst_before_pt_in_contact'],
                                                               time=self.time):
                    shall_we_test = True

            if shall_we_test:
                ltbi_test = self.individuals[contact_id].test_individual_for_ltbi(self.params)
                if ltbi_test:
                    self.provide_preventive_treatment(contact_id, delayed=True)
            else:  # provide pt without testing
                self.provide_preventive_treatment(contact_id, delayed=False)

    def add_recovery_to_programmed_recoveries(self, ind_id, recovery_date):
        if recovery_date in self.programmed_events['recovery'].keys():
            self.programmed_events['recovery'][recovery_date].append(ind_id)
        else:
            self.programmed_events['recovery'][recovery_date] = [ind_id]

    def trigger_programmed_recoveries(self):
        """
        The individuals listed in the programmed_recoveries dictionary for the times elapsed since the last iteration time
        have to recover from TB.
        """
        for recovery_time in [d for d in self.programmed_events['recovery'].keys() if d <= self.time]:
            for ind_id in self.programmed_events['recovery'][recovery_time]:
                self.tb_prevalence -= 1
                if self.individuals[ind_id].tb_strain == 'ds':
                    self.tb_prevalence_ds -= 1
                else:
                    self.tb_prevalence_mdr -= 1
                self.make_individual_recover(ind_id)
            del self.programmed_events['recovery'][recovery_time]

    def make_individual_recover(self, ind_id):
        """
        Make an individual recover and update the relevant dictionaries
        """
        self.individuals[ind_id].recover()
        self.active_cases = [ind for ind in self.active_cases if ind != ind_id]

    def program_tb_death(self, ind_id, date_of_tb_death):
        """
        Given that the individual activates TB at the current time and knowing that the individual will die from TB, we
        record the new time of programmed death.
        """

        previous_dOD = self.individuals[ind_id].programmed['death']
        # remove previously programmed death (natural death)
        self.programmed_events['death'][previous_dOD] = [inds for inds in self.programmed_events['death'][previous_dOD]
                                                if inds != ind_id]
        self.individuals[ind_id].programmed['death'] = date_of_tb_death

        self.add_event_to_programmed_events('tb_death', ind_id)  # re-define the date of death

    def infect_an_individual(self, ind_id, strain):
        """
        Make individual "ind_id" infected and define the time to potential activation
        """
        self.individuals[ind_id].infect_individual(self.time, self.params, strain)
        self.activation_stats['n_infections'] += 1
        if 'activation' in self.individuals[ind_id].programmed.keys():
            self.add_activation_to_programmed_activations(ind_id)
            self.activation_stats['n_activations'] += 1

    def add_activation_to_programmed_activations(self, ind_id):
        if self.individuals[ind_id].programmed['activation'] in self.programmed_events['activation'].keys():
            self.programmed_events['activation'][self.individuals[ind_id].programmed['activation']].append(ind_id)
        else:
            self.programmed_events['activation'][self.individuals[ind_id].programmed['activation']] = [ind_id]

    def provide_preventive_treatment(self, ind_id, delayed=False):
        """
        Provide preventive treatment to the individual ind_id.
        In case of effective treatment in an individual who was meant to progress to active tb, the programmed date of
        activation of the individual has to be removed from the programmed_events dictionary.
        """
        self.n_pt_provided += 1.
        pre_ltbi = copy.copy(self.individuals[ind_id].ltbi)
        date_prevented_activation = self.individuals[ind_id].get_preventive_treatment(self.params, time=self.time, delayed=delayed)
        if date_prevented_activation is not None:  # The treatment is successful and useful
            self.programmed_events['activation'][date_prevented_activation] = [ids for ids in \
                                                                                self.programmed_events['activation'][
                                                                                    date_prevented_activation]
                                                                                                  if ids != ind_id]
            self.n_useful_pt_provided += 1.
        if pre_ltbi and not self.individuals[ind_id].ltbi:
            self.ltbi_prevalence -= 1

    def generate_checkpoint_outputs(self):
        """
        Perform some measures that have to be done at checkpoint times only
        """
        # Age distribution
        if self.stopped_simulation:   # we don't want to store ages for individuals when model has stopped running.
            self.checkpoint_outcomes['ages'][self.time] = []
        else:
            self.checkpoint_outcomes['ages'][self.time] = [
                individual.get_age_in_years(self.time) for
                individual in self.individuals.values()]

        # Household size distribition
        self.checkpoint_outcomes['household_sizes'][self.time] = [h.size for h in self.households.values() if
                                                                  not h.retirement_home]

        # School and worplace sizes
        for group_type in ['school', 'workplace']:
            self.checkpoint_outcomes[group_type + '_sizes'][self.time] = [len(group) for group_id, group in
                                                                          self.groups.iteritems() if
                                                                          self.group_types[group_id] == group_type]

    def get_ages_by_household(self):
        """
        Not used at the moment
        """
        ages_by_household = []

        for i, hh in self.households.iteritems():
            ages_by_household.append([])
            for ind_id in hh.individual_ids:
                ages_by_household[i].append(round(self.individuals[ind_id].get_age_in_years(self.time)))
        return ages_by_household

    def average_timeseries(self):
        """
        Apply a moving average smoothing to the timeseries that are listed in self.timeseries_to_average
        """
        n_steps = int(round(self.moving_average_width / self.params['time_step']))
        if n_steps > 0:
            for series_name in self.timeseries_log.keys():
                if series_name in self.timeseries_to_average:
                    y = copy.copy(self.timeseries_log[series_name])
                    for i, time in enumerate(self.timeseries_log['times']):
                        ind_min = max(0, i - n_steps)
                        ind_max = i  # min(len(y)-1, i + n_steps)
                        self.timeseries_log[series_name][i] = np.mean(y[ind_min:ind_max+1])

    def clean_timeseries(self):
        """
        Clean some time-series such as incidence that has a peak at the very first time-step.
        """
        if 'tb_incidence' in self.timeseries_log.keys():
            self.timeseries_log['tb_incidence'][0] = 0.

    def turn_model_into_dict(self):
        m_dict = {}
        attributes_to_store = ['timeseries_log', 'params', 'checkpoint_outcomes', 'i_seed', 'scenario', 'i_run',
                               'contact_matrices', 'all_tb_ages', 'ltbi_age_stats', 'n_contacts', 'tb_prevalence_by_age',
                               'stopped_simulation']
        for attr in attributes_to_store:
            m_dict[attr] = getattr(self, attr)
        return m_dict

    def write_status_file(self):
        # create the main directory using the name entered through the console
        time = datetime.now().strftime('%d-%m-%Y %H:%M:%S')
        str_to_write = "Year " + str(self.last_year_completed) + ' (' + str(int(self.get_current_date())) +\
                       ") completed at time: " + time +\
                       " Nb of TB cases: " + str(int(self.tb_prevalence)) + ' (' + \
                          str(int(round(1.e5 * self.tb_prevalence / self.population))) + ' /100k)'
        base_path = os.path.join('outputs')
        dir_path = os.path.join(base_path, self.params['project_name'])
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        file_path = os.path.join(dir_path, 'progress_seed' + str(self.i_seed) + '_' + self.scenario + '_run' + str(self.i_run) + '.txt')

        # Implement manual interruption of the program. If status file is deleted, running is stopped.
        if not os.path.exists(file_path) and self.status_file_created:
            self.stop_running_model()

        file = open(file_path, 'a+')
        file.write(str_to_write + "\n")
        file.close()
        self.status_file_created = True

    def check_workplace_ages(self):
        """
        For debugging only. Check that all individuals in workplaces are "adult"
        """
        print self.time
        for w in self.groups_by_type['workplaces']:
            for ind_id in self.groups[w]:
                print "________"
                print str(ind_id) + ": " + str(self.individuals[ind_id].get_age_in_years(self.time))

    def check_school_ages(self):
        """
        For debugging only. Check that all individuals in schools are aged between 5 and 18
        """
        print self.time
        for w in self.groups_by_type['schools']:
            for ind_id in self.groups[w]:
                if self.individuals[ind_id].get_age_in_years(self.time)<5 or self.individuals[ind_id].get_age_in_years(self.time)>18:
                    print "////////////////////"
                    print str(ind_id) + ": " + str(self.individuals[ind_id].get_age_in_years(self.time))

    def check_calibration_targets(self):
        current_date = self.get_current_date()
        for y in self.remaining_calibration_targets.keys():
            if y <= current_date:
                for target in self.remaining_calibration_targets[y]:
                    if not self.is_target_verified(target):
                        self.stop_running_model()
                del self.remaining_calibration_targets[y]

    def is_target_verified(self, target):
        pass_test = False
        print "Check target for year " + str(target['year'])
        if 'category' in target.keys():
            if target['indicator'] == 'tb_prevalence':
                abs_prev = 0
                for ind_id in self.active_cases:
                    if '_smearpos' in target['category'] and self.individuals[ind_id].tb_organ != '_smearpos':
                        continue
                    if '_pulmonary' in target['category'] and self.individuals[ind_id].tb_organ == '_extrapulmonary':
                        continue
                    if 'more_than_15' in target['category'] and self.individuals[ind_id].get_age_in_years(self.time) < 15.:
                        continue
                    abs_prev += 1
                nb_kids = 0
                for age_cat in ['X_1', 'X_2', 'X_3']:
                    nb_kids += len(self.ind_by_agegroup[age_cat])
                deno = self.population - nb_kids
                model_measure = abs_prev * 1.e5 / deno
                print model_measure
        else:
            model_measure = self.timeseries_log[target['indicator']][-1]
        if target['min_accepted_value'] <= model_measure <= target['max_accepted_value']:
            pass_test = True
            print "Target checked"
        else:
            print "Model estimate for " + target['indicator'] + " in " + str(target['year']) + " must be in (" +\
                  str(target['min_accepted_value']) + "-" + str(target['max_accepted_value']) + "). Value: " +\
                  str(model_measure)
        return pass_test

    def store_calibrated_model(self):
        """
        If the calibration conditions are verified, the fully run model is stored
        :return:
        """
        print "!!!!!!!!!  We have a calibrated model !!!!!!!!!"
        print "Saving model run..."
        dir_path = os.path.join('calibrated_models')
        # dir_path = os.path.join(base_path, self.params['project_name'], self.scenario)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        second_dir_path = os.path.join(dir_path, self.params['project_name'])
        if not os.path.exists(second_dir_path):
            os.makedirs(second_dir_path)
        file_name_base = "pickled_calibrated_model_" + self.scenario + "_run" + str(self.i_run) + ".pickle"
        file_name = os.path.join(second_dir_path, file_name_base)
        if os.path.isfile(file_name):
            file_name_base = str(np.random.randint(1, 100)) + "_pickled_calibrated_model_" + self.scenario + "_run" + str(self.i_run) + ".pickle"
            file_name = os.path.join(second_dir_path, file_name_base)
        file_stream = open(file_name, "wb")

        del self.scale_up_functions
        del self.birth_numbers_function
        dill.dump(self, file_stream)
        file_stream.close()
        print "Calibrated model successfully saved for scenario: " + self.scenario
        self.has_been_stored = True
        if self.params['stop_running_after_calibration'] and self.params['running_mode'] != 'run_ks_based_calibration':
            print "Simulation will now be stopped as a calibrated model has been found."
            file_name = os.path.join('outputs', self.params['project_name'], 'keep_running.txt')
            os.remove(file_name)
            exit()

    def get_current_date(self):
        remaining_year = (self.age_pyramid_date - self.time)/365.25
        time_in_years = self.params['current_year'] - remaining_year
        return time_in_years


class TbModel(Model):
    def __init__(self, data, i_seed, scenario, i_run, initialised=True):

        Model.__init__(self, data, i_seed, scenario, i_run, initialised)

        self.active_cases = []
        self.all_tb_ages = []
        self.tb_prevalence_by_age = []
        self.transmission = False  # will become the same as the inputed parameter when demographic burning is done
        self.tb_has_started = False

        self.ltbi_prevalence = 0.
        self.tb_prevalence = 0.
        self.tb_prevalence_ds = 0.
        self.tb_prevalence_mdr = 0.
        self.tb_incidence = 0.  # Absolute number of new cases for current time-step. Reset at each new time-step.
        self.tb_deaths = 0.  # Absolute number of tb-deaths for current time-step. Reset at each new time-step.
        self.time_active = {'total_n_cases': 0, 'total_time_active': 0}

        self.n_pt_provided = 0  # absolute number of prev treatments provided during the current time-step
        self.n_useful_pt_provided = 0  # absolute number of prev treatments leading to infection cure

        self.programmed_events.update({'activation': {}, 'detection': {}, 'recovery': {}, 'tb_death': {}})

        # attributes pertaining to the series that need smoothing by moving average
        self.timeseries_to_average = ['tb_incidence', 'tb_deaths', 'n_pt_provided', 'n_useful_pt_provided']
        self.moving_average_width = 365.  # number of days considered before each time for averaging timeseries / hard-coded

        self.initialised = initialised
        if not initialised:
            self.initialise_model(data)

    def spread_infections(self):
        """
        Rules the whole transmission process.
        """
        for ind_id in self.active_cases:
            relative_infectiousness = self.individuals[ind_id].get_relative_infectiousness(self.params, self.time)
            if relative_infectiousness > 0.:  # the index case is infectious
                contact_dict = self.get_contacts_during_last_step(ind_id)  # returns a dictionary keyed with contact ids, valued with nb of contacts
                self.apply_transmission(contact_dict, relative_infectiousness, index_id=ind_id)

    def process_cdr(self):
        """
        The case detection rate (cdr) is provided as an input parameter. We need to calculate the rate of the exponential
        distribution associated with the time to detection which would yield the specified cdr.
        """
        for organ in ['_smearpos', '_closed_tb']:
            # cdr = self.params['perc_cdr' + organ] / 100.
            cdr = self.scale_up_functions_current_time['cdr_prop']
            mu = 0.01   # hard-coded natural mortality for CDR calculation
            if cdr > 0.95:
                print "WARNING: a CDR too close to 100% will lead to no contact identified as detection occurs very quickly"
            assert cdr <= 1., "Case detection must be <= 1"
            if cdr == 1.:
                self.params['lambda_timeto_detection' + organ] = 1.e9  # some big value
            elif cdr == 0.:
                self.params['lambda_timeto_detection' + organ] = 1. / 1.e9 # some tiny value
            else:
                self.params['lambda_timeto_detection' + organ] = \
                    (cdr / (1. - cdr)) * (self.params['rate_self_cure' + organ] +
                                          self.params['rate_tb_mortality' + organ] +
                                          mu)

    def process_organ_proportions(self):
        self.params['perc_smearpos'] = 100. * self.scale_up_functions_current_time['sp_prop']
        self.params['perc_extrapulmonary'] = 0.5 * (100. - self.params['perc_smearpos'])

    def adjust_attributes_after_calibration(self):
        """
        We need to reset some attributes to get the model ready for the recording / analysis phase
        """
        self.params['plot_contact_heatmap'] = True
        self.params['plot_all_tb_ages'] = True
        self.reset_recording_attributes()

    def reset_recording_attributes(self):
        self.all_tb_ages = []
        self.n_pt_provided = 0  # absolute number of prev treatments provided during the current time-step
        self.n_useful_pt_provided = 0  # absolute number of prev treatments leading to infection cure
        self.initialise_contact_matrices()
        self.ltbi_age_stats = {'ltbi_ages': [], 'ending_tb_ages': []}
        self.ltbi_age_stats_have_been_recorded = False

    def record_ltbi_ages(self):
        for individual in self.individuals.values():
            if individual.ltbi:
                age = individual.get_age_in_years(self.time)
                self.ltbi_age_stats['ltbi_ages'].append(age)
                if 'activation' in individual.programmed.keys():
                    self.ltbi_age_stats['ending_tb_ages'].append(age)

    def record_tb_prevalence_by_age(self):
        self.tb_prevalence_by_age = []
        age_breaks = [0., 5., 10., 15., 25., 35., 45., 55., 65.]
        age_cats = [['X_1'], ['X_2'], ['X_3'], ['X_4', 'X_5'], ['X_6', 'X_7'],  ['X_8', 'X_9'],  ['X_10', 'X_11'],
                    ['X_12', 'X_13'], ['X_14', 'X_15', 'X_16']]
        nb_cases = [0. for _ in range(len(age_breaks))]
        agegroup_size = []

        for i in range(len(age_breaks)):
            pop_size = 0.
            for age_cat in age_cats[i]:
                pop_size += len(self.ind_by_agegroup[age_cat])
            agegroup_size.append(pop_size)

        # calculate the absolute prevalence by age
        for tb_case_id in self.active_cases:
            tb_age = self.individuals[tb_case_id].get_age_in_years(self.time)
            if tb_age >= 65.:
                age_cat_index = 8
            else:
                age_cat_index = next(x[0] for x in enumerate(age_breaks) if tb_age <= x[1]) - 1
            nb_cases[age_cat_index] += 1.

        # make the prevalence as relative to age-group pop sizes
        for i in range(len(nb_cases)):
            if agegroup_size[i] > 0:
                self.tb_prevalence_by_age.append( 1.e5 * nb_cases[i]/agegroup_size[i])  # now /100,000
            else:
                self.tb_prevalence_by_age.append(0.)



