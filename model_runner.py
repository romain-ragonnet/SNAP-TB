import importData as imp
import numpy as np
import model
import copy
import dill
from scipy import stats
from multiprocessing import cpu_count
import os, shutil

class ModelRunner:
    """
    Runs the simulation. Allows several scenarios and several runs to be launched
    """
    def __init__(self, country=None, calibration_params=None, uncertainty_params=None):
        # country is only used for the purpose of multi-country analysis
        self.country = country
        self.calibration_params = calibration_params
        self.uncertainty_params = uncertainty_params

        self.data = imp.data(country, calibration_params, uncertainty_params) # read all the data
        self.base_path = os.path.join('outputs')
        self.model_diagnostics = {} # will store the model instances for the different scenarios and different runs
        self.m_init = []  # contains initialised models used for the different runs of the diff scenarios.
                          # list of dictionaries. One dict for each "loaded seed". Dict keys are scenarios.
                          # If no seed loaded, only one item in the list
        self.paths_to_calibrated_models = []
        self.nb_seeds = 1
        self.n_cpus = cpu_count()

        self.create_keep_running_file()
        self.initialise_simulation()
        print "########## The simulation has been successfully initialised  ##########"

    def clear_output_dir(self):
        """
        This method removes all files and directories contained in the directory outputs/<project_name>
        """
        folder = os.path.join(self.base_path, self.data.console['project_name'])

        if os.path.exists(folder):
            for the_file in os.listdir(folder):
                file_path = os.path.join(folder, the_file)
                if '.pickle' not in file_path and "keep_running" not in file_path and "lhs_values" not in file_path:
                    try:
                        if os.path.isfile(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path): shutil.rmtree(file_path)
                    except Exception as e:
                        print(e)

    def create_keep_running_file(self):
        """
        Create an empty file named "keep_running" in the project output directory. If this file disappear during
        simulation, all simulations related to the project are stopped.
        """
        folder = os.path.join(self.base_path, self.data.console['project_name'])
        if not os.path.exists(folder):
            os.makedirs(folder)
        file_path = os.path.join(folder, 'keep_running.txt')
        file = open(file_path, 'a+')
        file.write("Currently running...")
        file.close()

    def initialise_simulation(self):
        # create the project output directory

        if self.data.console['load_calibrated_models']:
            # we might have several seed models to load
            self.get_list_of_calibrated_models()
            self.nb_seeds = len(self.paths_to_calibrated_models)

        self.m_init = [{} for _ in range(self.nb_seeds)]

        # initialise m_init storage
        for seed_index in range(self.nb_seeds):
            if not self.data.console['different_init']:
                # initialise a common model that will be used as a base population for the different scenarios and different runs
                if self.data.console['load_root_models']:
                    m_init = self.load_model(self.data.scenarios.keys()[0])
                elif self.data.console['load_calibrated_models']:
                    print 'Loading seed ' + str(seed_index) + " ..."
                    m_init = self.load_model(self.data.scenarios.keys()[0], calibrated=True, seed_index=seed_index)
                    m_init.adjust_attributes_after_calibration()
                    print '... done'
                else:
                    m_init = model.TbModel(self.data, i_seed=seed_index, scenario='init', i_run=-1, initialised=False)
                for scenario in self.data.scenarios:
                    m_init.scenario = scenario
                    m_init.reset_params(self.data)
                    m_init.collect_params(self.data)
                    m_init.collect_scenario_specific_params(self.data)
                    if self.data.console['load_calibrated_models']:
                        m_init.adjust_attributes_after_calibration()
                    self.m_init[seed_index][scenario] = copy.deepcopy(m_init)
                    if self.data.console['store_root_models']:
                        self.store_model(m_init)
            else:         # If we need to initialise a different model for each scenario
                for scenario in self.data.scenarios:
                    if self.data.console['load_root_models']:
                        m_init = self.load_model(scenario)
                    elif self.data.console['load_calibrated_models']:
                        m_init = self.load_model(scenario, calibrated=True, seed_index=seed_index)
                    else:
                        m_init = model.TbModel(self.data, i_seed=seed_index, scenario=scenario, i_run=-1, initialised=False)

                    m_init.reset_params(self.data)
                    m_init.collect_params(self.data)
                    m_init.collect_scenario_specific_params(self.data)

                    if self.data.console['load_calibrated_models']:
                        print 'Loading seed ' + str(seed_index) + " ..."
                        m_init.adjust_attributes_after_calibration()
                        self.m_init[seed_index][scenario] = copy.deepcopy(m_init)
                        print "... done"
                    if self.data.console['store_root_models']:
                        self.store_model(m_init)

        # print message in the console
        string = ''
        for scenario in self.data.scenario_names:
            string += scenario + " "
        string += "will be run from " + str(self.nb_seeds) + " seeds, " + str(self.data.console['n_runs']) + \
                  " times each for " + str(self.data.console['n_years']) + " years."
        print string

        self.initialise_storage()  # initialise diagnostics storage

    def initialise_storage(self):
        """
        create storage variables to store scenario specific data
        """
        for scenario in self.data.scenario_names:
            self.model_diagnostics[scenario] = {'timeseries': {},
                                                'checkpoint_outcomes': {},
                                                'contact_matrices': {},
                                                'n_contacts': {},
                                                'all_tb_ages': [],
                                                'ltbi_age_stats': {'ltbi_ages': [], 'ending_tb_ages': []},
                                                'tb_prevalence_by_age': [],
                                                'n_accepted_runs': 0
                                                }
            for key in ['contact', 'transmission', 'transmission_end_tb']:
                self.model_diagnostics[scenario]['contact_matrices'][key] = {}
                self.model_diagnostics[scenario]['n_contacts'][key] = {}
                for location in ['school', 'workplace', 'household', 'community']:
                    self.model_diagnostics[scenario]['contact_matrices'][key][location] = np.zeros((101, 101))
                    self.model_diagnostics[scenario]['n_contacts'][key][location] = []

    def get_list_of_calibrated_models(self):
        dir_path = os.path.join('calibrated_models', self.data.console['calibrated_models_directory'] + "_" +
                                self.data.country)
        list_of_items = os.listdir(dir_path)
        self.paths_to_calibrated_models = [filename for filename in list_of_items if 'pickled' in filename]

    def store_a_model_run(self, m_dict):
        """
        Populate the diagnostics dictionaries of model_runner with the outputs of model m
        """

        # stored by run: timeseries/ nb of contacts and transmission events by location
        # stored by cumulating over runs: checkpoints / heatmaps / tb ages / ltbi ages stats

        # calculate row_index accounting for seed index and run_index
        row_index = m_dict['i_run'] + m_dict['i_seed'] * self.data.console['n_runs']

        # time series:
        for name, series in m_dict['timeseries_log'].iteritems():
            if name not in self.model_diagnostics[m_dict['scenario']]['timeseries'].keys():
                # initialize an array with the right dimensions
                n_row = m_dict['params']['n_runs'] * self.nb_seeds
                n_col = len(m_dict['timeseries_log']['times'])
                self.model_diagnostics[m_dict['scenario']]['timeseries'][name] = np.zeros(shape=(n_row, n_col))
            # store a time series

            self.model_diagnostics[m_dict['scenario']]['timeseries'][name][row_index, ] = series

        # chekpoint outcomes
        for name, series in m_dict['checkpoint_outcomes'].iteritems():
            if name not in self.model_diagnostics[m_dict['scenario']]['checkpoint_outcomes'].keys():
                # initialize a dictionary
                self.model_diagnostics[m_dict['scenario']]['checkpoint_outcomes'][name] = {}
                for checkpoint in m_dict['params']['checkpoints']:
                    self.model_diagnostics[m_dict['scenario']]['checkpoint_outcomes'][name][checkpoint] = {}
            for checkpoint in m_dict['params']['checkpoints']:
                if checkpoint in series.keys():
                    self.model_diagnostics[m_dict['scenario']]['checkpoint_outcomes'][name][checkpoint][row_index] = \
                        series[checkpoint]
                else:
                    print "Warning: checkpoint " + str(checkpoint) + " is not among the iteration times."

        # contact/transmission matrices
        for key in self.model_diagnostics[m_dict['scenario']]['contact_matrices'].keys():
            for location in self.model_diagnostics[m_dict['scenario']]['contact_matrices'][key].keys():
                self.model_diagnostics[m_dict['scenario']]['contact_matrices'][key][location] += \
                    m_dict['contact_matrices'][key][location]

        # all tb ages
        if m_dict['params']['plot_all_tb_ages']:
            self.model_diagnostics[m_dict['scenario']]['all_tb_ages'].append(m_dict['all_tb_ages'])

        # ltbi ages stats
        if m_dict['params']['plot_all_ltbi_ages']:
            for key in ['ltbi_ages', 'ending_tb_ages']:
                self.model_diagnostics[m_dict['scenario']]['ltbi_age_stats'][key].append(m_dict['ltbi_age_stats'][key])

        # nb of transmission events by location
        for key in self.model_diagnostics[m_dict['scenario']]['n_contacts'].keys():
            for location in self.model_diagnostics[m_dict['scenario']]['n_contacts'][key].keys():
                self.model_diagnostics[m_dict['scenario']]['n_contacts'][key][location].append(m_dict['n_contacts'][key][location])

        # prevalence by age
        self.model_diagnostics[m_dict['scenario']]['tb_prevalence_by_age'].append(m_dict['tb_prevalence_by_age'])

        # number of calibrated models
        if not m_dict['stopped_simulation']:
            self.model_diagnostics[m_dict['scenario']]['n_accepted_runs'] += 1

    def aggregate_scenario_results(self):
        """
        populate the model_diagnostics dictionary with the aggregated outputs for the different scenarios
        """
        for scenario in self.data.scenarios:
            # initialise storage
            self.model_diagnostics[scenario]['aggr_timeseries'] = {}
            self.model_diagnostics[scenario]['aggr_checkpoint_outcomes'] = {}
            # Timeseries
            for name, series_array in self.model_diagnostics[scenario]['timeseries'].iteritems():
                self.model_diagnostics[scenario]['aggr_timeseries'][name] = {'mean': [], 'low': [], 'high': []}
                mean = np.mean(series_array, axis=0)
                std = np.std(series_array, axis=0)
                # to avoid null std for samples that are duplicates
                std = [max(sd, 1.0e-9) for sd in std]
                intervals = stats.norm.interval(0.95, loc=mean, scale=std / np.sqrt(self.data.console['n_runs']*self.nb_seeds))
                self.model_diagnostics[scenario]['aggr_timeseries'][name]['mean'] = mean
                self.model_diagnostics[scenario]['aggr_timeseries'][name]['low'] = intervals[0]
                self.model_diagnostics[scenario]['aggr_timeseries'][name]['high'] = intervals[1]

            # Checkpoint outcomes
            # To be defined

    def open_output_directory(self):
        folder = os.path.join(self.base_path, self.data.console['project_name'])
        os.startfile(folder)

    def store_model(self, model_to_store):
        print "Storing root model for " + model_to_store.scenario + "..."
        file_name = "pickled_model_" + model_to_store.scenario + ".pickle"
        file_path = os.path.join(self.base_path, self.data.console['project_name'], file_name)
        file_stream = open(file_path, "wb")

        saved_scale_up_functions = copy.deepcopy(model_to_store.scale_up_functions)
        saved_birth_numbers_function = copy.deepcopy(model_to_store.birth_numbers_function)

        del model_to_store.scale_up_functions
        del model_to_store.birth_numbers_function

        dill.dump(model_to_store, file_stream)
        file_stream.close()

        model_to_store.scale_up_functions = saved_scale_up_functions
        model_to_store.birth_numbers_function = saved_birth_numbers_function
        print "Complete."

    def load_model(self, scenario, calibrated=False, seed_index=0):
        if calibrated:

            base_path = os.path.join('calibrated_models', self.data.console['calibrated_models_directory'] + "_" +
                                    self.data.country)
            file_name = self.paths_to_calibrated_models[seed_index]
            print "Loading calibrated model from " + self.data.console['calibrated_models_directory'] + "_" \
                  + self.data.country + "/" + file_name + " ..."
        else:
            print "Loading root model for " + scenario + "..."
            base_path = os.path.join('outputs', self.data.console['project_name'])
            file_name = "pickled_model_" + scenario + ".pickle"

        file_path = os.path.join(base_path, file_name)
        file_stream = open(file_path, "rb")
        loaded_model = dill.load(file_stream)
        loaded_model.scale_up_functions = self.data.scale_up_functions
        loaded_model.birth_numbers_function = self.data.birth_numbers_function
        file_stream.close()
        print "Complete."
        return loaded_model



if __name__ == "__main__":

    file_path = 'calibrated_models/test_LHScalibration_100s12r_Jul2019_Philippines/pickled_calibrated_model_lhs_sample_0_run11.pickle'
    file_stream = open(file_path, "rb")
    loaded_model = dill.load(file_stream)

    loaded_model.params['g_child']
    print "s"
