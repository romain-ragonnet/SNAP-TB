from importData import read_sheet, sheet_to_dict
import model_runner
import outputs
import copy
import time
import os
from multiprocessing import Pool
from numpy import random, linspace

max_ncores = 12

t_0 = time.time()
last_i_figure = 0
calibration_params = {None: {'proba_infection_per_contact': linspace(start=0.002, stop=0.0025, num=2)},
                        'India': {'proba_infection_per_contact': linspace(start=0.0014, stop=0.0014, num=1)},
                        'Indonesia': {'proba_infection_per_contact': linspace(start=0.0017, stop=0.0017, num=1)},
                        'China': {'proba_infection_per_contact': linspace(start=0.0014, stop=0.0014, num=1)},
                        'Philippines': {'proba_infection_per_contact': linspace(start=0.0021, stop=0.0021, num=1)},
                        'Pakistan': {'proba_infection_per_contact': linspace(start=0.0016, stop=0.0016, num=1)}
                      }

# for LHS calibration method:
uncertainty_params = {'proba_infection_per_contact': {'distri': 'uniform', 'pars': (0.0030, 0.0045)},
                      'infectiousness_switching_age': {'distri': 'triangular', 'pars': (10., 20.)},
                      'g_child': {'distri': 'beta', 'pars': (73.27, 40.46)},
                      'g_teen': {'distri': 'beta', 'pars': (147., 34.84)},
                      'g_adult': {'distri': 'beta', 'pars': (584.30, 28.53)},
                      'rate_sp_cure_smearpos': {'distri': 'triangular', 'pars': (0.16, 0.31)},
                      'rate_sp_cure_closed_tb': {'distri': 'triangular', 'pars': (0.05, 0.25)},
                      'rate_tb_mortality_smearpos': {'distri': 'triangular', 'pars': (0.31, 0.47)},
                      'rate_tb_mortality_closed_tb': {'distri': 'triangular', 'pars': (0.013, 0.039)},
                      'time_to_treatment': {'distri': 'triangular', 'pars': (0., 14.)},
                      'n_colleagues': {'distri': 'triangular', 'pars': (10., 30.)}
                      }


# read country(ies)
file_path = os.path.join('spreadsheets', 'console.xlsx')
sheet = read_sheet(file_path)
par_dict = sheet_to_dict(sheet)
countries = par_dict['country']
country_list = countries.split('/')
load_calibrated_models = par_dict['load_calibrated_models']
calibrated_models_directory = par_dict['calibrated_models_directory']
running_mode = par_dict['running_mode']
del par_dict

for country in country_list:
    if country is not None:
        print "******************************"
        print "Running model for " + country

    if country in calibration_params.keys():
        m_r = model_runner.ModelRunner(country=country, calibration_params=calibration_params[country],
                                       uncertainty_params=uncertainty_params)
    else:
        m_r = model_runner.ModelRunner(country=country, calibration_params=calibration_params[None],
                                       uncertainty_params=uncertainty_params)
    m_r.clear_output_dir()

    parallel_processing = False
    if m_r.data.console['n_runs'] * len(m_r.data.scenario_names) * m_r.nb_seeds > 1 and os.name != 'nt':
        parallel_processing = True

    def run_a_single_simulation(run_indices):
        """
        run_indices is a list (seed_index, scenario, i_run)
        """
        # Is keep_running.txt file still there?
        file_path = os.path.join(m_r.base_path, m_r.data.console['project_name'], 'keep_running.txt')
        if not os.path.exists(file_path):
            exit('Process exit from test.py: "keep_running" file was deleted')

        seed_index = run_indices[0]
        scenario = run_indices[1]
        i_run = run_indices[2]
        print "Running " + scenario + " run " + str(i_run)
        random.seed(i_run)
        m = copy.deepcopy(m_r.m_init[seed_index][scenario])
        m.i_seed = seed_index
        m.i_run = i_run
        m.initialised = True
        m.run()
        print "__________________________ " + scenario + " seed " + str(seed_index) + " run " + str(i_run) + " successfully run"

        if os.name != 'nt' and running_mode == 'run_lhs_calibration':
            mo_dict = {}
        else:
            mo_dict = m.turn_model_into_dict()

        del m
        return mo_dict

    run_indices = []
    nb_seeds = m_r.nb_seeds
    for seed_index in range(nb_seeds):
        for scenario in m_r.data.scenario_names:
            for i_run in range(m_r.data.console['n_runs']):
                run_indices.append((seed_index, scenario, i_run))

    if __name__ == '__main__':

        if parallel_processing:
            p = Pool(processes=min(m_r.n_cpus, m_r.data.console['n_runs'] * len(m_r.data.scenario_names) * m_r.nb_seeds,
                     max_ncores))
            output_models = p.map(func=run_a_single_simulation, iterable=run_indices)
        else:
            output_models = []
            for indices in run_indices:
                m_dict = run_a_single_simulation(indices)
                output_models.append(m_dict)

        if os.name != 'nt' and running_mode == 'run_lhs_calibration':
            print "No output created as not requested on remote machines"
            del m_r
            del output_models
        else:
            for m_dict in output_models:
                m_r.store_a_model_run(m_dict)

            m_r.aggregate_scenario_results()

            print 'Simulation completed'

            O = outputs.output(m_r, last_i_figure)
            del m_r
            del output_models

            if os.name == 'nt':
                O.make_graphs()
                O.write_timeseries_to_csv()

            last_i_figure = O.i_figure
            del O

elapsed = time.time() - t_0
print "The whole simulation took " + str(elapsed) + " seconds to complete."
