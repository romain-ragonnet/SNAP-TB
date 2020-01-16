import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as patches
from matplotlib.colors import LightSource
import matplotlib.ticker as ticker

import importData as imp

plt.style.use('ggplot')
import math
import numpy as np
import csv
import dill
import copy
import csv
if os.name == 'nt':
    from seaborn import violinplot
from scipy import stats

class output:

    def __init__(self, model_runner, last_i_figure=0):
        self.model_runner = model_runner
        self.base_path = os.path.join('outputs') # will be updated once directories are created
        self.i_figure = last_i_figure
        self.dictionary = {'ltbi_prevalence': 'LTBI prevalence (%)',
                            'tb_prevalence': 'TB prevalence (/100,000 population)',
                            'tb_prevalence_ds': 'DS-TB prevalence (/100,000 population)',
                            'tb_prevalence_mdr': 'MDR-TB prevalence (/100,000 population)',
                            'prop_mdr_prevalence': 'proportion of MDR-TB among prevalent TB cases',
                            'tb_incidence': 'TB incidence (/100,000 population/y)',
                            'tb_deaths': 'TB mortality (/100,000 population/y)',
                            'mean_age': 'mean age of the population',
                            'prop_under_5': 'poportion of population under 5 y.o.',
                            'n_pt_provided': 'number of prev. treat. provided',
                            'n_useful_pt_provided': 'number of useful prev. treat.',
                            'population': 'population',
                            'birth_rate': 'birth rate (/1,000 population/year)',
                            'proba_infection_per_contact': 'probability of infection per contact',
                            'cdr_prop': 'case detection proportion (%)',
                            'treatment_success_prop': 'treatment success rate (%)',
                            'bcg_coverage_prop': 'BCG coverage (%)',
                            'sp_prop': 'proportion of smear-positive cases (%)'
                           }
        self.timeseries_to_plot = self.model_runner.data.console['timeseries_to_record']

        self.scenario_colors = 'bgrcmykw'

        self.create_directories()

        if self.model_runner.data.console['running_mode'] != 'find_a_calibrated_model':
            self.store_outputs()

    def make_sensitivity_plot(self):
        param_vals = get_lhs_param_values('calibrated_models/test_LHScalibration_100s12r_Jul2019_' +
                                          self.model_runner.country + '/all_lhs_values.csv')

        print param_vals.keys()
        translate = {'proba_infection_per_contact': 'transmission probability', 'n_colleagues': '# of colleagues',
                     'infectiousness_switching_age': 'infectiousness switching age (a)',
                     'g_child': "g_child", 'g_teen': 'g_teen', 'g_adult': 'g_adult'}

        self.i_figure += 1
        plt.figure(self.i_figure, figsize=(3.7, (5 + 0.5)*0.9))

        h_ratios = [0.5, 5., 5., 5., 5., 5.]

        gs = gridspec.GridSpec(6, 4,
                               width_ratios=[0.7, 1., 1., 1.],
                               height_ratios=h_ratios
                               )
        gs.update(wspace=0.12, hspace=0.1)

        #  panel titles
        # titles = ['Contacts', 'Transmission', 'Transmission leading\nto active TB']
        # row = 0
        # for i, title in enumerate(titles):
        #     plt.subplot(gs[row, i+1])
        #     plt.text(x=0.5, y=0.9, s='', fontsize=5, horizontalalignment='center',
        #             verticalalignment='center',)
        #     plt.grid(False)
        #     plt.axis('off')

        output_names = ['% household', "% school", "% work", "% other", "% pediatric TB"]
        row = 0
        for output_name in output_names:
            row += 1

            # country name
            col = 0
            plt.subplot(gs[row, col])
            plt.text(x=-0.2, y=0.5, s=output_name, fontsize=5)
            plt.grid(False)
            plt.axis('off')

            j = 0
            #for param in ['proba_infection_per_contact', 'infectiousness_switching_age', 'n_colleagues']:
            for param in ['g_child', 'g_teen', 'g_adult']:


                j += 1
                # contacts heatmap
                col = j
                ax = plt.subplot(gs[row, col])
                self.make_a_sensitivity_plot(output_name, param, translate, param_vals, ax)

        filename = 'sensitivity_' + self.model_runner.country + '.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_a_sensitivity_plot(self, output_name, param, translate, param_vals, ax):

        x_vals = []
        y_vals = []

        for i_run in range(self.model_runner.nb_seeds):
            sample_id = self.model_runner.paths_to_calibrated_models[i_run].split("sample_")[1]
            sample_id = sample_id.split("_")[0]
            sample_id = int(sample_id)
            x_vals.append(float(param_vals[param][sample_id]))
            y_vals.append(self.get_output(i_run, output_name))

        plt.scatter(x_vals, y_vals, s=.6)
        plt.grid(False)

        ticks = [min(x_vals), max(x_vals)]
        ax.set_xticks(ticks)
        ax.set_xlim(ticks)

        show_x_ticks = 'off'
        if output_name == "% pediatric TB":
            show_x_ticks = 'on'
            plt.xlabel(translate[param],fontsize=5)

        show_y_ticks = 'off'
        if param == 'proba_infection_per_contact':
            show_y_ticks = 'on'

        plt.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=show_x_ticks,  # ticks along the bottom edge are off
                top='off',  # ticks along the top edge are off
                labelbottom=show_x_ticks,
                length=2.,
                pad=0.5,  # distance tick - label
                labelsize=4,
                colors='black'
            )  # labels along the bottom edge are off
        plt.tick_params(
                axis='y',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                left=show_y_ticks,  # ticks along the bottom edge are off
                right='off',  # ticks along the top edge are off
                labelleft=show_y_ticks,
                length=2.,
                pad=0.5,
                labelsize=4,
                colors='black'
            )  # labels along the bottom edge are off

    def get_output(self, i_run, output_name):

        if "diatric" not in output_name:
            contrib = self.model_runner.model_diagnostics['scenario_1']['n_contacts']['transmission_end_tb']
            num = 0.
            denom = 0.
            for key in contrib.keys():
                val = float(contrib[key][i_run])
                denom += val
                if "house" in output_name and key == 'household':
                    num = val
                elif "school" in output_name and key == 'school':
                    num = val
                elif "work" in output_name and key == 'workplace':
                    num = val
                elif "other" in output_name and key == 'community':
                    num = val
            if denom>0.:
                y = 100.*num/denom
            else:
                y = 0.
        else:
            tb_ages = self.model_runner.model_diagnostics['scenario_1']['all_tb_ages'][i_run]
            n_ped = len([a for a in tb_ages if a <= 15])
            y = float(n_ped) / float(len(tb_ages))

        return y


    def create_directories(self):
        """
        Create the directories where output files will be stored
        """
        # create the main directory using the name entered through the console
        dir_path = os.path.join(self.base_path, self.model_runner.data.console['project_name'])
        if not os.path.exists(dir_path):  # normally irrelevant. Must have been created during model-runner initialisation
            os.makedirs(dir_path)
        # create scenario-specific folders
        for scenario in self.model_runner.data.scenario_names:
            sc_dir_path = os.path.join(dir_path, scenario)
            if not os.path.exists(sc_dir_path):
                os.makedirs(sc_dir_path)
        # create a folder for the outputs displaying all scenarios

        all_dir_path = os.path.join(dir_path, 'all_scenarios')
        if not os.path.exists(all_dir_path):
            os.makedirs(all_dir_path)
        self.base_path = os.path.join(self.base_path, self.model_runner.data.console['project_name'])

    def fetch_timeseries_to_plot(self):
        for param in self.model_runner.data.console.keys():
            if 'plot_ts_' in param:
                if self.model_runner.data.console[param]:
                    self.timeseries_to_plot.append(param[8:])

    def convert_model_time_to_dates(self, times):
        """
        times contain the model times in days
        :param times:
        :return:
        """
        time_pyramid = self.model_runner.data.console['duration_burning_demo'] + self.model_runner.data.console[
            'duration_burning_tb']
        if self.model_runner.data.console['running_mode'] != 'manual':
            time_pyramid = self.model_runner.data.console['duration_burning_demo'] + self.model_runner.data.console[
                'duration_burning_tb'] + \
                           self.model_runner.data.console['n_years']
        remaining_years = [time_pyramid - g/365.25 for g in times]
        time_in_years = [self.model_runner.data.console['current_year'] - g for g in remaining_years]
        return time_in_years

    def make_graphs(self):
        """
        Generate the different graphs
        """
        if self.model_runner.data.console['running_mode'] == 'run_lhs_calibration':
            self.plot_lhs_calibration_results()
        self.make_scenario_specific_graphs()
        self.make_comparative_graphs()

        print str(self.i_figure) + " figures have been created."

    def make_scenario_specific_graphs(self):
        for scenario in self.model_runner.data.scenario_names:
            if self.model_runner.model_diagnostics[scenario]['n_accepted_runs'] == 0:
                continue

            if len(self.model_runner.data.console['years_checkpoints']) > 0:
                self.make_checkpoint_graphs_by_scenario(scenario)
            self.make_timeseries_graphs_by_scenario(scenario)
            if self.model_runner.data.console['plot_contact_heatmap'] or self.model_runner.data.console['load_calibrated_models']:
                self.plot_mixing_heatmaps(scenario, key='contact', interp='nearest')
            if self.model_runner.data.common_parameters['transmission']:
                self.plot_mixing_heatmaps(scenario, key='transmission', interp='nearest')
                self.plot_mixing_heatmaps(scenario, key='transmission_end_tb', interp='nearest')
            if self.model_runner.data.console['plot_all_tb_ages']:
                self.plot_all_tb_ages(scenario)
            if self.model_runner.data.console['plot_prem_like']:
                if self.model_runner.data.console['plot_contact_heatmap'] and len(self.model_runner.data.console['years_checkpoints']) > 0:
                    self.plot_prem_like(scenario)
                else:
                    print 'plot_prem_like requested but not possible.'
                    print 'missing heatmaps and/or age-pyramid'
            if self.model_runner.data.console['plot_transmission_by_location']:
                for key in ['contact', 'transmission', 'transmission_end_tb']:
                    self.plot_contribution_by_location(scenario, key_of_interest=key)

            if self.model_runner.data.console['running_mode'] == 'run_ks_based_calibration':
                self.make_calibration_graph(scenario, in_scenario_dir=True)

            if self.model_runner.data.console['plot_scale_up_functions']:
                self.plot_scale_up_functions(scenario)

            self.write_prevalence_by_age(scenario)

    def make_timeseries_graphs_by_scenario(self, scenario):
        """
        Generate all the timeseries graphs relevant to the scenario "scenario"
        """
        for series_name in self.timeseries_to_plot:
            self.plot_a_timeseries_variable(scenario, series_name)
            self.plot_an_aggregated_timeseries_variable(scenario, series_name)

    def plot_a_timeseries_variable(self, scenario, series_name):
        """
        Generate a spaghetti plot representing the timeseries obtained for the variable "series_name" and for the
        different runs of scenario "scenario"
        """
        data = self.model_runner.model_diagnostics[scenario]['timeseries'][series_name]

        # ymax = np.amax(data)
        self.i_figure += 1
        plt.figure(self.i_figure)
        for i in range(self.model_runner.data.console['n_runs'] * self.model_runner.nb_seeds):
            x = self.model_runner.model_diagnostics[scenario]['timeseries']['times'][i, ]
            x = self.convert_model_time_to_dates(x)
            y = data[i, ]
            plt.plot(x, y, color='black')
            x_max = math.ceil(x[-1] / 10.) * 10.
            plt.xlim((self.model_runner.data.console['start_plotting_time'], x_max))

        plt.xlabel('time (years)')
        plt.ylabel(self.dictionary[series_name])
        filename = scenario + '_' + series_name + '.' + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.base_path, scenario, filename)
        plt.savefig(file_path)
        plt.close()

    def plot_an_aggregated_timeseries_variable(self, scenario, series_name, ci=True):
        """
        Generate a plot representing the timeseries obtained for the variable "series_name" and the scenario "scenario".
        It draws a solid line representing the mean estimates and a shade for the asosciated confidence intervals.
        When ci is True, confidence intervals will be represented with shades around the average trend.
        """
        self.i_figure += 1
        plt.figure(self.i_figure)

        if scenario == 'All':
            scenarios = self.model_runner.data.scenario_names
        else:
            scenarios = [scenario]


        extra_delay = 0.
        if series_name == 'tb_incidence' and self.model_runner.data.console[
            'duration_burning_tb'] == 0:  # we want to ignore the first years of tb_incidence
            extra_delay = min(5., float(self.model_runner.data.console['n_years']))

        for i, sc in enumerate(scenarios):
            sc_name = self.model_runner.data.scenarios[sc]['scenario_title']
            data = self.model_runner.model_diagnostics[sc]['aggr_timeseries'][series_name]
            x = self.model_runner.model_diagnostics[sc]['timeseries']['times'][0, ]
            x = self.convert_model_time_to_dates(x)

            if i < len(self.scenario_colors):
                sc_color = self.scenario_colors[i]
            else:
                sc_color = 'b'

            if ci:
                y_low = data['low']
                y_high = data['high']
                y_mean = data['mean']
                plt.fill_between(x, y_low, y_high, color=sc_color, alpha=0.3, linewidth=0.)
            plt.plot(x, y_mean, color=sc_color, linewidth=2., label=sc_name)

        x_max = math.ceil(x[-1] / 10.) * 10.
        plt.xlim((self.model_runner.data.console['start_plotting_time'], x_max))
        plt.xlabel('time (years)')
        plt.ylabel(self.dictionary[series_name])

        filename = scenario + '_' + series_name + '_aggr.' + self.model_runner.data.console['graph_format']
        if scenario == 'All':
            plt.legend(loc=2)
            file_path = os.path.join(self.base_path, 'all_scenarios', filename)
        else:
            file_path = os.path.join(self.base_path, scenario, filename)

        plt.savefig(file_path)
        plt.close()

    def plot_mixing_heatmaps(self, scenario, key='contact', interp='nearest'):
        """
        key is one of 'contact', 'transmission', 'transmission_end_tb'
        """
        # create a new directory
        dir_path = os.path.join(self.base_path, scenario, key)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        # create a matrix with the sum of all contacts
        self.model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'] = np.zeros((101, 101))

        for location in self.model_runner.model_diagnostics[scenario]['contact_matrices'][key].keys():
            self.model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'] +=\
                self.model_runner.model_diagnostics[scenario]['contact_matrices'][key][location]

        # fix extreme values for the color range
        vmin = 0.
        vmax_all = np.amax(self.model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'])

        for location in self.model_runner.model_diagnostics[scenario]['contact_matrices'][key].keys():
            vmax = np.amax(self.model_runner.model_diagnostics[scenario]['contact_matrices'][key][location])
            self.i_figure += 1
            with plt.style.context(('seaborn-dark')):
                plt.figure(self.i_figure)

                a = self.model_runner.model_diagnostics[scenario]['contact_matrices'][key][location]
                plt.imshow(a.transpose(), cmap='hot', interpolation=interp, origin='lower', vmin=vmin, vmax=vmax)

                plt.xlabel('index age')
                plt.ylabel('contact age')
                # filename = scenario + '_mixing_pattern_' + key + '_' + location + '.' + self.model_runner.data.console['graph_format']  # makes the program crash
                filename = 'mixing_pattern_' + key + '_' + location + '.' + self.model_runner.data.console['graph_format']

                file_path = os.path.join(dir_path, filename)
                plt.savefig(file_path)
                plt.close()

    def get_contact_proportions(self, scenario, key='contact', age_min=0, age_max=100):

        # create a matrix with the sum of all contacts
        self.model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'] = np.zeros((101, 101))

        for location in self.model_runner.model_diagnostics[scenario]['contact_matrices'][key].keys():
            self.model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'] += \
                self.model_runner.model_diagnostics[scenario]['contact_matrices'][key][location]

        perc_contacts = {}
        for location in self.model_runner.model_diagnostics[scenario]['contact_matrices'][key].keys():
            a = self.model_runner.model_diagnostics[scenario]['contact_matrices'][key][location]
            index_min = max(0, age_min)
            index_max = min(101, age_max)
            a = a[index_min:index_max, ]
            perc_contacts[location] = np.sum(a)

        denominator = copy.copy(perc_contacts['all'])
        for location in self.model_runner.model_diagnostics[scenario]['contact_matrices'][key].keys():
            perc_contacts[location] = 100.*perc_contacts[location]/denominator

        print perc_contacts
        return perc_contacts

    def plot_all_tb_ages(self, scenario):
        self.i_figure += 1
        plt.figure(self.i_figure)

        x = [age for sublist in self.model_runner.model_diagnostics[scenario]['all_tb_ages'] for age in sublist]

        bottom_bounds = [0., 5., 15., 25., 35., 45., 55., 65.]
        values = [0 for i in range(len(bottom_bounds))]
        labs = ['' for i in range(len(bottom_bounds))]
        height = [0. for i in range(len(bottom_bounds))]
        for j, bottom_bound in enumerate(bottom_bounds):
            if bottom_bound == 0.:
                width = 5.
                height[j] = 5.
                labs[j] = '0-4'
            elif bottom_bound == 65.:
                width = 50.
                height[j] = 10.
                labs[j] = '65+'
            else:
                width = 10.
                height[j] = 10.
                labs[j] = str(int(bottom_bound)) + '-' + str(int(bottom_bound + 9))
            values[j] = sum(1 for age in x if bottom_bound <= age < (bottom_bound + width))

        center_ticks = [bottom_bounds[j] + 0.5*height[j] for j in range(len(bottom_bounds))]

        plt.barh(y=center_ticks, width=values, height=height, tick_label=labs, align='center')
        plt.tick_params(length=0.)
        ages_under_15 = [age for age in x if age < 15.]

        if len(x) > 0:
            prop_under_15 = float(len(ages_under_15))/ float(len(x))
            print "Country: " + self.model_runner.country
            print "Proportion of pediatric TB (<15 yo): " + str(prop_under_15)

        plt.xlabel('number of TB cases')
        plt.ylabel('age (years)')
        filename = scenario + '_all_tb_ages.' + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.base_path, scenario, filename)
        plt.savefig(file_path)
        plt.close()

    def make_checkpoint_graphs_by_scenario(self, scenario):
        """
        Generate all the checkpoint graphs relevant to the scenario "scenario"
        """
        self.plot_age_and_household_distributions(scenario)

    def plot_age_and_household_distributions(self, scenario, country=None, save=True, keys=None):
        """
        Plot histograms for the aggregated age distribution of the runs relevant to scenario "scenario"
        """
        if keys is None:
            data_keys = ['ages', 'household_sizes', 'school_sizes', 'workplace_sizes']
        else:
            data_keys = keys
        x_labs = {'ages': 'age distribution', 'household_sizes': 'household size',
                  'school_sizes': 'school size', 'workplace_sizes': 'workplace size'}
        bins = {'ages': range(100), 'household_sizes': range(20), 'school_sizes': 100, 'workplace_sizes': 100}
        orientation = {'ages': 'horizontal', 'household_sizes': 'vertical',
                       'school_sizes': 'vertical' , 'workplace_sizes': 'vertical'}

        for data_key in data_keys:
            for checkpoint in self.model_runner.model_diagnostics[scenario]['checkpoint_outcomes'][data_key].keys():
                self.i_figure += 1
                plt.figure(self.i_figure)
                values = [value for i_run in range(self.model_runner.data.console['n_runs'] * self.model_runner.nb_seeds)
                        for value in self.model_runner.model_diagnostics[scenario]['checkpoint_outcomes'][data_key][checkpoint][i_run]]

                if len(values) == 0:
                    print "No values to be plotted for " + data_key + " at year " + str(int(round(checkpoint/365.25)))
                    continue

                plt.hist(values, bins=bins[data_key], orientation=orientation[data_key], density=True)

                plt.xlabel(x_labs[data_key])
                if data_key == 'household_sizes':
                    plt.xlim((0., 18.))
                if country is not None:
                    plt.title(country, fontsize=40, y=1.03)
                    plt.ylabel('proportion', fontsize=35, color='black')
                    plt.xlabel(x_labs[data_key], fontsize=35, color='black')
                    plt.tick_params(axis='both', which='major', labelsize=30, labelcolor='black')
                    plt.ylim((0., 0.35))
                    x_tick_pos = np.linspace(1.5, 17.5, num=9)
                    x_tick_labs = np.linspace(1, 17, num=9)
                    x_tick_labs = [int(u) for u in x_tick_labs]

                    plt.xticks(x_tick_pos, x_tick_labs, fontsize=30, color='black')
                if save:
                    filename = scenario + '_' + data_key + '_at_year_' + str(round(checkpoint/365.25)) + '.' + self.model_runner.data.console['graph_format']
                    file_path = os.path.join(self.base_path, scenario, filename)
                    plt.savefig(file_path)
                    plt.close()

    def plot_prem_like(self, scenario):
        self.i_figure += 1
        plt.figure(self.i_figure, figsize=(6, 9))
        x_labs = {'ages': 'age distribution', 'household_sizes': 'household size distribution'}
        bins = {'ages': range(100), 'household_sizes': range(20)}
        orientation = {'ages': 'horizontal', 'household_sizes': 'vertical'}

        # Age pyramid and (Household sizes)
        for i, key in enumerate(['ages']): #, 'household_sizes']):
            plt.subplot(3, 2, i+1)
            checkpoint = self.model_runner.model_diagnostics[scenario]['checkpoint_outcomes'][key].keys()[-1]
            values = [value for i_run in range(self.model_runner.data.console['n_runs'] * self.model_runner.nb_seeds)
                      for value in
                      self.model_runner.model_diagnostics[scenario]['checkpoint_outcomes'][key][checkpoint][i_run]]

            plt.hist(values, bins=bins[key], orientation=orientation[key])
            plt.xlabel(x_labs[key])

        # create a matrix with the sum of all contacts
        self.model_runner.model_diagnostics[scenario]['contact_matrices']['contact']['all'] = np.zeros((101, 101))
        for location in self.model_runner.model_diagnostics[scenario]['contact_matrices']['contact'].keys():
            if location != 'all':
                self.model_runner.model_diagnostics[scenario]['contact_matrices']['contact']['all'] += \
                    self.model_runner.model_diagnostics[scenario]['contact_matrices']['contact'][location]

        # Heatmap
        # fix extreme values for the color range
        vmin = 0.
        vmax_all = np.amax(self.model_runner.model_diagnostics[scenario]['contact_matrices']['contact']['all'])

        contact_types = ["household", "workplace", "school", "community", "all"]
        titles = ["Home", "Work", "School", "Other locations", "All locations"]

        for i, location in enumerate(contact_types):
            vmax = np.amax(self.model_runner.model_diagnostics[scenario]['contact_matrices']['contact'][location])
            with plt.style.context(('seaborn-dark')):
                plt.subplot(3, 2, i + 2)
                a = self.model_runner.model_diagnostics[scenario]['contact_matrices']['contact'][location]
                plt.imshow(a.transpose(), cmap='hot', interpolation='nearest', origin='lower', vmin=vmin, vmax=vmax)

                plt.text(50., 90., titles[i], horizontalalignment='center', color='white')
                plt.xlabel('index age')
                plt.ylabel('contact age')

        filename = 'Prem_like_figure_' + scenario + '_' + self.model_runner.data.console['country'] + '.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.base_path, scenario, filename)
        plt.savefig(file_path)
        plt.close()

    def plot_contribution_by_location(self, scenario, key_of_interest):
        """
        key_of_interest is one of "contact", "transmission", "transmission_end_tb"
        """
        self.i_figure += 1
        plt.figure(self.i_figure)
        locations = ['home', 'school', 'work', 'other']
        keys = ['household', 'school', 'workplace', 'community']
        x_tick_pos = np.arange(1, 5)
        width = 0.5

        ylabs = {'contact': 'contribution to overall social mixing (%)',
                 'transmission': 'contribution to overall transmission (%)',
                 'transmission_end_tb': 'contribution to overall TB burden(%)'}

        perc_dict = copy.deepcopy(self.model_runner.model_diagnostics[scenario]['n_contacts'][key_of_interest])

        for i_run in range(self.model_runner.data.console['n_runs'] * self.model_runner.nb_seeds):
            sum_trans = 0
            for key in keys:
                sum_trans += perc_dict[key][i_run]
            if sum_trans > 0:
                for key in keys:
                 perc_dict[key][i_run] = 100.*float(perc_dict[key][i_run])/float(sum_trans)
        means = []
        stds = []
        for key in keys:
            means.append(np.mean(perc_dict[key]))
            if self.model_runner.data.console['n_runs'] > 1 or self.model_runner.nb_seeds > 1:
                stds.append(np.std(perc_dict[key]))
        if self.model_runner.data.console['n_runs'] == 1 and self.model_runner.nb_seeds == 1:
            stds = None

        plt.bar(x_tick_pos, means, width, yerr=stds, ecolor='black', align='center')

        plt.xlabel('location of contact')
        plt.ylabel(ylabs[key_of_interest])
        plt.xticks(x_tick_pos, locations)
        plt.xlim(0.5, 4.5)
        plt.ylim(0., 100.)

        filename = 'contributions_' + key_of_interest + '.' + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.base_path, scenario, filename)
        plt.savefig(file_path)
        plt.close()

    def plot_scale_up_functions(self, scenario):
        datasets = {'bcg_coverage_prop': self.model_runner.data.data_from_sheets['bcg'], 'treatment_success_prop':
            self.model_runner.data.data_from_sheets['outcomes']['c_new_tsr'], 'cdr_prop': self.model_runner.data.data_from_sheets['gtb_2016']['c_cdr']}

        # smear-positive proportion separately
        sp_number = self.model_runner.data.data_from_sheets['outcomes']['new_sp_coh']
        snep_number = self.model_runner.data.data_from_sheets['outcomes']['new_snep_coh']
        dataset_for_prop = {}
        for year in snep_number.keys():
            if year in sp_number.keys() and (snep_number[year] + sp_number[year]) > 0.:
                dataset_for_prop[year] = 100. * sp_number[year] / (sp_number[year] + snep_number[year])

        datasets.update({'sp_prop': dataset_for_prop})

        for key, func in self.model_runner.data.scale_up_functions.iteritems():
            self.i_figure += 1
            plt.figure(self.i_figure)

            time_in_years = np.linspace(1950., 2020., 1000.)
            y = [100.*func(val) for val in time_in_years]
            plt.plot(time_in_years, y, color='black')

            for year, value in datasets[key].iteritems():
                plt.plot(year, value, 'ro', markersize=3.)

            plt.xlim((1950., 2020.))
            plt.ylim((0., 100.))
            plt.xlabel('time (years)')
            plt.ylabel(self.dictionary[key])
            filename = scenario + '_scaleup_' + key + '.' + self.model_runner.data.console['graph_format']
            file_path = os.path.join(self.base_path, scenario, filename)
            plt.savefig(file_path)
            plt.close()

    def make_comparative_graphs(self):
        self.make_comparative_checkpoint_graphs()
        self.make_comparative_timeseries_graphs()
        if self.model_runner.data.console['running_mode'] == 'run_ks_based_calibration':
            self.make_calibration_graph()

    def make_comparative_checkpoint_graphs(self):
        pass

    def make_comparative_timeseries_graphs(self):
        """
        Create a graph representing a comparison between the different scenarios for all timeseries
        """
        for series_name in self.timeseries_to_plot:
            self.plot_an_aggregated_timeseries_variable(scenario='All', series_name=series_name, ci=True)

    def make_calibration_graph(self, scenario_to_plot=None, in_scenario_dir=False):
        if os.name != 'nt':  # seaborn package not installed on remote server
            return

        if scenario_to_plot is None:
            scenarios_to_consider = self.model_runner.data.scenario_names
        else:
            scenarios_to_consider = [scenario_to_plot]

        self.i_figure += 1
        plt.figure(self.i_figure)

        mu = self.model_runner.data.common_parameters['targetted_prevalence']
        width = self.model_runner.data.common_parameters['targetted_prevalence_high'] - self.model_runner.data.common_parameters['targetted_prevalence_low']
        sd = width / (2. * 1.96)
        random_samples = [np.random.normal(mu, sd, 100000) for _ in range(len(scenarios_to_consider))]
        positions = [float(i+1) for i in range(len(scenarios_to_consider))]
        labels = [scenario.split('_')[-1] for scenario in scenarios_to_consider]
        violin_widths = [0.6] * len(scenarios_to_consider)
        if scenario_to_plot is not None:
            violin_widths = [0.4]
        plt.violinplot(random_samples, showmeans=False, showmedians=True, showextrema=False, positions=positions,
                       widths=violin_widths)
        best_dist = 1.e10
        for i, scenario in enumerate(scenarios_to_consider):
            x = [float(i+1)] * self.model_runner.data.console['n_runs'] * self.model_runner.nb_seeds
            model_outputs_for_calibration = []
            target_indicator = 'tb_prevalence'  # hard-coded
            for i_run in range(self.model_runner.data.console['n_runs'] * self.model_runner.nb_seeds):
                model_output = self.model_runner.model_diagnostics[scenario]['timeseries'][target_indicator][i_run, -1]
                model_outputs_for_calibration.append(model_output)
            plt.plot(x, model_outputs_for_calibration, 'ro', markersize=3.)

            dist, p_value = stats.ks_2samp(model_outputs_for_calibration, random_samples[i])
            dist = round(dist, 4)
            plt.text(x=x[0]-0.4, y=self.model_runner.data.common_parameters['targetted_prevalence_high'],  s='D=' + str(dist),
                     rotation=90, verticalalignment='bottom')

            if dist < best_dist:
                best_dist = dist
                best_scenario = scenario

        plt.xlabel(self.dictionary[self.model_runner.data.calibration_params.keys()[0]])
        plt.ylabel(self.dictionary['tb_prevalence'])

        plt.xticks(positions, labels)

        filename = 'calibration_results.' + self.model_runner.data.console['graph_format']
        if scenario_to_plot is not None:
            plt.xlim(0., 2.)
            filename = 'calibration_results_param_' + scenario_to_plot.split('_')[-1] + '.' + self.model_runner.data.console['graph_format']
        if in_scenario_dir:
            file_path = os.path.join(self.base_path, scenario_to_plot, filename)
        else:
            file_path = os.path.join(self.base_path, 'all_scenarios', filename)
        plt.savefig(file_path)
        plt.close()

        if scenario_to_plot is None:
            self.make_calibration_graph(best_scenario)

    def perform_anova_test_comparing_seeds(self, scenario):
        all_model_outputs = {}
        target_indicator = 'tb_prevalence'  # hard-coded
        for seed_index in range(self.model_runner.nb_seeds):
            all_model_outputs[seed_index] = []
            for run_index in range(self.model_runner.data.console['n_runs']):
                i_run = seed_index * self.model_runner.data.console['n_runs'] + run_index
                model_output = self.model_runner.model_diagnostics[scenario]['timeseries'][target_indicator][i_run, -1]
                all_model_outputs[seed_index].append(model_output)

        # anova_results = stats.f_oneway()  # arguments have to be listed explicitly
        # print anova_results

    def write_prevalence_by_age(self, scenario):
        age_breaks = [0., 5., 10., 15., 25., 35., 45., 55., 65.]
        age_cats = [['X_1'], ['X_2'], ['X_3'], ['X_4', 'X_5'], ['X_6', 'X_7'], ['X_8', 'X_9'], ['X_10', 'X_11'],
                    ['X_12', 'X_13'], ['X_14', 'X_15', 'X_16', 'X_17']]
        prev_by_age = self.model_runner.model_diagnostics[scenario]['tb_prevalence_by_age']
        nb_runs = len(prev_by_age)

        for i, age_min in enumerate(age_breaks):
            prev = 0.
            prop_pop_in_agegroup = 0.
            for age_cat in age_cats[i]:
                prop_pop_in_agegroup += self.model_runner.data.age_pyramid[age_cat]
            nb_ind_in_agegroup = prop_pop_in_agegroup * self.model_runner.data.common_parameters['population'] # * nb_runs

            for j in range(nb_runs):
                prev += prev_by_age[j][i]
            prev /= nb_runs
            prev_as_prop = prev / 1.e5
            unc_gap = 1.96 * math.sqrt(prev_as_prop*(1. - prev_as_prop)/nb_ind_in_agegroup)
            unc_gap *= 1.e5
            print "Prevalence for age-group starting at " + str(age_min) + ": " + str(prev) + "( " + str(prev - unc_gap) +\
                " - " + str(prev + unc_gap) + ")"

    def write_prevalence_by_age_new(self, scenario):
        age_breaks = [0., 5., 10., 15., 25., 35., 45., 55., 65.]

        prev_by_age = self.model_runner.model_diagnostics[scenario]['tb_prevalence_by_age']
        for i, age_min in enumerate(age_breaks):
            prevs = [prev_by_age[j][i] for j in range(len(prev_by_age))]

            prevs_mean = np.mean(prevs)
            prev_sd = np.std(prevs)
            prevs_low = prevs_mean - 1.96*prev_sd/np.sqrt(float(len(prev_by_age)))
            prevs_high = prevs_mean + 1.96*prev_sd/np.sqrt(float(len(prev_by_age)))

            print "Prevalence for age-group starting at " + str(age_min) + ": " + str(prevs_mean) + " (" + str(prevs_low) +\
                " - " + str(prevs_high) + ")"

    def plot_lhs_calibration_results(self):
        self.i_figure += 1
        plt.figure(self.i_figure)

        for i, scenario in enumerate(self.model_runner.data.scenario_names):
            plt.plot(i, self.model_runner.model_diagnostics[scenario]['n_accepted_runs'], 'ro', markersize=3.)

        plt.xlabel('lhs sample #')
        plt.ylabel('n accepted')
        file_path = os.path.join(self.base_path, 'all_scenarios', 'lhs_outputs.' +
                                 self.model_runner.data.console['graph_format'])
        plt.savefig(file_path)
        plt.close()

        for param in self.model_runner.data.sampled_params.keys():
            self.i_figure += 1
            plt.figure(self.i_figure)
            for i, scenario in enumerate(self.model_runner.data.scenario_names):
                param_val = self.model_runner.data.sampled_params[param][i]
                n_accepted = self.model_runner.model_diagnostics[scenario]['n_accepted_runs']
                plt.plot(param_val, n_accepted, 'ro', markersize=3.)

            plt.xlabel(param)
            plt.ylabel('n accepted')

            file_path = os.path.join(self.base_path, 'all_scenarios', 'lhs_outputs_' + param + '.' +
                                     self.model_runner.data.console['graph_format'])
            plt.savefig(file_path)
            plt.close()

    def write_timeseries_to_csv(self):
        start_year = self.model_runner.data.console['duration_burning_demo'] + self.model_runner.data.console['duration_burning_tb']

        for sc in self.model_runner.data.scenario_names:
            times = list(self.model_runner.model_diagnostics[sc]['aggr_timeseries']['times']['mean'])
            times = [round(t - start_year*365.25) for t in times]
            start_index = next(t[0] for t in enumerate(times) if t[1] >= 0.)
            times = times[start_index:]

            file_path = os.path.join(self.base_path, sc, 'timeseries.csv')

            with open(file_path, "w") as f:
                writer = csv.writer(f, lineterminator='\n')
                row1 = times
                row1.insert(0, 'times')
                writer.writerow(row1)
                for timeseries_name in self.model_runner.model_diagnostics[sc]['aggr_timeseries'].keys():
                    if timeseries_name == 'times':
                        continue
                    for type in ['mean', 'low', 'high']:
                        name = timeseries_name + '_' + type
                        row = list(self.model_runner.model_diagnostics[sc]['aggr_timeseries'][timeseries_name][type])
                        row = row[start_index:]
                        row.insert(0, name)
                        writer.writerow(row)

    def store_outputs(self):
        print "Pickling the outputs object before producing figures..."
        os.path.join(self.base_path, "pickled_outputs.pickle")
        file_name = os.path.join(self.base_path, "pickled_outputs.pickle")
        file_stream = open(file_name, "wb")

        saved_data = copy.deepcopy(self.model_runner.data)

        # We must delete some attributes to allow pickling
        del self.model_runner.data
        del self.model_runner.m_init
        dill.dump(self, file_stream)
        file_stream.close()

        self.model_runner.data = saved_data

        print "Complete."

def load_outputs(file_path):
    print "Loading outputs from " + file_path + " ..."
    file_stream = open(file_path, "rb")
    loaded_outputs = dill.load(file_stream)
    file_stream.close()
    # we must rebuild the data attribute of model_runner as it had been discarded for pickling purpose
    loaded_outputs.model_runner.data = imp.data(loaded_outputs.model_runner.country,
                                                loaded_outputs.model_runner.calibration_params,
                                                loaded_outputs.model_runner.uncertainty_params
                                                )

    #
    # prov = {None: {'proba_infection_per_contact': np.linspace(start=0.002, stop=0.0025, num=2)},
    #                     'India': {'proba_infection_per_contact': np.linspace(start=0.0020, stop=0.0026, num=4)},
    #                     'Indonesia': {'proba_infection_per_contact': np.linspace(start=0.0012, stop=0.0018, num=4)},
    #                     'China': {'proba_infection_per_contact': np.linspace(start=0.007, stop=0.009, num=3)},
    #                     'Philippines': {'proba_infection_per_contact': np.linspace(start=0.0016, stop=0.0022, num=4)},
    #                     'Pakistan': {'proba_infection_per_contact': np.linspace(start=0.0034, stop=0.0040, num=4)}
    #                   }
    # loaded_outputs.model_runner.data = imp.data(loaded_outputs.model_runner.country,
    #                                             prov['Philippines'])

    loaded_outputs.base_path = os.path.join('outputs')
    loaded_outputs.create_directories()
    loaded_outputs.i_figure = 0

    # We need to retrieve the number of runs automatically as it may differ form the one currently written in the
    # console spreadsheet.
    scenario = loaded_outputs.model_runner.data.scenarios.keys()[0]
    loaded_outputs.model_runner.data.console['n_runs'] = \
        loaded_outputs.model_runner.model_diagnostics[scenario]['timeseries']['birth_rate'].shape[0] /\
        loaded_outputs.model_runner.nb_seeds

    print "Complete."
    return loaded_outputs

class multi_output:
    def __init__(self, countries=None, project_name='multi'):
        """
        :param countries: list of countries to plot (they must have been run before)
        :param project_name: beginning of the directory name for each country (ex "multi" for
         "multi_Viet Nam", "multi_Philippines", ...)
        """
        self.countries = countries
        self.project_name = project_name
        self.i_figure = 0
        self.multi_base_path = os.path.join('outputs', 'multi_country_outputs')  # will be updated once directories are created
        self.create_multi_directories()
        self.output_objects = {}
        self.load_output_objects()

    def create_multi_directories(self):
        if not os.path.exists(self.multi_base_path):
            os.makedirs(self.multi_base_path)

    def load_output_objects(self):
        """
        Load the previously saved outputs objects for the different countries. Populate self.output_objects
        """
        for country in self.countries:
            dir_name = self.project_name + "_" + country
            file_path = os.path.join('outputs', dir_name, 'pickled_outputs.pickle')
            self.output_objects[country] = load_outputs(file_path)
            # print country
            # print len(self.output_objects[country].model_runner.model_diagnostics['scenario_1']['checkpoint_outcomes']['ages'][73000.0])

    def make_country_graphs(self):
        for country in self.countries:
            self.output_objects[country].i_figure = self.i_figure
            self.output_objects[country].make_graphs()
            self.i_figure = self.output_objects[country].i_figure

    def make_multi_graphs(self):
        self.make_heatmaps_figure()
        self.make_double_pyramids()
        self.make_contribution_graphs()
        self.make_tb_burden_clocks()
        # self.make_multi_household_sizes()
        self.make_contact_heatmaps_figure()
        # self.make_multi_scale_up_functions()

    def make_heatmaps_figure(self):
        """
        Figure 2 of manuscript
        """
        self.i_figure += 1
        nb_countries = len(self.countries)
        plt.figure(self.i_figure, figsize=(3.7, (nb_countries + 0.5)*0.9))

        h_ratios = [0.5]
        for _ in self.countries:
            h_ratios.append(5.)

        gs = gridspec.GridSpec(nb_countries + 1, 4,
                               width_ratios=[0.7, 1., 1., 1.],
                               height_ratios=h_ratios
                               )
        gs.update(wspace=0.12, hspace=0.1)

        #  panel titles
        titles = ['Contacts', 'Transmission', 'Transmission leading\nto active TB']
        row = 0
        for i, title in enumerate(titles):
            plt.subplot(gs[row, i+1])
            plt.text(x=0.5, y=0.9, s=title, fontsize=5, horizontalalignment='center',
                    verticalalignment='center',)
            plt.grid(False)
            plt.axis('off')

        for country in self.countries:
            row += 1
            x_axis = False
            if country == self.countries[-1]:
                x_axis = True

            # country name
            col = 0
            plt.subplot(gs[row, col])
            plt.text(x=-0.2, y=0.5, s=country, fontsize=5)
            plt.grid(False)
            plt.axis('off')

            # contacts heatmap
            col = 1
            plt.subplot(gs[row, col])
            self.make_a_heatmap(country, 'contact', 'all', show_x_axis=x_axis, show_y_axis=True)

            # transmission heatmap
            col = 2
            plt.subplot(gs[row, col])
            self.make_a_heatmap(country, 'transmission', 'all', show_x_axis=x_axis)

            # tb heatmap
            col = 3
            plt.subplot(gs[row, col])
            self.make_a_heatmap(country, 'transmission_end_tb', 'all', show_x_axis=x_axis)

        filename = 'multi_country_heatmaps.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_a_heatmap(self, country, key, location, show_x_axis=False, show_y_axis=False):
        scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
        self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'] = np.zeros((101, 101))

        for loc in self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key].keys():
            self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'] += \
                self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key][loc]

        # fix extreme values for the color range
        vmin = 0.
        vmax = np.amax(self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all'])

        vmax = np.amax(self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key][location])
        with plt.style.context(('seaborn-dark')):
            a = self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key][location]
            plt.imshow(a.transpose(), cmap='hot', interpolation='nearest', origin='lower', vmin=vmin, vmax=vmax)

            plt.grid(False)
            show_x_ticks = 'off'
            show_y_ticks = 'off'
            if show_y_axis:
                show_y_ticks = 'on'
                plt.ylabel('contact age', fontsize=5, labelpad=0., color='black')
                plt.yticks([0., 20., 40., 60., 80.], [0, 20, 40, 60, 80])
            else:
                plt.ylabel('')
            if show_x_axis:
                show_x_ticks = 'on'
                plt.xlabel('index age', fontsize=5, labelpad=0., color='black')
                plt.xticks([0., 20., 40., 60., 80.], [0, 20, 40, 60, 80])
            else:
                plt.xlabel('')

            plt.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=show_x_ticks,  # ticks along the bottom edge are off
                top='off',  # ticks along the top edge are off
                labelbottom=show_x_ticks,
                length=2.,
                pad=0.5,  # distance tick - label
                labelsize=4,
                colors='black'
            )  # labels along the bottom edge are off
            plt.tick_params(
                axis='y',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                left=show_y_ticks,  # ticks along the bottom edge are off
                right='off',  # ticks along the top edge are off
                labelleft=show_y_ticks,
                length=2.,
                pad=0.5,
                labelsize=4,
                colors='black'
            )  # labels along the bottom edge are off


        # contribution of age-categories
        print "###############################"
        print key + ' in ' + self.output_objects[country].model_runner.country
        a = self.output_objects[country].model_runner.model_diagnostics[scenario]['contact_matrices'][key]['all']
        all_contacts = a.sum()
        S = 0
        # for i in range(17):
        #     age_min = i*5
        #     age_max = age_min + 5
        #     recipients = a[:, age_min:age_max].sum()
        #     index = a[age_min:age_max, :].sum()
        #     involving_either = recipients + index - a[age_min:age_max, age_min:age_max].sum()
        #
        #     if all_contacts > 0:
        #         S += index / all_contacts
        #         print "Proportion of all events involving " + str(age_min) + "-" + str(age_max) + " as index: " + str(index / all_contacts)
        #         # print "Proportion of contacts involving 15-20 as either recipient or index: " + str(
        #         #     contacts_involving_15_20 / all_contacts)

        contacts_recipients_15_20 = a[:, 15:21].sum()
        contacts_index_15_20 = a[15:21, :].sum()
        contacts_involving_15_20 = contacts_recipients_15_20 + contacts_index_15_20 - a[15:21, 15:21].sum()
        contacts_involving_15_20_only = a[15:21, 15:21].sum()

        #
        if all_contacts > 0 and key != 'contact':

            print "Perc of contacts involving 15-20 as index: " + str(round(100.*contacts_index_15_20 / all_contacts))
            print "Perc of contacts involving 15-20 as recipient: " + str(
                round(100.*contacts_recipients_15_20 / all_contacts))
            print "Perc of contacts between 15-20 : " + str(
                round(100.*contacts_involving_15_20_only / all_contacts))
            print "Perc of contacts involving 15-20 as either recipient or index: " + str(
                round(100.*contacts_involving_15_20 / all_contacts))

        print "S = " + str(S)
        # contribution of 15-20 years old
        # print "#########"
        # print key
        #
        # contacts_recipients_15_20 = a[:, 15:21].sum()
        # contacts_index_15_20 = a[15:21, :].sum()
        # contacts_involving_15_20 = contacts_recipients_15_20 + contacts_index_15_20 - a[15:21, 15:21].sum()
        #
        # if all_contacts > 0:
        #     print "Proportion of contacts involving 15-20 as recipients: " + str(contacts_recipients_15_20 / all_contacts)
        #     print "Proportion of contacts involving 15-20 as index: " + str(contacts_index_15_20 / all_contacts)
        #     print "Proportion of contacts involving 15-20 as either recipient or index: " + str(contacts_involving_15_20/all_contacts)


        # 10-14   20-24

    def make_contact_heatmaps_figure(self):
        self.i_figure += 1
        nb_countries = len(self.countries)
        plt.figure(self.i_figure, figsize=(5.7, (nb_countries + 0.5) * 0.9))

        h_ratios = [0.5]
        for _ in self.countries:
            h_ratios.append(5.)

        gs = gridspec.GridSpec(nb_countries + 1, 6,
                               width_ratios=[0.7, 1., 1., 1., 1., 1.],
                               height_ratios=h_ratios
                               )
        gs.update(wspace=0.0, hspace=0., left=0.25, right=0.9, bottom=0.2)

        #  panel titles
        locations = ['household', 'school', 'workplace', 'other', 'all']
        row = 0
        for i, location in enumerate(locations):
            plt.subplot(gs[row, i + 1])
            plt.text(x=0.5, y=0.9, s=location, fontsize=5, horizontalalignment='center',
                     verticalalignment='center', )
            plt.grid(False)
            plt.axis('off')

        for country in self.countries:
            row += 1
            # country name
            plt.subplot(gs[row, 0])
            plt.text(x=-0.4, y=0.5, s=country, fontsize=5)
            plt.grid(False)
            plt.axis('off')

            # contacts heatmaps
            for i, location in enumerate(['household', 'school', 'workplace', 'community', 'all']):
                plt.subplot(gs[row, i+1])
                y_axis = False
                if i == 0:
                    y_axis = True
                x_axis = False
                if row == len(self.countries):
                    x_axis = True
                self.make_a_heatmap(country, 'contact', location, show_x_axis=x_axis, show_y_axis=y_axis)

        filename = 'multi_country_contact_heatmaps.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_double_pyramids(self):
        self.i_figure += 1
        plt.figure(self.i_figure, figsize=(40., 20.))
        gs = gridspec.GridSpec(2, 3)
        gs.update(hspace=0.3)

        colors = {'left': (132./255., 186./255., 91./255.), 'right': (211./255., 94./255., 96./255.)}

        for i, country in enumerate(self.countries):
            col = i%3
            row = 0
            if i > 2:
                row = 1
            ax = plt.subplot(gs[row, col])
            self.make_age_tb_age_pyramid(ax, country, colors)

            age_min, age_max = 15., 24.
            self.get_age_specific_incidence(country, age_min, age_max)
            if country == 'Philippines':
                self.get_age_specific_incidence(country, 10., 19.)

        # legend
        ax = plt.subplot(gs[1, 2])
        ax.add_patch(patches.Rectangle(xy=(0., 0.45), width=0.16, height=0.1, color=colors['left']))
        ax.add_patch(patches.Rectangle(xy=(0., 0.25), width=0.16, height=0.1, color=colors['right']))

        ax.text(x=0.21, y=0.45, s='Population age pyramid', fontsize=40, verticalalignment='bottom')
        ax.text(x=0.21, y=0.25, s='Age distribution of TB cases', fontsize=40, verticalalignment='bottom')
        ax.axis('off')

        filename = 'multi_country_pyramids.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_age_tb_age_pyramid(self, ax, country, colors):
        scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
        if len(self.output_objects[country].model_runner.data.console['years_checkpoints']) == 0:
            print "No checkpoint found. tb_age_pyramid could not be generated."
            return
        checkpoint = self.output_objects[country].model_runner.model_diagnostics[scenario]['checkpoint_outcomes']['ages'].keys()[-1]
        print len(self.output_objects[country].model_runner.model_diagnostics[scenario]['checkpoint_outcomes']['ages'][checkpoint])
        x_left = [value for i_run in range(self.output_objects[country].model_runner.data.console['n_runs']  * self.output_objects[country].model_runner.nb_seeds)
                  for value in
                  self.output_objects[country].model_runner.model_diagnostics[scenario]['checkpoint_outcomes']['ages'][checkpoint][i_run]]

        tb_ages = [age for sublist in self.output_objects[country].model_runner.model_diagnostics[scenario]['all_tb_ages'] for age in sublist]

        bins = [5.*x for x in range(21)]
        weights_left = [-1. for _ in range(len(x_left))]

        if len(tb_ages) == 0:
            return

        multiplier = len(x_left) / len(tb_ages)
        weights_right = [multiplier for _ in range(len(tb_ages))]

        right = ax.hist(tb_ages, orientation='horizontal', weights=weights_right, bins=bins, color=(1., 1., 1., 0.))  # hidden
        left = ax.hist(x_left, orientation='horizontal', weights=weights_left, bins=bins, color=colors['left'])

        plt.title(country, fontsize=40, y=1.03)

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom='on',  # ticks along the bottom edge are off
            top='off',  # ticks along the top edge are off
            labelbottom='on',
            length=8.)  # labels along the bottom edge are off
        plt.tick_params(
            axis='y',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            left='off',  # ticks along the bottom edge are off
            right='off',  # ticks along the top edge are off
            labelleft='off')  # labels along the bottom edge are off

        ax.axvline(x=0, lw=2., color='black')

        y_lims = ax.axes.get_xlim()
        w = y_lims[1] - y_lims[0]
        x_min = -0.007*w
        x_max = -x_min
        cpt = 0

        for h in bins[2:-1]:
            if cpt % 2. == 0:
                l = mlines.Line2D([x_min, x_max], [h, h], linewidth=1.5, color='black')
                ax.add_line(l)
                ax.text(- 0.015*w, h, str(int(h)), horizontalalignment='right', verticalalignment='center', fontsize=25)
            cpt += 1

        # ax.plot(right[0][0], 0.5*(right[1][0] + right[1][1]), color='black', marker='o')

        # confidence intervals
        print "######################################################"
        print "Quantitative results for the TB age distribution in " + country + ":"
        my_grey = (83./255., 81./255., 84./255.)
        coef = multiplier * self.output_objects[country].model_runner.data.console['n_runs'] * self.output_objects[country].model_runner.nb_seeds
        tb_ages_by_run = self.output_objects[country].model_runner.model_diagnostics[scenario]['all_tb_ages']
        average_nb_tb_cases = float(len(tb_ages)) / \
                               float(self.output_objects[country].model_runner.data.console['n_runs'] * self.output_objects[country].model_runner.nb_seeds)

        max_prop_mean = 0.
        for age_low in bins[0:-1]:
            age_high = age_low+5.
            if age_low == 95.:
                age_high = 150.
            nb_in_bin = [len([x for x in sublist if age_low <= x < age_high]) for sublist in tb_ages_by_run]
            deno_in_bin = [len(tb_ages_by_run[j]) for j in range(len(tb_ages_by_run))]
            if 0 in deno_in_bin:
                deno_in_bin = [max(deno_in_bin[h], 1) for h in range(len(deno_in_bin))]  # prov
            prop_in_bin = [float(nb_in_bin[j])/float(deno_in_bin[j]) for j in range(len(tb_ages_by_run))]
            prop_mean = np.mean(prop_in_bin)
            max_prop_mean = max(max_prop_mean, prop_mean)

            ax.add_patch(patches.Rectangle(xy=(0., age_low), width=prop_mean*coef*average_nb_tb_cases, height=5.,
                                           color=colors['right']))
            if self.output_objects[country].model_runner.data.console['n_runs'] > 1 or\
                            self.output_objects[country].model_runner.nb_seeds > 1:
                h = age_low + 2.5
                # sd = np.std(nb_in_bin)
                sd_prop = np.std(prop_in_bin)
                sample_size = len(tb_ages_by_run)

                low_prop = prop_mean - 1.96 * sd_prop / math.sqrt(float(sample_size))
                low_prop = max([0., low_prop])
                high_prop = prop_mean + 1.96 * sd_prop / math.sqrt(float(sample_size))

                l = mlines.Line2D([low_prop*coef*average_nb_tb_cases, high_prop*coef*average_nb_tb_cases], [h, h], linewidth=2.5, color=my_grey)
                ax.add_line(l)

                # ax.plot(prop_mean*coef*average_nb_tb_cases, h, 'bs')

                perc_low = round(100. * low_prop, 2)
                perc_high = round(100. * high_prop, 2)

            string = "Ages " + str(int(age_low)) + " to " + str(int(age_high)) + ": " + str(round(100.*prop_mean, 2)) + "%"
            if self.output_objects[country].model_runner.data.console['n_runs'] > 1 or\
                            self.output_objects[country].model_runner.nb_seeds > 1:
                string += " (" +\
                  str(perc_low) + " - " + str(perc_high) + ")"
            print string

        # x axis
        x_props = np.linspace(-0.1, 0.1, num=5)
        if country == 'China':
            x_props = np.linspace(-0.1, 0.15, num=6)
            plt.xlim((-0.1, 0.18*coef*average_nb_tb_cases))

        x_labs = [str(int(abs(100.*x))) for x in x_props]
        x_locs = [x * coef*average_nb_tb_cases for x in x_props]
        plt.xticks(x_locs, x_labs, size=20.)
        plt.xlabel('%', size=30.)

        # x_tick_pos = np.linspace(1.5, 17.5, num=9)
        # x_tick_labs = np.linspace(1, 17, num=11)
        # x_tick_labs = [int(u) for u in x_tick_labs]
        #
        # plt.xticks(x_tick_pos, x_tick_labs, fontsize=30, color='black')

    def get_age_specific_incidence(self, country, age_min, age_max):
        scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
        if len(self.output_objects[country].model_runner.data.console['years_checkpoints']) == 0:
            print "No checkpoint found. tb_age_pyramid could not be generated."
            return
        checkpoint = \
        self.output_objects[country].model_runner.model_diagnostics[scenario]['checkpoint_outcomes']['ages'].keys()[-1]
        ages = [value for i_run in range(
            self.output_objects[country].model_runner.data.console['n_runs'] * self.output_objects[
                country].model_runner.nb_seeds)
                  for value in
                  self.output_objects[country].model_runner.model_diagnostics[scenario]['checkpoint_outcomes']['ages'][
                      checkpoint][i_run]]

        tb_ages = [age for sublist in
                   self.output_objects[country].model_runner.model_diagnostics[scenario]['all_tb_ages'] for age in
                   sublist]

        ages_in_category = [a for a in ages if age_min <= a <= age_max]
        tb_ages_in_category = [a for a in tb_ages if age_min <= a <= age_max]

        if len(ages_in_category) > 0:
            prop = float(len(tb_ages_in_category)) / float(len(ages_in_category))
        else:
            print "nobody lives in this age category"
            return None

        incidence = prop * 1.e5 / 5.  # hard-coded  divided by 5 because 5 years
        print "TB incidence among " + str(int(age_min)) + "-" + str(int(age_max)) + " in " + country + ": " +\
              str(round(incidence, 2)) + " /100,000/y"

    def make_contribution_graphs(self):
        self.i_figure += 1
        plt.figure(self.i_figure, figsize=(40., 20.))
        gs = gridspec.GridSpec(2, 3)
        gs.update(hspace=0.3)

        colors = {'contact': (114./256., 147./256., 203./256.), 'transmission': (225./256., 151./256., 76./256.),
                  'transmission_end_tb': (211./256., 94./256., 96./256.) }

        for i, country in enumerate(self.countries):
            col = i % 3
            row = 0
            if i > 2:
                row = 1
            ax = plt.subplot(gs[row, col])
            self.make_a_contribution_graph(ax, country, colors)

        # legend
        ax = plt.subplot(gs[1, 2])
        ax.add_patch(patches.Rectangle(xy=(0., 0.6), width=0.16, height=0.08, color=colors['contact']))
        ax.add_patch(patches.Rectangle(xy=(0., 0.4), width=0.16, height=0.08, color=colors['transmission']))
        ax.add_patch(patches.Rectangle(xy=(0., 0.2), width=0.16, height=0.08, color=colors['transmission_end_tb']))

        ax.text(x=0.21, y=0.6, s='Social contacts', fontsize=40, verticalalignment='bottom')
        ax.text(x=0.21, y=0.4, s='Transmission events', fontsize=40, verticalalignment='bottom')
        ax.text(x=0.21, y=0.2, s='Transmission leading to active TB', fontsize=40, verticalalignment='bottom')

        ax.axis('off')

        filename = 'multi_country_contributions.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_a_contribution_graph(self, ax, country, colors):
        scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
        locations = ['home', 'school', 'work', 'other']
        keys = ['household', 'school', 'workplace', 'community']
        event_types = ['contact', 'transmission', 'transmission_end_tb']
        x_tick_pos = np.arange(1, 5)
        x_positions = {'contact': np.arange(0.75, 4.75),
                      'transmission': np.arange(1, 5),
                      'transmission_end_tb': np.arange(1.25, 5.25)
                      }
        width = 0.25

        ylabs = {'contact': 'contribution to overall social mixing (%)',
                 'transmission': 'contribution to overall transmission (%)',
                 'transmission_end_tb': 'contribution to overall TB burden(%)'}

        plt.title(country, fontsize=40, y=1.03)

        print "###############################"
        print "Quantitative contributions for " + country + ":"
        for i, key_of_interest in enumerate(event_types):
            print "********** " + key_of_interest + " ***********"
            perc_dict = copy.deepcopy(self.output_objects[country].model_runner.model_diagnostics[scenario]['n_contacts'][key_of_interest])

            for i_run in range(self.output_objects[country].model_runner.data.console['n_runs'] * self.output_objects[country].model_runner.nb_seeds):
                sum_trans = 0
                for key in keys:
                    sum_trans += perc_dict[key][i_run]
                if sum_trans > 0:
                    for key in keys:
                        perc_dict[key][i_run] = 100. * float(perc_dict[key][i_run]) / float(sum_trans)
            means = []
            stds = []
            for key in keys:
                means.append(np.mean(perc_dict[key]))
                if self.output_objects[country].model_runner.data.console['n_runs']*self.output_objects[country].model_runner.nb_seeds > 1:
                    stdd = np.std(perc_dict[key])
                    stds.append(stdd)
            if self.output_objects[country].model_runner.data.console['n_runs'] == 1 and\
                            self.output_objects[country].model_runner.nb_seeds == 1:
                stds = None

            b = ax.bar(x_positions[key_of_interest], means, width, yerr=stds, ecolor='black', align='center',
                   color=colors[key_of_interest])
            for j, location in enumerate(keys):
                mean = round(b[j]._height, 2)
                string = location + ": " + str(mean) +"%"
                if self.output_objects[country].model_runner.data.console['n_runs'] > 1 or\
                            self.output_objects[country].model_runner.nb_seeds > 1:
                    low = round(b.errorbar.lines[2][0]._paths[j]._vertices[0:2][0][1], 2)
                    high = round(b.errorbar.lines[2][0]._paths[j]._vertices[0:2][1][1], 2)
                    string += " (" + str(low) + " - " + str(high) + ")"
                print string

        plt.xlabel('')
        plt.ylabel('contribution (%)', fontsize=35, color='black')
        plt.xticks(x_tick_pos, locations, fontsize=35, color='black')
        plt.xlim(0.5, 4.5)
        plt.ylim(0., 100.)

        plt.tick_params(
            axis='y',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            left='on',
            right='off',
            labelsize=27,
            labelcolor='black')
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            top='off',
            bottom='on')

    def make_tb_burden_clocks(self):
        self.i_figure += 1
        plt.figure(self.i_figure, figsize=(40., 25.))
        gs = gridspec.GridSpec(2, 3)
        colors = {'population': (200./255., 200./255., 200./255.), 'ltbi': (132./255., 186./255., 91./255.), 'tb': (107./255., 76./255., 154./255.)}
        for i, country in enumerate(self.countries):
            col = i % 3
            row = 0
            if i > 2:
                row = 1
            ax = plt.subplot(gs[row, col]) # , projection='3d')
            self.make_a_tb_burden_clock(ax, country, colors)

        # legend
        ax = plt.subplot(gs[1, 2])

        ax.add_patch(patches.Circle(xy=(0.15, 0.55), radius=0.08, color=colors['population']))
        ax.add_patch(patches.Circle(xy=(0.15, 0.35), radius=0.06, color=colors['ltbi']))
        ax.add_patch(patches.Circle(xy=(0.15, 0.15), radius=0.04, color=colors['tb']))

        ax.text(x=0.3, y=0.55, s='Population size', fontsize=40, verticalalignment='center')
        ax.text(x=0.3, y=0.35, s='LTBI prevalence', fontsize=40, verticalalignment='center')
        ax.text(x=0.3, y=0.15, s='Active TB disease', fontsize=40, verticalalignment='center')

        plt.axis('off')

        filename = 'multi_country_clocks.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_a_tb_burden_clock(self, ax, country, colors):
        scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
        ltbi_age_stats = self.output_objects[country].model_runner.model_diagnostics[scenario]['ltbi_age_stats']

        radius_scaling = 3.  # 2. for disc effect  / 3. for sphere effect

        n_cases_by_age_group = self.process_ltbi_age_stats(ltbi_age_stats, country)

        ax.text(0., 0., country, fontsize=40, verticalalignment='center', horizontalalignment='center')
        r_max = 0.35
        r_clock_min = 0.7
        r_clock_max = 0.8
        r_clock_text = 0.6
        max_population = max(n_cases_by_age_group['population'])
        # max_n_ltbi = max(n_cases_by_age_group['n_ltbi_in_age_group'])
        print "############################"
        print "LTBI contributions (clock diagram) for " + country + ":"
        for i in range(8):
            age_min = float(i) * 10.
            angle = math.pi / 2 - ((1. + 2. * float(i)) * math.pi / 8.)
            x = math.cos(angle)
            y = math.sin(angle)

            population = float(n_cases_by_age_group['population'][i])
            rel_population = population / float(max_population)
            radius_population = r_max * (rel_population) ** (1. / radius_scaling)

            n_ltbi = n_cases_by_age_group['n_ltbi_in_age_group'][i]

            ax.add_patch(patches.Circle(xy=(x, y), radius=radius_population, color=colors['population']))
            if n_ltbi > 0:
                n_tb = n_cases_by_age_group['n_tb_in_age_group'][i]

                rel_ltbi = n_ltbi / population
                radius_ltbi = radius_population * (rel_ltbi ** (1. / radius_scaling))
                rel_tb = float(n_tb) / float(n_ltbi)
                radius_tb = radius_ltbi * (rel_tb ** (1. / radius_scaling))

                # rel_ltbi = float(n_ltbi)/float(max_n_ltbi)
                # rel_tb = float(n_tb) / float(n_ltbi)

                # radius_ltbi = r_max*math.sqrt(rel_ltbi)
                # radius_tb = radius_ltbi * math.sqrt(rel_tb)

                # radius_ltbi = r_max * (rel_ltbi ** (1./3.))
                # radius_tb = radius_ltbi * (rel_tb ** (1./3.))

                ax.add_patch(patches.Circle(xy=(x, y), radius=radius_ltbi, color=colors['ltbi']))
                ax.add_patch(patches.Circle(xy=(x, y), radius=radius_tb, color=colors['tb']))

            # draw a clock
            angle_clock = math.pi / 2 - (float(i) * math.pi / 4.)
            l = mlines.Line2D([r_clock_min*math.cos(angle_clock), r_clock_max*math.cos(angle_clock)],
                              [r_clock_min*math.sin(angle_clock), r_clock_max*math.sin(angle_clock)], color='black')
            ax.add_line(l)
            ax.text(r_clock_text * math.cos(angle_clock), r_clock_text * math.sin(angle_clock), str(int(age_min)),
                    fontsize=30, color='black', va='center', ha='center')

            # print
            age_max = age_min + 10. if age_min < 70 else 150
            string = "Ages " + str(int(age_min)) + " to " + str(int(age_max)) + ": "
            if n_ltbi > 0:
                ltbi_prev_perc = round(100. * float(n_ltbi) / float(population), 2)
                perc_of_all_ltbi = round(100. * float(n_ltbi) / float(sum(n_cases_by_age_group['n_ltbi_in_age_group'])),
                                         2)
                perc_tb_among_ltbi = round(100. * float(n_tb) / float(n_ltbi), 2)
                if n_tb > 0:
                    perc_tb_among_all_tb = round(
                        100 * float(n_tb) / float(sum(n_cases_by_age_group['n_tb_in_age_group'])), 2)
                else:
                    perc_tb_among_all_tb = 0
            else:
                ltbi_prev_perc = 0.
                perc_of_all_ltbi = 0
                perc_tb_among_ltbi = 0
                perc_tb_among_all_tb = 0
            string += "LTBI prev of " + str(ltbi_prev_perc) + "% / "
            string += "contributing " + str(perc_of_all_ltbi) + "% of all LTBI / among which " + str(
                perc_tb_among_ltbi) + "% will ever activate TB"
            string += "  /  represents " + str(perc_tb_among_all_tb) + "% of future TB burden"
            print string
        ax.text(0., 0.9, '100', fontsize=30, color='black', va='center', ha='center')

        plt.xlim((-1.3, 1.3))
        plt.ylim((-1.3, 1.3))
        plt.axis('off')

    def make_a_tb_burden_clock_3d(self, ax, country, colors):
        scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
        ltbi_age_stats = self.output_objects[country].model_runner.model_diagnostics[scenario]['ltbi_age_stats']

        n_cases_by_age_group = self.process_ltbi_age_stats(ltbi_age_stats, country)

        ax.text(0., 0., 0., country, fontsize=40, verticalalignment='center', horizontalalignment='center')
        r_max = 0.35
        r_clock_min = 0.7
        r_clock_max = 0.8
        r_clock_text = 0.6
        max_population = max(n_cases_by_age_group['population'])
        # max_n_ltbi = max(n_cases_by_age_group['n_ltbi_in_age_group'])
        print "############################"
        print "LTBI contributions (clock diagram) for " + country + ":"
        for i in range(8):
            age_min = float(i) * 10.
            angle = math.pi/2 - ((1. + 2. * float(i)) * math.pi/8.)
            x = math.cos(angle)
            y = math.sin(angle)

            population = float(n_cases_by_age_group['population'][i])
            rel_population = population / float(max_population)
            radius_population = r_max*(rel_population)**(1./3.)

            n_ltbi = n_cases_by_age_group['n_ltbi_in_age_group'][i]

            # ax.add_patch(patches.Circle(xy=(x, y), radius=radius_population, color=colors['population']))
            ax.scatter(x, y, 0., s=10000., color='g', alpha=0.5)

            if n_ltbi > 0:
                n_tb = n_cases_by_age_group['n_tb_in_age_group'][i]

                rel_ltbi = n_ltbi / population
                radius_ltbi = radius_population * (rel_ltbi**(1./3.))
                rel_tb = float(n_tb) / float(n_ltbi)
                radius_tb = radius_ltbi * (rel_tb ** (1./3.))

                # rel_ltbi = float(n_ltbi)/float(max_n_ltbi)
                # rel_tb = float(n_tb) / float(n_ltbi)

                # radius_ltbi = r_max*math.sqrt(rel_ltbi)
                # radius_tb = radius_ltbi * math.sqrt(rel_tb)

                # radius_ltbi = r_max * (rel_ltbi ** (1./3.))
                # radius_tb = radius_ltbi * (rel_tb ** (1./3.))

                # ax.add_patch(patches.Circle(xy=(x, y), radius=radius_ltbi, color=colors['ltbi']))
                # ax.add_patch(patches.Circle(xy=(x, y), radius=radius_tb, color=colors['tb']))



            # draw a clock
            angle_clock = math.pi/2 - (float(i) * math.pi/4.)
            # l = mlines.Line2D([r_clock_min*math.cos(angle_clock), r_clock_max*math.cos(angle_clock)],
            #                   [r_clock_min*math.sin(angle_clock), r_clock_max*math.sin(angle_clock)], color='black')
            # ax.add_line(l)
            ax.text(r_clock_text*math.cos(angle_clock), r_clock_text*math.sin(angle_clock), 0., str(int(age_min)),
                     fontsize=30, color='black', va='center', ha='center')

            # print
            age_max = age_min + 10. if age_min < 70 else 150
            string = "Ages " + str(int(age_min)) + " to " + str(int(age_max)) + ": "
            if n_ltbi > 0:
                ltbi_prev_perc = round(100. * float(n_ltbi) / float(population), 2)
                perc_of_all_ltbi = round(100.*float(n_ltbi)/float(sum(n_cases_by_age_group['n_ltbi_in_age_group'])), 2)
                perc_tb_among_ltbi = round(100. * float(n_tb) / float(n_ltbi), 2)
                if n_tb > 0:
                    perc_tb_among_all_tb = round(100*float(n_tb)/float(sum(n_cases_by_age_group['n_tb_in_age_group'])), 2)
                else:
                    perc_tb_among_all_tb = 0
            else:
                ltbi_prev_perc = 0.
                perc_of_all_ltbi = 0
                perc_tb_among_ltbi = 0
                perc_tb_among_all_tb = 0
            string += "LTBI prev of " + str(ltbi_prev_perc) + "% / "
            string += "contributing " + str(perc_of_all_ltbi) + "% of all LTBI / among which " + str(perc_tb_among_ltbi) + "% will ever activate TB"
            string += "  /  represents " + str(perc_tb_among_all_tb) + "% of future TB burden"
            print string
        ax.text(0., 0.9, 0., '100', fontsize=30, color='black', va='center', ha='center')

        plt.xlim((-1.3, 1.3))
        plt.ylim((-1.3, 1.3))
        ax.set_zlim(-1., 1.)
        ax.view_init(elev=90., azim=0.)
        ls = LightSource(azdeg=315, altdeg=45)

        # plt.axis('off')

    def process_ltbi_age_stats(self, ltbi_age_stats, country):
        out_dict = {'population': [], 'n_ltbi_in_age_group': [], 'n_tb_in_age_group': []}  # index 0 for age cat 0-10
        scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
        checkpoint = self.output_objects[country].model_runner.model_diagnostics[scenario]['checkpoint_outcomes']['ages'].keys()[-1]

        ages = [value for i_run in range(
            self.output_objects[country].model_runner.data.console['n_runs'] * self.output_objects[
                country].model_runner.nb_seeds)
                for value in
                self.output_objects[country].model_runner.model_diagnostics[scenario]['checkpoint_outcomes']['ages'][
                    checkpoint][i_run]]

        for i in range(8):
            age_min = float(i) * 10.
            age_max = age_min + 10. if age_min < 70. else 100.


            out_dict['population'].append(sum(1 for age in ages if age_min <= age < age_max))

            out_dict['n_ltbi_in_age_group'].append(sum(1 for j in range(len(ltbi_age_stats['ltbi_ages'])) for age in ltbi_age_stats['ltbi_ages'][j] if age_min <= age < age_max))

            out_dict['n_tb_in_age_group'].append(sum(1 for j in range(len(ltbi_age_stats['ending_tb_ages'])) for age in ltbi_age_stats['ending_tb_ages'][j] if age_min <= age < age_max))

        return out_dict

    def make_multi_household_sizes(self):
        self.i_figure += 1
        plt.figure(self.i_figure, figsize=(40., 25.))
        gs = gridspec.GridSpec(2, 3)
        gs.update(hspace=0.3)

        for i, country in enumerate(self.countries):
            col = i % 3
            row = 0
            if i > 2:
                row = 1
            plt.subplot(gs[row, col])
            scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
            self.output_objects[country].plot_age_and_household_distributions(scenario=scenario, country=country,
                                                                              save=False, keys=['household_sizes'])

        filename = 'multi_household_sizes.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_multi_scale_up_functions_prev(self):
        self.i_figure += 1
        plt.figure(self.i_figure, figsize=(30., 20.))
        gs = gridspec.GridSpec(2, 2)

        colors = {'Philippines': (57. / 256., 106. / 256., 177. / 256.),
                  'China': (204. / 256., 37. / 256., 41. / 256.),
                  'India': (218. / 256., 124. / 256., 48. / 256.),
                  'Indonesia': (83. / 256., 81. / 256., 84. / 256.),
                  'Pakistan': (62. / 256., 150. / 256., 81. / 256.)
                  }

        time_in_years = np.linspace(1950., 2020., 1000.)
        for i, key in enumerate(['bcg_coverage_prop', 'cdr_prop', 'treatment_success_prop']):
            col = i % 2
            row = 0
            if i > 1:
                row = 1
            ax = plt.subplot(gs[row, col])
            plt.xlim((1950., 2020.))
            plt.ylim((0., 100.))
            plt.xlabel('time (years)')
            plt.ylabel(self.output_objects[self.countries[0]].dictionary[key])
            for i, country in enumerate(self.countries):
                func = self.output_objects[country].model_runner.data.scale_up_functions[key]
                if key == 'bcg_coverage_prop':
                    dataset = self.output_objects[country].model_runner.data.data_from_sheets['bcg']
                elif key == 'cdr_prop':
                    dataset = self.output_objects[country].model_runner.data.data_from_sheets['gtb_2016']['c_cdr']
                elif key == 'treatment_success_prop':
                    dataset = self.output_objects[country].model_runner.data.data_from_sheets['outcomes']['c_new_tsr']
                else:
                    dataset = {}

                y = [100. * func(val) for val in time_in_years]
                plt.plot(time_in_years, y, color=colors[country], linewidth=3.)

                for year, value in dataset.iteritems():
                    plt.plot(year, value, marker='o', color=colors[country], markersize=6.)

        filename = 'multi_country_scale_ups.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def make_multi_scale_up_functions(self):
        self.i_figure += 1
        nb_countries = len(self.countries)
        plt.figure(self.i_figure, figsize=(3.7, (nb_countries + 0.5) * 0.9))

        h_ratios = [0.5]
        for _ in self.countries:
            h_ratios.append(5.)

        gs = gridspec.GridSpec(nb_countries + 1, 4,
                               width_ratios=[0.7, 1., 1., 1.],
                               height_ratios=h_ratios
                               )
        gs.update(wspace=0.12, hspace=0.1)

        #  panel titles
        keys = ['bcg_coverage_prop', 'cdr_prop', 'treatment_success_prop']  #, 'sp_prop']
        row = 0

        short_dic = {'cdr_prop': 'case detection', 'treatment_success_prop': 'treatment success',
                    'bcg_coverage_prop': 'BCG coverage', 'sp_prop': 'smear-positive'}
        for i, key in enumerate(keys):
            plt.subplot(gs[row, i + 1])
            plt.text(x=0.5, y=0.9, s=short_dic[key], fontsize=5,
                     horizontalalignment='center', verticalalignment='center', )
            plt.grid(False)
            plt.axis('off')

        time_in_years = np.linspace(1950., 2020., 1000.)
        for country in self.countries:
            row += 1
            x_axis = False
            if country == self.countries[-1]:
                x_axis = True

            # country name
            col = 0
            plt.subplot(gs[row, col])
            plt.text(x=-0.2, y=0.5, s=country, fontsize=5)
            plt.grid(False)
            plt.axis('off')

            for i, key in enumerate(keys):
                col = i + 1
                ax = plt.subplot(gs[row, col])
                ax.grid(linewidth=0.3)
                plt.xlim((1950., 2020.))
                plt.ylim((0., 100.))
                if x_axis:
                    plt.xlabel('time (years)', fontsize=5, labelpad=0., color='black')
                y_axis = False
                if i == 0:
                    y_axis = True
                    plt.ylabel('%', fontsize=5, labelpad=0., color='black')

                func = self.output_objects[country].model_runner.data.scale_up_functions[key]
                if key == 'bcg_coverage_prop':
                    dataset = self.output_objects[country].model_runner.data.data_from_sheets['bcg']
                elif key == 'cdr_prop':
                    dataset = self.output_objects[country].model_runner.data.data_from_sheets['gtb_2016']['c_cdr']
                elif key == 'treatment_success_prop':
                    dataset = self.output_objects[country].model_runner.data.data_from_sheets['outcomes']['c_new_tsr']
                elif key == 'sp_prop':
                    # smear-positive proportion separately
                    sp_number = self.output_objects[country].model_runner.data.data_from_sheets['outcomes']['new_sp_coh']
                    snep_number = self.output_objects[country].model_runner.data.data_from_sheets['outcomes']['new_snep_coh']
                    dataset = {}
                    for year in snep_number.keys():
                        if year in sp_number.keys() and (snep_number[year] + sp_number[year]) > 0.:
                            dataset[year] = 100. * sp_number[year] / (sp_number[year] + snep_number[year])
                else:
                    dataset = {}

                y = [100. * func(val) for val in time_in_years]
                plt.plot(time_in_years, y, color='black', linewidth=0.8)

                for year, value in dataset.iteritems():
                    plt.plot(year, value, 'ro', markersize=0.5)

                plt.tick_params(
                    axis='x',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    bottom=x_axis,  # ticks along the bottom edge are off
                    top='off',  # ticks along the top edge are off
                    labelbottom=x_axis,
                    length=2.,
                    pad=0.5,  # distance tick - label
                    labelsize=4,
                    colors='black'
                )  # labels along the bottom edge are off
                plt.tick_params(
                    axis='y',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    left=y_axis,  # ticks along the bottom edge are off
                    right='off',  # ticks along the top edge are off
                    labelleft=y_axis,
                    length=2.,
                    pad=0.5,
                    labelsize=4,
                    colors='black'
                )  # labels along the bottom edge are off

        filename = 'multi_country_scale_ups.pdf'  # + self.model_runner.data.console['graph_format']
        file_path = os.path.join(self.multi_base_path, filename)
        plt.savefig(file_path)
        plt.close()

    def write_ltbi_prevalences(self, year=2014.):
        for country in self.countries:
            scenario = self.output_objects[country].model_runner.data.scenarios.keys()[0]
            ltbi_data = self.output_objects[country].model_runner.model_diagnostics[scenario]['timeseries']['ltbi_prevalence']
            times = self.output_objects[country].model_runner.model_diagnostics[scenario]['timeseries']['times'][0, :]
            times = self.output_objects[country].convert_model_time_to_dates(times)

            time_index = 0
            for time in times:
                if time < year:
                    time_index += 1
                else:
                    break

            ltbi_prevs = ltbi_data[:, time_index]

            mean_prev = np.mean(ltbi_prevs)
            sd = np.std(ltbi_prevs)

            # up = mean_prev + 1.96*sd/np.sqrt(float(len(ltbi_prevs)))
            # low = mean_prev - 1.96*sd/np.sqrt(float(len(ltbi_prevs)))

            low = np.percentile(ltbi_prevs, 2.5)
            up = np.percentile(ltbi_prevs, 97.5)

            print "LTBI prevalence (%) in " + str(year) + " for " + country + ":"
            print str(mean_prev) + " (" + str(low), "- " + str(up) + ")"


def get_lhs_param_values(path):

    with open(path) as csv_data:
        reader = csv.reader(csv_data)

        # eliminate blank rows if they exist
        rows = [row for row in reader if row]
        headings = rows[0] # get headings

        param_dict = {}
        for row in rows[1:]:
            # append the dataitem to the end of the dictionary entry
            # set the default value of [] if this key has not been seen
            for col_header, data_column in zip(headings, row):
                param_dict.setdefault(col_header, []).append(data_column)

    return param_dict



if __name__ == "__main__":

    # for country in ['India','Indonesia','Philippines', 'Pakistan']:
    #     Out = load_outputs('outputs/analysis_100s12r_Jul2019_' + country + '/pickled_outputs.pickle')
    #     Out.make_sensitivity_plot()
    #
    # Out = load_outputs('outputs/test_20K_pop_Indonesia/pickled_outputs.pickle')
    # Out.make_graphs()
    # #
    # exit()

    # ['India', 'Indonesia', 'China', 'Philippines', 'Pakistan']


    countries = ['India', 'Indonesia', 'China', 'Philippines', 'Pakistan']
    # countries = ['Indonesia']
    MO = multi_output(countries=countries, project_name='analysis_100s12r_Jul2019')
    MO.make_heatmaps_figure()
    # # MO.make_double_pyramids()
    #
    # MO.make_multi_graphs()
    # MO.make_country_graphs()  # check the folder path in local console spreadsheet
    # MO.write_ltbi_prevalences()
    #
    # for country in countries:
    #     print "###############################"
    #     print country
    #     scenario = MO.output_objects[country].model_runner.data.scenarios.keys()[0]
    #     MO.output_objects[country].write_prevalence_by_age_new(scenario)
    #     print ""
