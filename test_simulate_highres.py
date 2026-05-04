#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

import os
import pandas as pd

from simulator.code.SimulationRunner import SimulationRunner
from simulator.code.utils.read_input_files import read_sim_settings
import simulator.magic_values.etkidney_simulator_settings as es
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


if __name__ == '__main__':
    pd.set_option('display.max_rows', 500)

    sim_set = read_sim_settings(
        os.path.join(
            es.DIR_SIM_SETTINGS,
            '2026-05-04',
            'CurrentETKAS_actual_1_1.yml'
        )
    )
    sim_set.PATH_PATIENTS = 'data/final_data_all_high_res_african_random.csv'
   # sim_set.PATH_DONORS = 'data/fake_donors_for_etkidney_simulator_hr.csv'

    # Read in simulation settings
    simulator = SimulationRunner(
        sim_set=sim_set,
        verbose=0
    )

    simulator.simulate_allocation(
        verbose=False
    )

    simulator.sim_results.save_outcomes_to_file(
        patients=simulator.patients,
        cens_date=sim_set.SIM_END_DATE
    )
