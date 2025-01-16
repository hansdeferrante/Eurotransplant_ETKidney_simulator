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
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


if __name__ == '__main__':
    pd.set_option('display.max_rows', 500)

    sim_set = read_sim_settings(
        os.path.join(
            es.DIR_SIM_SETTINGS,
            '2025-01-16',
            'CurrentETKAS_vPRA_sliding_scale_b2_100p_1_1.yml'
        )
    )
    # Read in simulation settings
    simulator = SimulationRunner(
        sim_set=sim_set,
        verbose=0
    )

    ml = simulator.simulate_allocation(
        verbose=False
    )

    simulator.sim_results.save_outcomes_to_file(
        patients=simulator.patients,
        cens_date=sim_set.SIM_END_DATE
    )
