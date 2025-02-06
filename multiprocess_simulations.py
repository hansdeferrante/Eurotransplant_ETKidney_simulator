#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu 20 10 2022

@author: H.C. de Ferrante
"""

from __future__ import division
from multiprocessing import Pool, Lock
from time import time, sleep
from typing import Tuple, List, Optional
import os
import simulator.magic_values.etkidney_simulator_settings as es
from pathlib import Path
import re
from random import shuffle
import warnings
import pandas as pd
import argparse
import logging
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

# Configure logging
logging.basicConfig(filename='simulation_errors.log', level=logging.ERROR,
                    format='%(asctime)s %(levelname)s:%(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Script to configure n_workers.")
    parser.add_argument(
        "-n", "--n_workers",
        type=int,
        default=8,
        help="Number of workers to use. Default is 8."
    )
    parser.add_argument(
        "-t", "--timeout",
        type=int,
        default=0,
        help=(
            "Number of seconds to sleep before starting the simulation. "
            "Default is 0."
        )
    )
    return parser.parse_args()


def list_full_paths(directories: Tuple[str, List[str]]) -> Optional[List[str]]:
    if isinstance(directories, str):
        return [
            os.path.join(directories, file)
            for file in os.listdir(directories)
        ]
    elif isinstance(directories, list):
        return [
            os.path.join(directory, file)
            for directory in directories
            for file in os.listdir(directory)
        ]
    else:
        return None


def simulate_allocation(sim_path, skip_if_exists: bool = False):

    from simulator.code.SimulationRunner import SimulationRunner
    from simulator.code.utils.read_input_files import read_sim_settings

    # Read in simulation settings
    with lock_init:

        try:

            sim_set = read_sim_settings(
                sim_path
            )

            patient_file = Path(
                f'{sim_set.RESULTS_FOLDER}/{sim_set.PATH_FINAL_PATIENT_STATUS}'
            )
            transplant_file = Path(
                f'{sim_set.RESULTS_FOLDER}/{sim_set.PATH_TRANSPLANTATIONS}'
            )
            if patient_file.is_file() and transplant_file.is_file() and skip_if_exists:
                print(f'{sim_path} already processed... Continuing with next one')
                return 0
            elif not Path(sim_set.PATH_STATUS_UPDATES).is_file():
                print(
                    f"Status update file '{sim_set.PATH_STATUS_UPDATES} "
                    f"does not exist. Skipping...")
                return 0
            else:
                print(f'Working on {sim_path}')

            # Initialize simulation
            print(f"Lock acquired for {sim_path}")
            simulator = SimulationRunner(
                sim_set=sim_set,
                verbose=False
            )

        except Exception as e:
            print('\n\n***********')
            msg = f'An error occurred when initializing for {sim_path}: {e}'
            print(msg)
            logging.exception(msg)
            print('\n\n***********')
            return 0

    print(f"Lock released for {sim_path}")
    try:
        simulator.simulate_allocation(
            verbose=False,
            print_progress_every_k_days=365
        )
    except Exception as e:
        print('\n\n***********')
        msg = f'An error occurred when simulating for {sim_path}: {e}'
        print(msg)
        logging.exception(msg)
        print('\n\n***********')
        return 0

    with lock_gzip:
        try:

            simulator.sim_results.save_outcomes_to_file(
                patients=simulator.patients,
                cens_date=sim_set.SIM_END_DATE
            )
            return 1
        except:
            print('\n\n***********')
            msg = f'An error occurred when saving outcomes for {sim_path}: {e}'
            print(msg)
            logging.exception(msg)
            print('\n\n***********')
            return 0

    return 1


def init(il, gl):
    global lock_gzip
    lock_gzip = gl
    global lock_init
    lock_init = il


if __name__ == '__main__':

    args = parse_arguments()

    timeout = args.timeout      # Get the timeout duration

    # Sleep for the specified timeout duration
    print(f"Sleeping for {timeout} seconds before starting the simulation...")
    sleep(timeout)

    n_workers = args.n_workers
    print(f'Multi-processing on {n_workers} cores')

    paths = list_full_paths(
        [
            os.path.join(es.DIR_SIM_SETTINGS, path)
            for path in ['2025-02-06']
        ]
    )
    shuffle(paths)

    pattern = re.compile('.*validation.*\\.yml')
    pattern = re.compile('.*\\.yml')
    paths = [path for path in paths if pattern.match(path)]

    paths.sort(key=lambda x: -len(x.split('validation')))

    init_lock = Lock()
    gzip_lock = Lock()

    start = time()
    print(f'Processing {len(paths)} files')

    with Pool(
        processes=n_workers,
        initializer=init,
        initargs=(init_lock, gzip_lock)
    ) as p:

        for i, _ in enumerate(p.imap(simulate_allocation, paths)):
            print('\rdone {0:%}'.format(i / len(paths)))

    end = time()
    print(f'{len(paths)} jobs executed in {(end - start) / 60:.1f} minutes')
