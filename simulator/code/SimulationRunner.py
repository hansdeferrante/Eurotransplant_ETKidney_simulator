#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 16-08-2022

@author: H.C. de Ferrante
"""

import os
from datetime import timedelta, datetime
from typing import Optional, Dict, List, ValuesView
import time
import typing

from simulator.code.utils.utils import DotDict, round_to_decimals
import simulator.code.utils.read_input_files as rdr
from simulator.code.matchlist.MatchListETKASandESP import \
    MatchListCurrentETKAS, \
    MatchListESP, PatientMatchRecord
from simulator.code.matchlist.MatchList import MatchDistanceCache
from simulator.code.matchlist.ScoringFunction import MatchPointFunction
from simulator.code.AcceptanceModule import AcceptanceModule
from simulator.code.PostTransplantPredictor import PostTransplantPredictor
from simulator.code.events.EventQueue import EventQueue
from simulator.code.events.Event import Event
from simulator.code.utils.load_simulation_entities import \
    preload_profiles, preload_status_updates, load_patients, \
    load_donors, load_retransplantations, \
    load_balances, load_nonetkasesp_balances
from simulator.code.entities import (
    Patient, Donor
)
from simulator.code.BalanceSystem import BalanceSystem
from simulator.code.HLA.HLASystem import HLASystem
from simulator.code.HLA.MMPSystem import MMPSystem
from simulator.code.utils.SimResults import SimResults
import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es

if typing.TYPE_CHECKING:
    from simulator.code.matchlist import MatchList


class ActiveList:
    """
        Class which tracks actively listed patients, per blood type.
        The order of the list is on dialysis time (if it exists).
        Retrieving patients from the active lists prevents having to
        iterate over all patients, and speeds up sorting.
    """

    def __init__(self, init_patients: Dict[int, Patient],
                 bal_system: BalanceSystem,
                 match_point_fun: MatchPointFunction):

        # Keep track of necessary attributes
        self.bal_system = bal_system
        self.match_point_fun = match_point_fun

        # Order active lists by blood type
        self._active_lists: Dict[str, Dict[int, Patient]] = {
            bg: dict() for bg in es.ALLOWED_BLOODGROUPS
        }
        for id_reg, pat in init_patients.items():
            if pat.active:
                self._active_lists[pat.bloodgroup][id_reg] = pat

        # Initialize sort
        self.resort_lists(match_time=0)

    def get_active_list(self, bg: str) -> ValuesView[Patient]:
        """Returning patients, ordered with priority for bg"""
        return self._active_lists[bg].values()

    def add_to_active_list(self, identifier: int, pat: Patient) -> None:
        """Adding an identifier to the active list"""
        self._active_lists[pat.bloodgroup][identifier] = pat

    def pop_from_active_lists(self, identifier: int) -> None:
        """Remove a patient identifier from active list"""
        for active_list in self._active_lists.values():
            active_list.pop(identifier, None)

    def is_active(self, identifier: int) -> bool:
        for active_list in self._active_lists.values():
            if identifier in active_list:
                return True
            else:
                return False
        return False

    def resort_lists(self, match_time: float):
        """
            Sort active list based on patient match points.
            This starting order makes sorting faster for the match
            lists.
        """
        for bg, active_list in self._active_lists.items():
            match_point_dict = {
                id: PatientMatchRecord(
                    pr,
                    calc_points=self.match_point_fun,
                    bal_system=self.bal_system,
                    match_time=match_time
                ).total_match_points
                for id, pr in active_list.items()
            }
            # Sort dictionary by patient match points
            self._active_lists[bg] = dict(
                sorted(
                    active_list.items(),
                    key=lambda x: match_point_dict[x[0]],
                    reverse=True
                )
            )

    def return_all_active_patients(self) -> Dict[int, Patient]:
        """Return all active patients"""
        return {
            id_reg: patient
            for _, patients in self._active_lists.items()
            for id_reg, patient in patients.items()
        }


class SimulationRunner:
    """Class which implements an SimulationRunner simulator

    Attributes   #noqa
    ----------
    patients: List[Patient]
        list of patients
    sim_set: DotDict
        simulation settings
    max_n_patients: Optional[int]
        maximum number of patients to read in
    max_n_statuses: Optional[int]
        maximum number of statuses to read in
    sim_time: float
        simulation time
    ...
    """

    def __init__(
        self,
        sim_set: DotDict,
        max_n_patients: Optional[int] = None,
        verbose: Optional[int] = True
    ):

        # Set up (empty) event queue
        self.event_queue = EventQueue()

        # Read in simulation settings
        self.sim_set = sim_set
        self.sim_start_date = sim_set['SIM_START_DATE']
        self.verbose = verbose
        self.sim_rescue = sim_set.get('SIMULATE_RESCUE', False)
        self.hla_system = HLASystem(sim_set)
        self.mmp_system = MMPSystem(sim_set, self.hla_system)

        self.match_distance_cache = MatchDistanceCache(
            sim_set
        )

        self.bal_system = load_balances(sim_set)

        self.nonetkas_esp_transplants = load_nonetkasesp_balances(sim_set)
        self.allow_discards = sim_set.get('ALLOW_DISCARDS', False)

        for bal_id, rcrd in self.nonetkas_esp_transplants.items():
            self.event_queue.add(
                Event(
                    type_event=cn.BAL,
                    event_time=(
                        rcrd[cn.D_DATE] - self.sim_set.SIM_START_DATE
                    ) / timedelta(days=1),
                    identifier=bal_id
                )
            )

        # Set up times, and resort list every 14 days
        self.sim_time = 0
        self._time_active_list_resorted = 0
        self._resort_every_k_days = 14

        self.max_sim_time = (
            (
                self.sim_set.SIM_END_DATE - self.sim_start_date
            ) / timedelta(days=1)
        )

        # Add snapshot events if they are to be saved
        if self.sim_set.SNAPSHOT_FREQUENCY:
            for ss_time in (
                range(
                    0,
                    int(self.max_sim_time),
                    int(self.sim_set.SNAPSHOT_FREQUENCY)
                )
            ):
                self.event_queue.add(
                    Event(
                        type_event=cn.SNAPSHOT,
                        event_time=ss_time,
                        identifier=0
                    )
                )

        # Read in travel time dictionary
        self.dict_travel_times = rdr.read_travel_times()

        # Initialize the acceptance module
        self.acceptance_module = AcceptanceModule(
            seed=self.sim_set.SEED,
            center_acc_policy=str(sim_set.CENTER_ACC_POLICY),
            patient_acc_policy=str(sim_set.PATIENT_ACC_POLICY),
            max_recipient_driven_offers_per_center=(
                sim_set.MAX_OFFERS_PER_CENTER
            ),
            simulate_rescue=self.sim_rescue,
            verbose=verbose,
            simulate_random_effects=sim_set.SIMULATE_RANDOM_EFFECTS,
            simulate_joint_random_effects=sim_set.JOINT_RANDOM_EFFECTS,
            re_variance_components=sim_set.VARCOMPS_RANDOM_EFFECTS
        )
        if self.verbose:
            print('Loaded acceptance module')

        # Load donors
        self.donors = load_donors(
            sim_set=self.sim_set,
            hla_system=self.hla_system
        )

        # Load patients, and initialize a DataFrame with
        # blood group compatibilities.
        if self.verbose:
            print('Loading patient registrations, statuses, and profiles')
        self.patients = load_patients(
            sim_set=self.sim_set,
            nrows=max_n_patients,
            hla_system=self.hla_system,
            mmp_system=self.mmp_system
        )

        # Initialize random effects for known levels.
        self.acceptance_module._initialize_random_effects(
            re_levels={
                cn.ID_REGISTRATION: set(
                    p.__dict__[cn.ID_REGISTRATION]
                    for p in self.patients.values()
                ),
                cn.RECIPIENT_CENTER: set(
                    p.__dict__[cn.RECIPIENT_CENTER]
                    for p in self.patients.values()
                ),
                cn.ID_DONOR: set(
                    d.__dict__[cn.ID_DONOR]
                    for d in self.donors.values()
                ),
            }
        )

        # Preload status updates
        preload_status_updates(
            patients=self.patients,
            sim_set=sim_set
        )

        # Preload profile information
        preload_profiles(
            patients=self.patients,
            sim_set=sim_set
        )

        # If simulating re-transplantations, load in retransplantations
        # Note that re-transplantations future to simulation start time
        # are not loaded in by load_patients if SIM_RETX.
        if sim_set.SIM_RETX:

            # Load retransplantations
            if self.verbose:
                print('Loading retransplantations')
            self.retransplantations = load_retransplantations(
                sim_set=self.sim_set,
                hla_system=self.hla_system,
                nrows=max_n_patients,
                mmp_system=self.mmp_system
            )
            self.initialize_retransplantations()

            preload_status_updates(
                patients=self.retransplantations,
                sim_set=sim_set,
                end_date_col='LOAD_RETXS_TO'
            )
            self.remove_patients_without_statusupdates(
                pat_dict_name='retransplantations'
            )

        # Trigger all historic updates for the real patient list.
        for pat in self.patients.values():
            pat.trigger_historic_updates()

        # Remove patients from real registrations that are
        # not initialized / never have any status update.
        self.remove_patients_without_statusupdates(
            pat_dict_name='patients',
            verbosity=1
        )

        # Remove patients without valid HLA
        self.remove_patients_without_hlas(
            pat_dict_name='patients',
            verbosity=1
        )

        # Initialize EventQueue with patient updates
        # for all patients, active or not
        for pat_id, pat in self.patients.items():
            get_next_update_time = pat.get_next_update_time()
            if get_next_update_time:
                self.event_queue.add(
                    Event(
                        type_event=cn.PAT,
                        event_time=get_next_update_time,
                        identifier=pat_id
                    )
                )

        # Initialize queue of active patients
        self.active_lists = ActiveList(
            init_patients=self.patients,
            bal_system=self.bal_system,
            match_point_fun=self.sim_set.calc_etkas_score
        )

        if self.verbose:
            print('Initialized patients')

        # Initialize the post-transplant module.
        if sim_set.SIM_RETX:
            self.ptp = PostTransplantPredictor(
                seed=self.sim_set.SEED,
                offset_ids_transplants=max(
                    max(self.patients.keys()),
                    max(self.retransplantations.keys())
                ),
                retransplants=self.retransplantations,
                discrete_match_vars=es.POSTTXP_DISCRETE_MATCH_VARS,
                cvars_trafos=es.POSTTXP_TRANSFORMATIONS,
                cvars_caliper=es.POSTTXP_MATCH_CALIPERS,
                continuous_match_vars_rec=(
                    es.POSTTXP_CONTINUOUS_MATCH_VARS[cn.RETRANSPLANT]
                ),
                continuous_match_vars_off=(
                    es.POSTTXP_CONTINUOUS_MATCH_VARS[cn.OFFER]
                ),
                min_matches=es.POSTTXP_MIN_MATCHES
            )

        # Schedule events for donors.
        for don_id, don in self.donors.items():
            self.event_queue.add(
                Event(
                    type_event=cn.DON,
                    event_time=don.arrival_at(
                        self.sim_start_date
                    ),
                    identifier=don_id
                )
            )

        # Initialize simulation results & path to match list file.
        self.sim_results = SimResults(
            cols_to_save_exit=es.OUTPUT_COLS_EXITS,
            cols_to_save_discard=es.OUTPUT_COLS_DISCARDS,
            cols_to_save_patients=es.OUTPUT_COLS_PATIENTS,
            cols_to_save_snapshot=es.OUTPUT_COLS_SNAPSHOT,
            sim_set=self.sim_set
        )
        self.match_list_file = (
            str(self.sim_set.RESULTS_FOLDER) +
            str(self.sim_set.PATH_MATCH_LISTS)
        )

    def initialize_retransplantations(self):
        for retx in self.retransplantations.values():
            if (retx_dial := retx.get_dial_time_at_listing()) is not None:
                retx.__dict__[cn.YEARS_ON_DIAL] = (
                    -1 * retx_dial / es.DAYS_PER_YEAR_FLOAT
                    if retx_dial < 0
                    else 0
                )

    def remove_patients_without_statusupdates(
        self, pat_dict_name: str, verbosity: int = 0
    ) -> None:
        """Removes patients from dict without any updates."""

        n_removed = 0
        for id_reg, pat in list(self.__dict__[pat_dict_name].items()):
            if (
                not pat.is_initialized() or
                pat.__dict__[cn.EXIT_STATUS] is not None
            ):
                n_removed += 1

                del self.__dict__[pat_dict_name][id_reg]
        if verbosity:
            print(
                f'Removed {n_removed} patients from self.{pat_dict_name} '
                f"that already exited or who don't have any status updates "
            )

    def remove_patients_without_hlas(
        self, pat_dict_name: str, verbosity: int = 0
    ) -> None:
        """Removes patients from dict without any updates."""

        n_removed = 0
        for id_reg, pat in list(self.__dict__[pat_dict_name].items()):
            if not pat.hla.required_loci_known:
                n_removed += 1

                del self.__dict__[pat_dict_name][id_reg]

        if verbosity:
            print(
                f'Removed {n_removed} patients from self.{pat_dict_name} '
                f"that have no valid HLA"
            )

    def get_active_list(self, bg: str) -> ValuesView[Patient]:
        return self.active_lists.get_active_list(bg)

    def process_nonacceptance(
        self, donor: Donor,
        match_records: List['MatchList.MatchRecord'],
        n_discards: int,
        current_date: datetime
    ):
        """ Process in case the graft was declined by all patients/centers.
            We either (i) force allocation to a FP-candidate or (ii) save
            the graft as a discard.

            Note that this happens very rarely (e.g. for HCV-positive donors.)
        """

        if match_records is None or self.allow_discards:
            self.sim_results.save_discard(
                donor=donor,
                matchl_=match_records
            )
        elif match_records is not None:
            # Deduplicate list (in case candidate occurs in both ESP and ETKAS,
            # take only ETKAS record).
            filtered_match_records = list({
                off.__dict__[cn.ID_RECIPIENT]: off
                for off in match_records
                if
                (
                    off.__dict__.get(cn.PROB_ACCEPT_P) is not None or
                    off.__dict__.get(cn.PROB_ACCEPT_C) is not None
                ) and
                off.patient.__dict__.get(cn.EXIT_STATUS) != cn.FU
            }.values())

            # First offer to candidates who turned down,
            # then to centers who turned down.
            filtered_match_records.sort(
                key=lambda x: (
                    x.__dict__.get(cn.PROB_ACCEPT_P) is not None,
                    x.__dict__.get(cn.PROB_ACCEPT_P),
                    x.__dict__.get(cn.PROB_ACCEPT_C) is not None,
                    x.__dict__.get(cn.PROB_ACCEPT_C)
                ),
                reverse=True
            )
            for mr in filtered_match_records:
                if n_discards == 0:
                    break
                mr.__dict__[cn.DKT] = 0
                mr.set_acceptance(
                    cn.T3
                )
                self.process_accepted_mr(
                    current_date=current_date,
                    donor=donor,
                    accepted_mr=mr
                )
                n_discards -= 1
            if n_discards != 0:
                self.sim_results.save_discard(
                    donor=donor,
                    matchl_=match_records
                )

    def simulate_allocation(
            self,
            verbose: bool = False,
            print_progress_every_k_days=90
    ) -> Optional[MatchListCurrentETKAS]:
        """Simulate the allocation algorithm"""

        # Set seed and maintain count for how many days were simulated.
        next_k_days = 1
        start_time = time.time()
        print(
            f'Simulating ETKAS from {self.sim_start_date.strftime("%Y-%m-%d")}'
            f' to {self.sim_set.SIM_END_DATE.strftime("%Y-%m-%d")} '
            f'with ETKAS-allocation based on:\n'
            f'\t{self.sim_set["calc_etkas_score"]}'
        )

        # Start with an empty match list file.
        if self.sim_set.SAVE_MATCH_LISTS:
            if os.path.exists(self.match_list_file):
                os.remove(
                    self.match_list_file
                )
            if os.path.exists(self.match_list_file + '.gzip'):
                os.remove(
                    self.match_list_file + '.gzip'
                )

        # Simulate until simulation is finished.
        while (
            not self.event_queue.is_empty() and
            (self.sim_time < self.max_sim_time)
        ):
            # Print simulation progress
            if self.sim_time / print_progress_every_k_days >= next_k_days:
                current_date = (
                    self.sim_start_date + timedelta(days=self.sim_time)
                )
                print(f"Simulated up to {current_date.strftime('%Y-%m-%d')}")
                next_k_days += 1

            # Resort active list every k days.
            # This speeds up match list sorting.
            if (
                (self.sim_time - self._time_active_list_resorted) >
                self._resort_every_k_days
            ):
                self.active_lists.resort_lists(match_time=self.sim_time)
                self._time_active_list_resorted = self.sim_time

            # Progress to next event
            event = self.event_queue.next()
            self.sim_time = event.event_time

            if event.type_event == cn.PAT:
                if event.identifier not in self.patients:
                    continue
                # Update patient information, and schedule future event
                # if patient remains active.
                self.patients[event.identifier].do_patient_update(
                    sim_results=self.sim_results
                )
                # If patient does not exit, schedule a new future event
                # time at his next update.
                if self.patients[event.identifier].exit_status is None:
                    event.update_event_time(
                        new_time=self.patients[
                            event.identifier
                        ].get_next_update_time()
                    )
                    self.event_queue.add(
                        event
                    )

                # If patient becomes active, add to active list.
                # If inactive, pop from active list.
                if (
                    not self.active_lists.is_active(event.identifier) and
                    self.patients[event.identifier].active
                ):
                    self.active_lists.add_to_active_list(
                        event.identifier,
                        self.patients[event.identifier]
                    )
                if not self.patients[event.identifier].active:
                    self.active_lists.pop_from_active_lists(
                        event.identifier
                    )

            elif event.type_event == cn.DON:
                donor = self.donors[event.identifier]
                donor_dict = donor.__dict__

                current_date: datetime = (
                    self.sim_start_date + timedelta(days=self.sim_time)
                )

                # ESP-allocation. First do ESP allocation,
                # then continue with rescue allocation
                if donor.__dict__[cn.DONOR_AGE] >= donor.age_esp_eligible:
                    match_list_esp = MatchListESP(
                        patients=(
                            p for p in self.get_active_list(
                                donor_dict[cn.D_BLOODGROUP]
                            ) if (
                                (
                                    donor_dict[cn.D_DCD] == 0 or p.dcd_country
                                ) and
                                p.get_esp_eligible(s=self.sim_time)
                            ) and
                            not p.am and      # Do not allow TXP to AM cand.
                            p.valid_pra == 1  # Require a valid PRA
                        ),
                        donor=donor,
                        match_date=current_date,
                        sim_start_date=self.sim_start_date,
                        hla_system=self.hla_system,
                        bal_system=self.bal_system,
                        calc_points=self.sim_set.calc_esp_score,
                        store_score_components=(
                            self.sim_set.STORE_SCORE_COMPONENTS
                        ),
                        travel_time_dict=self.dict_travel_times
                    )
                    match_list_etkas = None

                    # Now determine which patient accepts the organ
                    acc_mrs_esp = (
                        self.acceptance_module.simulate_esp_allocation(
                            match_list_esp,
                            n_kidneys_available=donor.__dict__[
                                cn.N_KIDNEYS_AVAILABLE
                            ]
                        )
                    )
                    if self.sim_set.SAVE_MATCH_LISTS:
                        self.sim_results.save_match_list(match_list_esp)

                    if acc_mrs_esp:
                        for mr in acc_mrs_esp:
                            self.process_accepted_mr(
                                accepted_mr=mr,
                                current_date=current_date,
                                donor=donor
                            )
                        n_accepted = sum(
                            2 if mr.__dict__[cn.DKT]
                            else 1
                            for mr in acc_mrs_esp
                        )
                    else:
                        n_accepted = 0
                else:
                    # Allocation via ETKAS. Construct an match-list
                    # with all patients that are ETKAS eligible,
                    # not AM, and in the same blood type.
                    match_list_etkas = MatchListCurrentETKAS(
                        patients=(
                            p for p in self.get_active_list(
                                donor_dict[cn.D_BLOODGROUP]
                            ) if (
                                p.get_etkas_eligible(s=self.sim_time) and
                                # Do not allow TXP for AM candidates
                                not p.am and
                                # Only candidates with not yet expired PRA
                                # screening are selected in ETKAS
                                p.valid_pra == 1 and
                                (p.dcd_country or donor_dict[cn.D_DCD] == 0)

                            )
                        ),
                        donor=donor,
                        match_date=current_date,
                        sim_start_date=self.sim_start_date,
                        hla_system=self.hla_system,
                        bal_system=self.bal_system,
                        calc_points=self.sim_set.calc_etkas_score,
                        store_score_components=(
                            self.sim_set.STORE_SCORE_COMPONENTS
                        ),
                        travel_time_dict=self.dict_travel_times,
                        distance_cache=self.match_distance_cache
                    )
                    match_list_esp = None

                    # Now determine which patient accepts the organ
                    accepted_mrs_etkas = (
                        self.acceptance_module.simulate_etkas_allocation(
                            match_list_etkas,
                            n_kidneys_available=(
                                donor.__dict__[cn.N_KIDNEYS_AVAILABLE]
                            )
                        )
                    )

                    # Donors of ESP-age are allocated in extended allocation.
                    # We do not trigger rescue allocation for these ESP-aged
                    # donors imm., but give national / regional priority.
                    # This is a 2014 recommendation, implemented in 2021.
                    if self.sim_set.SAVE_MATCH_LISTS:
                        self.sim_results.save_match_list(match_list_etkas)

                    if accepted_mrs_etkas:
                        for mr in accepted_mrs_etkas:
                            self.process_accepted_mr(
                                accepted_mr=mr,
                                current_date=current_date,
                                donor=donor
                            )
                        n_accepted = sum(
                            2 if mr.__dict__[cn.DKT]
                            else 1
                            for mr in accepted_mrs_etkas
                        )
                    else:
                        n_accepted = 0

                if (
                    n_discards := donor.__dict__[
                        cn.N_KIDNEYS_AVAILABLE
                    ] - n_accepted
                ) != 0:
                    match_records = []
                    if match_list_esp:
                        match_records += match_list_esp.return_match_list()
                    if match_list_etkas:
                        match_records += match_list_etkas.return_match_list()

                    if match_records:
                        self.process_nonacceptance(
                            donor=donor,
                            match_records=match_records,
                            n_discards=n_discards,
                            current_date=current_date
                        )

            elif event.type_event == cn.BAL:
                self.bal_system.add_balance_from_txp(
                    rcrd=self.nonetkas_esp_transplants[event.identifier],
                    txp_time=event.event_time,
                    expiry_time=(
                        event.event_time + self.sim_set.WINDOW_FOR_BALANCE
                    )
                )
            elif event.type_event == cn.SNAPSHOT:
                self.sim_results.save_waitlist_snapshot(
                    patients=self.active_lists.return_all_active_patients(),
                    time=self.sim_time
                )
            else:
                raise ValueError(
                    f"{event.type_event} is not a valid event type."
                )
        print(
            "--- Finished simulation in {0} seconds ---" .
            format(round_to_decimals(time.time() - start_time, 1))
        )

    def process_accepted_mr(
            self,
            current_date: datetime,
            donor: Donor,
            accepted_mr: 'MatchList.MatchRecord'
    ):
        # Process an accepted match record. For this,
        # we first simulate whether
        # the transplant is an en bloc transplantation.
        accepted_mr.patient.set_transplanted(
            tx_date=current_date,
            donor=donor,
            match_record=accepted_mr,
            sim_results=self.sim_results
        )
        self.active_lists.pop_from_active_lists(
            accepted_mr.patient.id_registration
        )

        # Add transplantation to balance if it is international
        if accepted_mr.__dict__[cn.MATCH_INTERNATIONAL]:
            self.bal_system.add_balance_from_txp(
                rcrd=accepted_mr.__dict__,
                txp_time=self.sim_time,
                expiry_time=self.sim_time + self.sim_set.WINDOW_FOR_BALANCE
            )

        # Also remove future registrations for this patient if
        # they were transplanted. This is needed to remove
        # re-registrations for candidates who were delisted
        # in their primary registration.
        future_reg_ids = tuple(
            id_reg
            for id_reg, p in self.patients.items()
            if p.__dict__[cn.ID_RECIPIENT] == (
                accepted_mr.patient.__dict__[cn.ID_RECIPIENT] and
                p.listing_offset >= self.sim_time
            )
        )
        for id_reg in future_reg_ids:
            self.patients.pop(id_reg)

        # Simulate post-transplant survival, i.e. simulate
        # patient / graft failure, and add a synthetic
        # reregistration if graft fails before end
        # of simulation period.
        if self.sim_set.SIM_RETX:
            if hasattr(accepted_mr, '_initialize_posttxp_information'):
                accepted_mr._initialize_posttxp_information(
                    ptp=self.ptp
                )
            self.simulate_posttxp(
                txp=accepted_mr
            )

    def simulate_posttxp(self, txp: 'MatchList.MatchRecord'):
        """Code to simulate post-transplant outcomes for transplanted
            patients. This simulates failure and relisting dates, but
            also adds a synthetic re-registration to the simulation
            in case the re-listing occurs during the sim period.
        """

        # Simulate post-transplant survival
        date_fail, date_relist, cause_fail = self.ptp.simulate_posttransplant(
            offer=txp,
            current_date=(
                self.sim_start_date +
                timedelta(days=self.sim_time)
            )
        )

        self.sim_results.save_posttransplant(
            date_failure=date_fail,
            date_relist=date_relist,
            cens_date=self.sim_set.SIM_END_DATE,
            matchr_=txp,
            rereg_id=(
                int(
                    self.ptp.offset_ids_transplants +
                    self.ptp.synth_regs
                )
                if date_relist and date_relist < self.sim_set.SIM_END_DATE
                else None
            )
        )

        if (
            cause_fail == cn.PATIENT_RELISTING and
            date_relist and
            date_relist < self.sim_set.SIM_END_DATE
        ):
            # If time-to-event is longer in matched patient than non-matched
            # patient. But then only match to patients with
            # longer time-to-events.
            synth_reg = self.ptp.generate_synthetic_reregistration(
                offer=txp,
                relist_date=date_relist,
                fail_date=date_fail,
                curr_date=(
                    self.sim_start_date +
                    timedelta(days=self.sim_time)
                ),
                verbosity=0
            )

            self.patients[synth_reg.__dict__[cn.ID_REGISTRATION]] = synth_reg
            next_update_time = synth_reg.get_next_update_time()

            if next_update_time:
                self.event_queue.add(
                    Event(
                        type_event=cn.PAT,
                        event_time=next_update_time,
                        identifier=synth_reg.__dict__[cn.ID_REGISTRATION]
                    )
                )
