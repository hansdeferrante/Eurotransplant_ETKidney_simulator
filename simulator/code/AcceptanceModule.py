#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13

@author: H.C. de Ferrante
"""

from typing import Optional, Tuple, Dict, Union, List, Any
from math import isnan
from statistics import mean
from operator import attrgetter, and_
from collections import defaultdict
from functools import reduce
import numpy as np
import pandas as pd

import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
from simulator.code.utils.utils import (
    round_to_decimals, construct_piecewise_term
)
from simulator.code.matchlist.MatchListETKASandESP import \
    MatchListCurrentETKAS, MatchRecordETKAS, \
    MatchListESP, MatchRecordESP
from simulator.code.matchlist.MatchList import MatchRecord, MatchList
from simulator.code.entities import Donor
from simulator.code.utils.read_input_files import read_rescue_baseline_hazards

import numpy as np


def boolkey_to_str(x):
    """
    Convert a boolean key to a string (0/1) representation.

    Parameters
    ----------
    x : bool or None
        The boolean key to convert.

    Returns
    -------
    str or None
        The string representation of the boolean key,
        or None if the input is None.
    """
    if x is None:
        return None
    return str(int(x)) if type(x) is bool else str(x)


def inv_logit(logit_: float) -> float:
    """
    Calculate the probability from a logit (log-odds).

    Parameters
    ----------
    logit_ : float
        The logit value.

    Returns
    -------
    float
        The probability corresponding to the logit.
    """
    return np.exp(logit_) / (1 + np.exp(logit_))


class AcceptanceModule:
    """
    Class to implement the graft offering module. This simulates
    graft offer acceptance behavior based on a logistic regression

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(
        self,
        seed,
        patient_acc_policy: str,
        center_acc_policy: str,
        dict_paths_coefs: Dict[str, str] = es.ACCEPTANCE_PATHS,
        simulate_rescue: bool = False,
        path_coefs_rescueprobs: Optional[str] = None,
        path_basehaz_rescueprobs: Optional[str] = None,
        verbose: Optional[int] = None,
        max_recipient_driven_offers_per_center: Optional[Dict[str, int]] = {
            mgr.ESP: 10,
            mgr.ETKAS: 5
        },
        simulate_random_effects: bool = True,
        simulate_joint_random_effects: Optional[bool] = False,
        re_variance_components: Optional[Dict[str, Dict[str, float]]] = None
    ):
        # Initialize random number generators
        rng = np.random.default_rng(seed=seed)
        seeds = rng.choice(999999, size=10, replace=False)
        self.rng_center = np.random.default_rng(seed=seeds[1])
        self.rng_rescue = np.random.default_rng(seed=seeds[3])
        self.rng_random_eff = np.random.default_rng(seed=seeds[4])

        self.verbose = verbose if verbose else 0

        # Set patient acceptance policy.
        if patient_acc_policy.lower() == 'LR'.lower():
            self.determine_patient_acceptance = self._patient_accept_lr
        elif patient_acc_policy.lower() == 'Always'.lower():
            self.determine_patient_acceptance = self._patient_accept_always
        else:
            raise ValueError(
                f'Patient acceptance policy should be one of '
                f'{", ".join(es.PATIENT_ACCEPTANCE_POLICIES)}, '
                f'not {patient_acc_policy}'
            )

        # If simulating rescue allocation, also read in
        # fixed effects for Cox model
        self.simulate_rescue = simulate_rescue
        if simulate_rescue:
            if path_basehaz_rescueprobs is None:
                path_basehaz_rescueprobs = es.PATH_RESCUE_COX_BH
            if path_coefs_rescueprobs is None:
                path_coefs_rescueprobs = es.PATH_RESCUE_COX_COEFS

            self.rescue_bh, self.rescue_bh_stratavars = (
                read_rescue_baseline_hazards(path_basehaz_rescueprobs)
            )
            dict_paths_coefs.update(
                {'coxph_rescue': path_coefs_rescueprobs}
            )

        # Initialize coefficients for the logistic regression
        if (
            (patient_acc_policy == 'LR') |
            (center_acc_policy == 'LR')
        ):
            self._initialize_lr_coefs(
                dict_paths_coefs=dict_paths_coefs
            )

        # Set center acceptance policy.
        if center_acc_policy.lower() == 'LR'.lower():
            self.determine_center_acceptance = self._center_accept_lr
        elif center_acc_policy.lower() == 'Always'.lower():
            self.determine_center_acceptance = self._center_accept_always
        else:
            raise ValueError(
                f'Patient acceptance policy should be one of '
                f'{", ".join(es.CENTER_ACCEPTANCE_POLICIES)}, '
                f'not {center_acc_policy}'
            )

        self.calculate_prob_patient_accept = self._calc_prob_accept

        self.simulate_random_effects = simulate_random_effects
        self.simulate_joint_random_effects = simulate_joint_random_effects
        self.max_patient_rejections_per_center = (
            max_recipient_driven_offers_per_center if
            max_recipient_driven_offers_per_center is not None
            else {
                mgr.ESP: 10,
                mgr.ETKAS: 5
            }
        )

        if self.simulate_random_effects:

            if re_variance_components is None:
                raise Exception(
                    f'Cannot simulate random effects, '
                    f'without variance components'
                )
            else:
                self.re_varcomps = re_variance_components

            if self.simulate_joint_random_effects:
                self.joint_re_vars = set(
                    reduce(
                        and_,
                        (sd.keys() for _, sd in self.re_varcomps.items())
                    )
                )
                if self.joint_re_vars:
                    self.joint_res = {
                        joint_re: mean(
                            sd[joint_re]
                            for sd in self.re_varcomps.values()
                        )
                        for joint_re in self.joint_re_vars
                    }
                self.realization_joint_random_effects = {
                    k: defaultdict(dict)
                    for k in self.joint_res
                }
            else:
                self.joint_re_vars = {}
                self.joint_res = None
                self.realization_joint_random_effects = None

            self.realizations_random_effects = {
                k: defaultdict(dict)
                for k in self.re_varcomps.keys()
                if self.joint_res is None or k not in self.joint_res
            }

    def _generate_rescue_eventcurve(
            self, donor: Donor,
            verbose: Optional[int] = 0) -> Tuple[
        np.ndarray,
        np.ndarray
    ]:
        """
        Generate the rescue event curve for a donor.

        Parameters
        ----------
        donor : Donor
            The donor object.
        verbose : int, optional
            Verbosity level.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            The number of offers till rescue and the event probabilities.
        """
        if verbose:
            print('****** generate rescue event curve')
        lp = self._calculate_lp(
            item=donor.__dict__,
            which='coxph_rescue',
            realization_intercept=0,
            verbose=verbose
        )

        cbh = self.rescue_bh
        if self.rescue_bh_stratavars is not None:
            for strata_var in self.rescue_bh_stratavars:
                if cn.N_OFFERS_TILL_RESCUE not in cbh:
                    cbh = cbh[donor.__dict__[strata_var]]
        else:
            cbh = cbh[np.nan]

        ind_cbh = cbh[cn.CBH_RESCUE] * np.exp(lp)
        return cbh[cn.N_OFFERS_TILL_RESCUE], 1 - np.exp(-ind_cbh)

    def simulate_offers_until_nonstandard_alloc(
        self, donor: Donor,
        r_prob: Optional[float] = None,
        stratum: Optional[str] = None,
        verbose: Optional[int] = 0
    ) -> int:
        """
        Simulate the number of offers made until allocation switches
        to non-standard allocation.

        Parameters
        ----------
        donor : Donor
            The donor object.
        r_prob : float, optional
            The rescue probability.
        stratum : str, optional
            The stratum.
        verbose : int, optional
            Verbosity level.

        Returns
        -------
        int
            The number of rejections.
        """
        """ Sample the number of rejections made at triggering rescue/
            extended allocation from the empirical distribution per country
        """

        n_offers, event_probs = self._generate_rescue_eventcurve(
            donor=donor,
            verbose=verbose
        )

        if r_prob is None:
            r_prob = self.rng_rescue.random()

        if all(r_prob <= event_probs):
            return int(0)
        elif any(r_prob <= event_probs):
            which_n_offers = np.argmax(
                event_probs >= r_prob
            )
            kth_offer = n_offers[which_n_offers] - 1
        else:
            which_n_offers = len(n_offers) - 1
            kth_offer = n_offers[which_n_offers]

        if verbose:
            print(f'{kth_offer} due to prob {r_prob}')

        return int(kth_offer)

    def predict_rescue_prob(
            self, donor: Donor, kth_offer: int,
            verbose: Optional[int] = 0) -> float:
        """
        Predict the probabity that rescue is triggered
        after a given number of offers for a specific donor

        Parameters
        ----------
        donor : Donor
            The donor object.
        kth_offer : int
            The offer number.
        verbose : int, optional
            Verbosity level.

        Returns
        -------
        float
            The rescue probability.
        """
        n_offers, event_probs = self._generate_rescue_eventcurve(
            donor=donor,
            verbose=verbose
        )

        if all(n_offers < kth_offer):
            which_prob = len(n_offers) - 1
        else:
            which_prob = np.argmax(
                n_offers > kth_offer
            ) - 1

        return event_probs[which_prob]

    def _patient_accept_always(
        self, match_record: Union[MatchRecord, MatchRecordETKAS]
    ) -> float:
        """
        Policy to always accept the offer if the profile permits.

        Parameters
        ----------
        match_record : Union[MatchRecord, MatchRecordETKAS]
            The match record object.

        Returns
        -------
        float
            Acceptance status.
        """
        """Policy to always accept the offer, if profile permits."""
        if match_record.patient.profile is not None:
            if match_record.offerable:
                match_record.set_acceptance(
                    reason=cn.T3 if match_record.donor.rescue else cn.T1
                )
                return True
        match_record.set_acceptance(reason=cn.FP)
        return False

    def _patient_accept_lr(
        self, match_record: Union[MatchRecord, MatchRecordETKAS]
    ) -> float:
        """
        Simulate whether patient accepts the offer using
        a logistic model.

        Parameters
        ----------
        match_record : Union[MatchRecord, MatchRecordETKAS]
            The match record object.

        Returns
        -------
        float
            Acceptance status.
        """
        """Check acceptance with LR if profile acceptable or not checked."""
        if match_record.offerable:
            if (
                self.calculate_prob_patient_accept(
                    offer=match_record,
                    verbose=self.verbose
                ) >= match_record.patient.get_acceptance_prob()
            ):
                match_record.set_acceptance(
                    reason=cn.T3 if match_record.donor.rescue else cn.T1
                )
                return True
            else:
                match_record.set_acceptance(
                    reason=cn.RR
                )
        else:
            match_record.set_acceptance(
                reason=cn.FP
            )
        return False

    def _center_accept_always(
        self, center_offer: MatchRecord,
        k_previous_center_rejections: int,
        verbose: Optional[int] = None
    ) -> bool:
        """
        Policy to always accept the offer at the center level.

        Parameters
        ----------
        center_offer : MatchRecord
            The center offer object.
        k_previous_center_rejections : int
            Number of previous center rejections.
        verbose : int, optional
            Verbosity level.

        Returns
        -------
        bool
            Acceptance status.
        """
        center_offer.set_acceptance(
            reason=cn.T3 if center_offer.donor.rescue else cn.T1
        )
        return True

    def _center_accept_lr(
        self, center_offer: MatchRecord,
        k_previous_center_rejections: int,
        verbose: Optional[int] = None
    ) -> bool:
        """
        Simulate whether the center accepts using logistic regression.

        Parameters
        ----------
        center_offer : MatchRecord
            The center offer object.
        k_previous_center_rejections : int
            Number of previous center rejections.
        verbose : int, optional
            Verbosity level.

        Returns
        -------
        bool
            Acceptance
        """
        """Check whether the center accepts."""

        ESP = type(center_offer) is MatchRecordESP

        if ESP:
            center_offer.__dict__[cn.K_PREVIOUS_CENTER_REJECTIONS] = (
                str(k_previous_center_rejections)
                if k_previous_center_rejections < 4
                else '4+'
            )
        else:
            center_offer.__dict__[cn.K_PREVIOUS_CENTER_REJECTIONS] = (
                str(k_previous_center_rejections)
                if k_previous_center_rejections < 8
                else '8+'
            )
        center_offer.__dict__[cn.DRAWN_PROB_C] = self.rng_center.random()
        if (
            self.calc_prob_center_willing_to_accept(
                offer=center_offer,
                verbose=verbose,
                selected_model='esp_cd' if ESP else 'etkas_cd'
            ) >= center_offer.__dict__[cn.DRAWN_PROB_C]
        ):
            center_offer.set_acceptance(
                reason=cn.T3 if center_offer.donor.rescue else cn.T1
            )
            return True
        return False

    def simulate_random_intercept(
        self,
        offerdict: Dict[str, Any],
        selected_model: str
    ) -> float:
        """
        Simulate the random intercept for a given model.

        Parameters
        ----------
        offerdict : Dict[str, Any]
            The offer dictionary.
        selected_model : str
            The selected model.

        Returns
        -------
        float
            The simulated random intercept.
        """
        # If simulating joint random effects, take from joint res only.
        realization_intercept = 0
        for var, sd in self.re_varcomps[selected_model].items():
            if var in self.joint_re_vars and self.joint_res is not None:

                if (
                    re := self.realization_joint_random_effects[
                        var
                    ].get(offerdict[var])
                ) is not None:
                    realization_intercept += re
                else:
                    re = self.rng_random_eff.normal(
                        loc=0,
                        scale=self.joint_res[var]
                    )
                    self.realization_joint_random_effects[
                        var
                    ][offerdict[var]] = re
                    realization_intercept += re
            elif (
                re := self.realizations_random_effects[
                    selected_model
                ][var].get(offerdict[var])
            ) is not None:
                realization_intercept += re
            else:
                re = self.rng_random_eff.normal(
                    loc=0,
                    scale=sd
                )
                self.realizations_random_effects[selected_model][
                    var
                ][offerdict[var]] = re
                realization_intercept += re
        return realization_intercept

    def _calc_prob_accept(
        self, offer: MatchRecord, verbose: Optional[int] = None
    ):
        """
        Calculate the probability of acceptance with separate
        models for ETKAS and ESP.

        Parameters
        ----------
        offer : MatchRecord
            The offer object.
        verbose : int, optional
            Verbosity level.

        Returns
        -------
        float
            The probability of acceptance.
        """
        """
        Calculate probability acceptance with
        separate models for ETKAS and ESP"""
        if verbose is None:
            verbose = self.verbose
        if verbose > 1:
            print('*******')

        selected_model = (
            'esp_rd' if type(offer) is MatchRecordESP
            else 'etkas_rd'
        )
        if self.simulate_random_effects:
            realization_intercept = self.simulate_random_intercept(
                offerdict=offer.__dict__,
                selected_model=selected_model
            )
        else:
            realization_intercept = 0

        return self._calculate_logit(
            offer=offer,
            which=selected_model,
            verbose=verbose,
            realization_intercept=realization_intercept
        )

    def calc_prob_dkt(
        self, offer: MatchRecord, verbose: Optional[int] = None
    ):
        """
        Calculate the probability of dual kidney transplantation.

        Parameters
        ----------
        offer : MatchRecord
            The offer object.
        verbose : int, optional
            Verbosity level.

        Returns
        -------
        float
            The probability of dual kidney transplantation.
        """
        """Calculate probability dual kidney transplantation"""
        if verbose is None:
            verbose = self.verbose
        if verbose > 1:
            print('*******')
        return self._calculate_logit(
            offer=offer,
            which='dkt',
            verbose=verbose
        )

    def calc_prob_center_willing_to_accept(
            self, offer: MatchRecord,
            verbose: Optional[int] = 0,
            selected_model: str = 'etkas_cd'
    ) -> float:
        """
        Calculate the probability that the center is willing
        to accept the offer.

        Parameters
        ----------
        offer : MatchRecord
            The offer object.
        verbose : int, optional
            Verbosity level.
        selected_model : str, optional
            The selected model.

        Returns
        -------
        float
            The probability that the center is willing to accept the offer.
        """
        """Calculate probability center accepts"""
        if verbose and verbose > 1:
            print('*******')

        if self.simulate_random_effects:
            realization_intercept = self.simulate_random_intercept(
                offerdict=offer.__dict__,
                selected_model=selected_model
            )
        else:
            realization_intercept = 0
        if verbose:
            print(
                self._calculate_logit(
                    offer=offer,
                    which=selected_model,
                    realization_intercept=realization_intercept
                )
            )

        return self._calculate_logit(
            offer=offer,
            which=selected_model,
            verbose=verbose,
            realization_intercept=realization_intercept
        )

    def _calculate_lp(
            self, item: Dict,
            which: str, verbose: Optional[int] = None,
            realization_intercept: Optional[float] = None
    ):
        """
        Calculate the linear predictor for a given item and model.

        Parameters
        ----------
        item : Dict
            The item dictionary.
        which : str
            The selected model.
        verbose : int, optional
            Verbosity level.
        realization_intercept : float, optional
            The realization intercept.

        Returns
        -------
        float
            The linear predictor.
        """
        # Realization of random intercept
        if realization_intercept:
            lp = realization_intercept
        else:
            lp = 0

        for key, fe_dict in self.fixed_effects[which].items():
            slogit_b4 = lp
            var2 = None
            sel_coefs = None
            if isinstance(fe_dict, dict):
                if es.REFERENCE in fe_dict:
                    var_slope = 0
                    for variable, subgroup in fe_dict.items():
                        if variable == es.REFERENCE:
                            var_slope += subgroup
                        else:
                            var_slope += fe_dict[variable].get(
                                str(item[variable]),
                                0
                            )
                    lp += var_slope * item[key]
                else:
                    # If it is a regular item, it is not a slope.
                    # Simply add the matching coefficient.
                    sel_coefs = fe_dict.get(
                        boolkey_to_str(item[key]),
                        None
                    )
                    if isinstance(sel_coefs, dict):
                        for var2, dict2 in sel_coefs.items():
                            if var2 == es.REFERENCE:
                                lp += dict2
                            else:
                                lp += dict2.get(
                                    boolkey_to_str(item[var2]),
                                    0
                                )
                                val2 = dict2.get(boolkey_to_str(item[var2]), 0)
                                if (verbose > 1) and val2 != 0:
                                    print(
                                        f'{key}-{item[key]}:'
                                        f'{var2}-{item[var2]}: '
                                        f'{val2}'
                                    )
                    elif sel_coefs is None:
                        newkey = (
                            str(int(item[key]))
                            if type(item[key]) is bool
                            else str(item[key])
                        )
                        if self.reference_levels[which][key] is None:
                            fe_dict[newkey] = 0
                            self.reference_levels[which][key] = newkey
                        else:
                            raise Exception(
                                f'Multiple reference levels for {key}:\n'
                                f'\t{self.reference_levels[which][key]} and '
                                f'{newkey}\n are both assumed reference '
                                f'levels.\nExisting keys are:\n'
                                f'{self.fixed_effects[which][key]}'
                            )
                    else:
                        lp += sel_coefs

            elif key == 'intercept':
                lp += fe_dict
            else:
                lp += item[key] * fe_dict

            if (
                (slogit_b4 != lp) & (verbose > 1) and
                (not isinstance(sel_coefs, dict))
            ):
                if key in item:
                    if var2:
                        print(
                            f'{key}-{item[key]}:'
                            f'{var2}-{item[var2]}: '
                            f'{lp - slogit_b4}'
                        )
                    else:
                        print(
                            f'{key}-{item[key]}: '
                            f'{lp - slogit_b4}'
                        )
                else:
                    print(f'{key}: {lp - slogit_b4}')

        for orig_var, fe_dict in (
            self.continuous_transformations[which].items()
        ):
            for coef_to_get, trafo in fe_dict.items():
                if (value := item[orig_var]) is not None:
                    contr = (
                        trafo(value) *
                        self.continuous_effects[which][orig_var][coef_to_get]
                    )
                    lp += contr

                    if (contr != 0) & (verbose > 1):
                        print(f'{coef_to_get}-{value}: {contr}')
                else:
                    print(f'{orig_var} yields None for {item}')

        return lp

    def _initialize_random_effects(self, re_levels: Dict[str, Any]) -> None:
        """
        Initialize the random effects for the model.

        Parameters
        ----------
        re_levels : Dict[str, Any]
            The random effect levels.

        Returns
        -------
        None
        """
        for selected_model, varcomps in self.re_varcomps.items():
            for var, sd in varcomps.items():
                if not (
                    self.simulate_joint_random_effects and
                    var in self.joint_re_vars
                ):
                    if (lvls := re_levels.get(var)):
                        d = dict(
                            zip(
                                lvls,
                                self.rng_random_eff.normal(
                                    loc=0,
                                    scale=sd,
                                    size=len(lvls)
                                )
                            )
                        )
                        self.realizations_random_effects[
                            selected_model
                        ][var].update(d)
                    else:
                        print(f'Random effects for {var} not initialized')

        if self.joint_res:
            for var, sd in self.joint_res.items():
                if (lvls := re_levels.get(var)):
                    d = dict(
                        zip(
                            lvls,
                            self.rng_random_eff.normal(
                                loc=0,
                                scale=sd,
                                size=len(lvls)
                            )
                        )
                    )
                    self.realization_joint_random_effects[var].update(d)
                else:
                    print(f'Random effects for {var} not initialized')

        return

    def _calculate_logit(
            self, offer: MatchRecord,
            which: str, verbose: Optional[int] = None,
            realization_intercept: Optional[float] = None
    ) -> float:
        """
        Calculate the logit for a given offer and model.

        Parameters
        ----------
        offer : MatchRecord
            The offer object.
        which : str
            The selected model.
        verbose : int, optional
            Verbosity level.
        realization_intercept : float, optional
            The realization intercept.

        Returns
        -------
        float
            The logit value.
        """
        """Calculate probability patient accepts"""
        if verbose is None:
            verbose = self.verbose

        slogit = self._calculate_lp(
            item=offer.__dict__,
            which=which,
            verbose=verbose,
            realization_intercept=realization_intercept
        )

        if 'cd' in which:
            offer.__dict__[cn.PROB_ACCEPT_C] = round_to_decimals(
                inv_logit(slogit), 3
            )
        elif 'rd' in which:
            offer.__dict__[cn.PROB_ACCEPT_P] = round_to_decimals(
                inv_logit(slogit), 3
            )
        elif which == 'dkt':
            offer.__dict__[cn.PROB_DKT] = round_to_decimals(
                inv_logit(slogit), 3
            )
        else:
            print(
                f'{which} is not a valid option for prediction '
                f'of acceptance with the ETKidney simulator'
            )
        return inv_logit(slogit)

    def simulate_esp_allocation(
        self, match_list: MatchListESP,
        n_kidneys_available=int
    ) -> Tuple[
        Optional[List[MatchRecord]],
        Dict[str, bool],
        Dict[str, int],
        bool
    ]:
        """
        Simulate the ESP allocation process.

        Parameters
        ----------
        match_list : MatchListESP
            The match list object.
        n_kidneys_available : int
            The number of kidneys available.

        Returns
        -------
        Tuple[
            Optional[List[MatchRecord]],
            Dict[str, bool], Dict[str, int], bool
        ]
            The accepting match records, center willingness to accept,
            rejections per center, and rescue status.
        """

        # Do not simulate rescue in ESP. This was rarely triggered
        # before 2021, and it is not really possible to distinguish
        # between rescue allocation and when rescue was triggered because
        # of an empty / short filtered ESP match list.
        offers_till_rescue = 99999

        # Find accepting match records, center willingness to accept, and
        # whether rescue was triggered for ESP allocation
        acc_mrs, ctr_will, rej_p_ctr_esp, rescue_triggered = (
            self._find_accepting_matchrecords(
                match_list,
                max_offers=offers_till_rescue,
                max_patient_rejections_per_center=(
                    self.max_patient_rejections_per_center[mgr.ESP]
                ),
                n_kidneys_available=n_kidneys_available
            )
        )

        if acc_mrs:
            n_kidneys_accepted = sum(
                2 if mr.__dict__[cn.DKT]
                else 1
                for mr in acc_mrs
            )
        else:
            n_kidneys_accepted = 0

        # If rescue was triggered before an acceptance, continue
        # allocation
        if (
            n_kidneys_accepted < n_kidneys_available
        ):
            match_list.initialize_extalloc_priorities(
                MatchRecordESP
            )
            match_list.donor.rescue = True

            if acc_mrs is not None:
                n_kidneys_available -= n_kidneys_accepted

            # If allocation was unsuccessful on the normal
            # match list, re-try allocation with ESP rescue triggered.
            acc_matchrecords_rescue, _, _, _ = (
                self._find_accepting_matchrecords(
                    match_list,
                    center_willing_to_accept=ctr_will,
                    n_kidneys_available=n_kidneys_available,
                    rej_per_center=rej_p_ctr_esp,
                    max_patient_rejections_per_center=(
                        self.max_patient_rejections_per_center[mgr.ESP]
                    )
                )
            )

            if acc_matchrecords_rescue is not None:
                if acc_mrs is None:
                    return acc_matchrecords_rescue
                else:
                    return acc_mrs + acc_matchrecords_rescue

        if rescue_triggered:
            match_list.donor.rescue = True

        return acc_mrs

    def simulate_etkas_allocation(
        self, match_list: MatchListCurrentETKAS,
        n_kidneys_available: int,
        esp_rescue: bool = False,
        ctr_will: Optional[Dict[str, bool]] = None,
        rej_p_ctr_esp: Optional[Dict[str, int]] = None
    ) -> Optional[List[MatchRecord]]:
        """
        Simulate the ETKAS allocation process.

        Parameters
        ----------
        match_list : MatchListCurrentETKAS
            The match list object.
        n_kidneys_available : int
            The number of kidneys available.
        esp_rescue : bool, optional
            Whether ESP rescue is enabled.
        ctr_will : Dict[str, bool], optional
            Center willingness to accept.
        rej_p_ctr_esp : Dict[str, int], optional
            Rejections per center.

        Returns
        -------
        Optional[List[MatchRecord]]
            The accepting match records.
        """
        """ Iterate over a list of match records. If we model rescue alloc,
            we simulate when rescue will be triggered and terminate
            recipient-driven allocation. In that case, we simulate further
            allocation with the donor identified as a rescue donor, and
            prioritize locally in Belgium / regionally in Germany.
        """
        # When not simulating rescue, offer to all patients on the match list
        if esp_rescue:
            offers_till_rescue = 0
        else:
            if self.simulate_rescue:
                offers_till_rescue = (
                    self.simulate_offers_until_nonstandard_alloc(
                        match_list.donor,
                        stratum=mgr.ETKAS,
                        verbose=0
                    )
                )
            else:
                offers_till_rescue = 9999

        # Allocate until `offers_till_rescue` rescue offers have been made.
        # This returns the match record for the accepting patient
        # (`acc_matchrecord`) if the graft was accepted, and a list of
        # centers willing to accept the organ (determined at the
        # center level)
        acc_mrs, ctr_will, etkas_reg_rej_per_center, _ = (
            self._find_accepting_matchrecords(
                match_list,
                max_offers=offers_till_rescue,
                max_patient_rejections_per_center=(
                    self.max_patient_rejections_per_center[mgr.ETKAS]
                ),
                n_kidneys_available=n_kidneys_available,
                center_willing_to_accept=ctr_will,
                rej_per_center=rej_p_ctr_esp
            )
        )
        if acc_mrs:
            n_kidneys_accepted = sum(
                2 if mr.__dict__[cn.DKT]
                else 1
                for mr in acc_mrs
            )
        else:
            n_kidneys_accepted = 0

        # If rescue was triggered before an acceptance, continue
        # allocation with prioritization for candidates with
        # rescue priority (local in BE, regional in DE).
        if (
            n_kidneys_accepted < n_kidneys_available
        ):
            match_list.initialize_extalloc_priorities(
                MatchRecordETKAS
            )
            match_list.donor.rescue = True

            if acc_mrs is not None:
                n_kidneys_available -= n_kidneys_accepted

            # If rescue is triggered, prioritize remaining waitlist
            # on rescue priority.
            acc_matchrecords_rescue, _, _, _ = (
                self._find_accepting_matchrecords(
                    match_list,
                    center_willing_to_accept=ctr_will,
                    n_kidneys_available=n_kidneys_available,
                    rej_per_center=etkas_reg_rej_per_center,
                    max_patient_rejections_per_center=(
                        self.max_patient_rejections_per_center[mgr.ETKAS]
                    )
                )
            )

            if acc_matchrecords_rescue is not None:
                if acc_mrs is None:
                    return acc_matchrecords_rescue
                else:
                    return acc_mrs + acc_matchrecords_rescue

        return acc_mrs

    def _find_accepting_matchrecords(
            self,
            match_list: Union[MatchListCurrentETKAS, MatchListESP, MatchList],
            n_kidneys_available: int,
            max_offers: Optional[int] = 9999,
            max_patient_rejections_per_center: int = 9999,
            center_willing_to_accept: Optional[Dict[str, bool]] = None,
            rej_per_center: Optional[Dict[str, int]] = None
    ) -> Tuple[
        Optional[List[MatchRecord]],
        Dict[str, bool],
        Dict[str, int],
        bool
    ]:
        """
        Iterate over all match records in the match list and
        simulate if the match object accepts the graft offer.

        Parameters
        ----------
        match_list : Union[MatchListCurrentETKAS, MatchListESP, MatchList]
            The match list object.
        n_kidneys_available : int
            The number of kidneys available.
        max_offers : int, optional
            The maximum number of offers that will be made
            in standard allocation
        max_patient_rejections_per_center : int, optional
            The maximum number of offers per center that counts
            towards triggering rescue allocation.
        center_willing_to_accept : Dict[str, bool], optional
            Dictionary of center names with booleans whether they
            are willing to consider the graft or not.
        rej_per_center : Dict[str, int], optional
            Rejections per center.

        Returns
        -------
        Tuple[
            Optional[List[MatchRecord]],
            Dict[str, bool], Dict[str, int], bool
        ]
            The accepting match records, center willingness to accept,
            rejections per center, and rescue status.
        """
        """ Iterate over all match records in the match list, and simulate
            if the match object accepts the graft offer.

        Parameters
        ------------
        match_list: Union[MatchListCurrentETKAS, MatchListESP, MatchList]
            match list to iterate over, and make offers to
        max_offers: int
            maximum number of offers that will be made
            (rescue allocation is triggered after it)
        max_patient_rejections_per_center: int
            maximum number of offers per center that counts towards
            triggering rescue allocation. Recipients will not be offered
            the graft after this many offers.
        center_willing_to_accept: Optional[Dict[str, bool]]
            Dictionary of center names with booleans whether they are
            willing to consider the graft or not.
        """

        count_rejections_total = 0
        rescue_triggered: bool = False
        n_rejections_per_center = defaultdict(int)
        if center_willing_to_accept is None:
            center_willing_to_accept = {}

        if rej_per_center is not None:
            n_rejections_per_center.update(rej_per_center)

        if match_list == 0:
            return (
                None, center_willing_to_accept,
                n_rejections_per_center, rescue_triggered
            )

        type_match_list = (
            mgr.ETKAS if type(match_list) is MatchListCurrentETKAS
            else mgr.ESP
        )

        # If extended allocation has been triggered, order in offer
        # of extended allocation.
        if (
            match_list.ext_alloc_priority
        ):
            match_records_list = list(
                sorted(
                    match_list.return_match_list(),
                    key=attrgetter(cn.EXT_ALLOC_PRIORITY),
                    reverse=True
                )
            )
        else:
            match_records_list = match_list.return_match_list()

        match_objects: List[MatchRecord] = list()
        n_kidneys_accepted = 0
        # Make offers to match objects. We count an offer as made if
        # (i) it is the first offer to a center,
        # (ii) it is an offer to a filtered match list candidate
        # not among the first 5 offers.
        for match_object in match_records_list:

            object_dict = match_object.__dict__
            rec_center = object_dict[cn.RECIPIENT_CENTER]

            # If more than max number of offers are made, break
            # to initiate rescue.
            if (
                type(max_offers) is int and
                count_rejections_total >= max_offers
            ):
                rescue_triggered = True
                break

            # Skip patients who already rejected the offer.
            # This can happen in case a candidate rejected an
            # offer prior to rescue / extended allocation.
            if (
                cn.ACCEPTANCE_REASON in object_dict and
                object_dict[cn.ACCEPTANCE_REASON] != cn.FP
            ):
                continue

            # Check whether offer is profile compatible.
            if not match_object.donor.rescue and not match_object.offerable:
                match_object.set_acceptance(reason=cn.FP)
                continue

            # For profile compatible offers, initialize
            # acceptance information
            if hasattr(match_object, '_initialize_acceptance_information'):
                match_object._initialize_acceptance_information()

            # 1) determine whether center is willing to accept
            if not (rec_center in center_willing_to_accept):
                if self.determine_center_acceptance(
                    match_object,
                    k_previous_center_rejections=sum(
                        not value
                        for _, value
                        in center_willing_to_accept.items()
                    )
                ):
                    center_willing_to_accept[rec_center] = True
                else:
                    center_willing_to_accept[rec_center] = False
                    count_rejections_total += 1

            # 2) if center finds patient acceptable, make a patient-
            #       driven offer
            if center_willing_to_accept[rec_center]:

                # Skip pediatric candidates if the donor is not pediatric
                # Such acceptances happen almost never
                if (
                    type_match_list == mgr.ETKAS and
                    object_dict[cn.R_PED] and
                    match_object.donor.rescue
                ):
                    match_object.set_acceptance(reason=cn.FP)
                    continue

                # Skip offers to non-ESP aged patients in case
                # extended or rescue was not triggered.
                if (
                    type_match_list == mgr.ESP and
                    not object_dict[cn.PATIENT_ESP_AGED] and
                    not match_object.donor.rescue
                ):
                    continue

                # Make offer if the center has received fewer than
                # `max_patient_rejections_per_center` offers
                if (
                    n_rejections_per_center[rec_center] <
                    max_patient_rejections_per_center
                ):
                    if self.determine_patient_acceptance(match_object):
                        match_objects.append(match_object)
                        # if two kidneys are available, and none were
                        # accepted yet, simulate whether candidate accepts
                        # kidneys en bloc
                        if (
                            n_kidneys_available == 2 and
                            n_kidneys_accepted == 0
                        ):
                            if (
                                self.calc_prob_dkt(
                                    offer=match_object,
                                    verbose=self.verbose
                                ) >
                                match_object.patient.
                                get_acceptance_prob()
                            ):
                                object_dict[cn.DKT] = 1
                                n_kidneys_accepted += 2
                            else:
                                object_dict[cn.DKT] = 0
                                n_kidneys_accepted += 1
                        else:
                            object_dict[cn.DKT] = 0
                            n_kidneys_accepted += 1
                        if n_kidneys_accepted == n_kidneys_available:
                            return (
                                match_objects, center_willing_to_accept,
                                n_rejections_per_center, rescue_triggered
                            )
                        elif n_kidneys_accepted > n_kidneys_available:
                            raise Exception(
                                f'{n_kidneys_accepted} candidate(s)'
                                f'have accepted the kidney, while only '
                                f'{n_kidneys_available} were available.'
                            )
                    elif (
                        object_dict[cn.ACCEPTANCE_REASON] != cn.FP
                    ):
                        n_rejections_per_center[rec_center] += 1
                        count_rejections_total += 1
                elif (
                    n_rejections_per_center[rec_center] >=
                    max_patient_rejections_per_center
                ):
                    # Else record center rejection (CR)
                    match_object.set_acceptance(reason=cn.CR)
                    center_willing_to_accept[rec_center] = False
            else:
                # Else record center rejection (CR)
                match_object.set_acceptance(
                    reason=cn.CR
                )

        # If too little acceptances in regular allocation, return
        # the accepting match record.
        if match_objects:
            return (
                match_objects, center_willing_to_accept,
                n_rejections_per_center, rescue_triggered
            )

        # In case of no acceptance in regular allocation, return
        # center willingness to accept (to trigger rescue)
        return (
            None, center_willing_to_accept,
            n_rejections_per_center, rescue_triggered
        )

    def fe_coefs_to_dict(self, level: pd.Series, coef: pd.Series):
        """
        Process coefficients to dictionary.

        Parameters
        ----------
        level : pd.Series
            The level series.
        coef : pd.Series
            The coefficient series.

        Returns
        -------
        dict
            The processed coefficients dictionary.
        """
        """Process coefficients to dictionary"""
        fe_dict = {}
        for lev, val in zip(level, coef):
            if ':' in str(lev):
                l1, lv2 = lev.split(':')
                try:
                    var2, l2 = lv2.split('-')
                except Exception as e:
                    raise Exception(f'Cannot split {lv2}')
                if not l1:
                    if var2 in fe_dict:
                        fe_dict[var2].update(
                            {l2: val}
                        )
                    else:
                        fe_dict[var2] = {l2: val}
                elif l1 in fe_dict:
                    if var2 in fe_dict[l1]:
                        fe_dict[l1][var2].update(
                            {l2: val}
                        )
                    else:
                        fe_dict[l1].update(
                            {var2: {l2: val}}
                        )
                else:
                    fe_dict.update(
                        {l1: {var2: {l2: val}}}
                    )
            else:
                if (
                    not pd.isnull(level).all() and
                    any(level.str.contains(':', regex=False))
                ):
                    if type(lev) is float and isnan(lev):
                        fe_dict = {es.REFERENCE: val}
                    else:
                        fe_dict[lev] = {es.REFERENCE: val}
                else:
                    fe_dict.update(
                        {lev: val}
                    )
        return fe_dict

    def _initialize_lr_coefs(
        self,
        dict_paths_coefs: Dict[str, str]
    ):
        """
        Initialize the logistic regression coefficients for
        recipient (rd) and center-driven (cd) allocation.

        Parameters
        ----------
        dict_paths_coefs : Dict[str, str]
            Dictionary of paths to coefficient files.

        Returns
        -------
        None
        """
        """Initialize the logistic regression coefficients for
            recipient (rd) and center-driven (cd) allocation
        """

        for k, v in dict_paths_coefs.items():
            self.__dict__[k] = pd.read_csv(v, dtype='str')
            self.__dict__[k]['coef'] = self.__dict__[k]['coef'].astype(float)
            self.__dict__[k].fillna({'coef': 0}, inplace=True)

        coef_keys = list(dict_paths_coefs.keys())

        # Create dictionary for fixed effects
        self.fixed_effects = {}
        self.reference_levels = {}
        for dict_name in coef_keys:
            self.fixed_effects[dict_name] = (
                self.__dict__[dict_name].loc[
                    ~ self.__dict__[dict_name].variable_transformed.notna()
                ].groupby('variable').
                apply(
                    lambda x: self.fe_coefs_to_dict(x['level'], x['coef']),
                    include_groups=False
                ).to_dict()
            )
            self.reference_levels[dict_name] = defaultdict(lambda: None)

        for _, td in self.fixed_effects.items():
            for k, v in td.items():
                if (
                    isinstance(list(v.keys())[0], float) and
                    isnan(list(v.keys())[0])
                ):
                    td[k] = list(v.values())[0]

        self.continuous_transformations = {}
        self.continuous_effects = {}
        for dict_name in coef_keys:
            self.continuous_transformations[dict_name] = (
                self.__dict__[dict_name].loc[
                    self.__dict__[dict_name]
                    .variable_transformed.notna()
                ].groupby('variable').
                apply(
                    lambda x: {
                        k: construct_piecewise_term(v)
                        for k, v in zip(
                            x['variable_transformed'],
                            x['variable_transformation']
                        )
                    },
                    include_groups=False
                ).to_dict()
            )

            self.continuous_effects[dict_name] = (
                self.__dict__[dict_name].loc[
                    self.__dict__[dict_name]
                    .variable_transformed.notna()
                ].groupby('variable').
                apply(
                    lambda x: dict(zip(x['variable_transformed'], x['coef'])),
                    include_groups=False
                ).to_dict()
            )
