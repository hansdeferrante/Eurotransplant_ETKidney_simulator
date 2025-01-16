#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

from datetime import timedelta, datetime
from typing import List, Any, Dict, Optional, Type, Tuple, Generator
from itertools import count
from copy import deepcopy
import pandas as pd

import simulator.magic_values.column_names as cn
import simulator.magic_values.column_groups as cg
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
from simulator.code.utils.utils import round_to_int, DotDict, RuleTuple
from simulator.code.matchlist.ScoringFunction import MatchPointFunction
from simulator.code.entities import (
    Patient, Donor, Profile
)
from simulator.code.HLA.HLASystem import HLASystem
from simulator.code.BalanceSystem import BalanceSystem


class MatchDistanceCache:
    """Cache for determining distances between patient and donor"""

    def __init__(self, sim_set: DotDict):
        # Cache to store computed match distances
        # based on (donor_center, recipient_center)
        self._cache = {}

        self.COUNTRIES_WITH_REGIONS = sim_set.COUNTRIES_WITH_REGIONS

    def get_or_compute(
        self, donor_center: str, recipient_center: str,
        match_record: 'MatchRecord'
    ) -> dict:
        # Create a unique key based on the donor and recipient centers
        cache_key = (donor_center, recipient_center)

        # Check if the result is already cached
        if cache_key in self._cache:
            return self._cache[cache_key]

        # If not cached, compute the match distance
        result = self._compute_match_distance(
            match_record
        )

        # Store the result in the cache
        self._cache[cache_key] = result
        return result

    def _compute_match_distance(
        self,
        match_record: 'MatchRecord'
    ) -> Dict[str, Any]:
        """Determine match criterium for allocation"""

        match_record_dict = match_record.__dict__
        d_alloc_country = match_record_dict[cn.D_ALLOC_COUNTRY]
        pat_dict = match_record_dict['patient'].__dict__
        result_dict = {}
        if (
            pat_dict[cn.RECIPIENT_COUNTRY] !=
            d_alloc_country
        ):
            result_dict[cn.GEOGRAPHY_MATCH] = cn.A
            result_dict[cn.MATCH_DISTANCE] = cn.INT
        else:
            # Determine match criteria & geography in case of
            # national allocation
            if (
                pat_dict[cn.RECIPIENT_CENTER] ==
                match_record_dict[cn.D_ALLOC_CENTER]
            ):
                result_dict[cn.GEOGRAPHY_MATCH] = cn.L
            elif match_record.allocation_regional:
                result_dict[cn.GEOGRAPHY_MATCH] = mgr.R
            elif match_record.allocation_national:
                result_dict[cn.GEOGRAPHY_MATCH] = cn.H

            # Determine the match distances
            if d_alloc_country in self.COUNTRIES_WITH_REGIONS:
                if (
                    pat_dict[cn.RECIPIENT_REGION] ==
                    match_record_dict[cn.D_ALLOC_REGION]
                ):
                    result_dict[cn.MATCH_DISTANCE] = cn.L
                else:
                    result_dict[cn.MATCH_DISTANCE] = cn.N
            else:
                result_dict[cn.MATCH_DISTANCE] = cn.L

        # Convert to information necessary for allocation
        result_dict[cn.MATCH_LOCAL] = (
            1 if result_dict[cn.MATCH_DISTANCE] == cn.L else 0
        )
        result_dict[cn.MATCH_NATIONAL] = (
            1 if result_dict[cn.MATCH_DISTANCE] == cn.N else 0
        )
        result_dict[cn.MATCH_INTERNATIONAL] = (
            1 if result_dict[cn.MATCH_DISTANCE] == cn.INT else 0
        )
        return result_dict


class MatchRecord:
    """Class which implements an MatchList
    ...

    Attributes   #noqa
    ----------
    patient: Patient
        Patient
    donor: Donor
        Donor info
    match_date: datetime
        Date that match list is generated.

    Methods
    -------

    """
    match_time: float
    store_score_components: bool
    center_travel_times: Dict[str, Dict[str, Any]]
    total_match_points: float
    date_match: datetime
    sc: Optional[Dict[str, float]]
    donor: Donor
    patient: Patient

    _ids = count(0)

    def __init__(
        self,
        patient: Patient,
        sim_set: DotDict,
        common_info: Dict[Any, Any],
        hla_system: HLASystem,
        calc_points: MatchPointFunction,
        bal_system: BalanceSystem,
        distance_cache: Optional[MatchDistanceCache] = None
    ) -> None:
        "Common info is information copied over to the match record."

        self.distance_cache = distance_cache
        self.sim_set = sim_set

        # Copy over needed information for match record.
        # This includes information that is the same for
        # the entire match list, and information which
        #
        self.__dict__.update(common_info)
        self.__dict__.update(patient.needed_match_info)

        if not self.attr_order_match:
            self.attr_order_match = (cn.TOTAL_MATCH_POINTS,)

        self.patient = patient

        # Add match age and age group
        self.unrounded_match_age = (
            (
                self.patient.age_days_at_listing +
                self.__dict__['match_time'] -
                self.patient.listing_offset
            ) / es.DAYS_PER_YEAR_FLOAT
        )
        self.__dict__[cn.R_MATCH_AGE] = int(
            self.unrounded_match_age
        )

        if self.unrounded_match_age < 18:
            self.__dict__[cn.R_AGE_GROUP] = 'lt18'
        elif self.unrounded_match_age < 50:
            self.__dict__[cn.R_AGE_GROUP] = '18_to_49'
        elif self.unrounded_match_age < 65:
            self.__dict__[cn.R_AGE_GROUP] = '50_to_64'
        else:
            self.__dict__[cn.R_AGE_GROUP] = '65p'

        # Initial attributes
        self._other_profile_compatible = None
        self._mq_compatible = None
        self._no_unacceptables = None
        self._alloc_nat = None
        self._alloc_reg = None
        self._match_tuple = None
        self.total_match_points = 0

        # HLA matches
        try:
            self.__dict__.update(
                hlas := hla_system.count_mismatches(
                    d_hla=self.donor.hla,
                    p_hla=self.patient.hla
                )
            )
            for k in hla_system.loci_zero_mismatch:
                if hlas[k] != 0:
                    self.__dict__[cn.ZERO_MISMATCH] = False
                    break
            else:
                self.__dict__[cn.ZERO_MISMATCH] = True
        except Exception as e:
            raise Exception(
                f'Cannot calculate HLAs for donor {self.donor.id_donor} '
                f'and patient {self.patient.id_recipient}\n'
                f'with donor hla: {self.donor.hla}\n'
                f'with pat.  hla: {self.patient.hla}'
            )

    def add_patient_rank(self, rnk: int) -> None:
        """
        Add patient rank in current sorted list,
        and add it as tie-breaker.
        """
        self.__dict__[cn.PATIENT_RANK] = int(rnk + 1)
        self.attr_order_match += (cn.PATIENT_RANK, )

    def set_acceptance(self, reason: str):
        if reason in cg.ACCEPTANCE_CODES:
            self.__dict__[cn.ACCEPTANCE_REASON] = reason
            if reason == cn.T1 or reason == cn.T3:
                self.__dict__[cn.ACCEPTED] = 1
            else:
                self.__dict__[cn.ACCEPTED] = 0
        else:
            raise ValueError(
                f'{reason} is not a valid acceptance reason.'
            )

    def determine_mismatch_string(self):
        mms = (
            str(
                self.__dict__.get(
                    es.MATCH_TO_SPLITS[loc],
                    self.__dict__.get(es.MATCH_TO_BROADS[loc], ' ')
                )
            ) for loc in mgr.ALL_HLA_LOCI)
        return ''.join(mm for mm in mms if mm != ' ')

    def return_match_info(
        self, cols: Optional[Tuple[str, ...]] = None
    ) -> Dict[str, Any]:
        """Return relevant match information"""
        if cols is None:
            cols = es.MATCH_INFO_COLS

        result_dict = {
            cn.TYPE_RECORD: type(self).__name__
        }

        for key in cols:
            if key in self.__dict__:
                result_dict[key] = self.__dict__[key]
            elif key in self.donor.__dict__:
                result_dict[key] = self.donor.__dict__[key]
            elif key in self.patient.__dict__:
                result_dict[key] = self.patient.__dict__[key]

        return result_dict

    def return_match_tuple(self):
        """Return a match tuple"""
        return tuple(
            self.__dict__[attr] for attr in self.attr_order_match
        )

    def _determine_match_distance(self):
        """ Determine match criterium for allocation.
            Use the DistanceCache if available
        """

        self_dict = self.__dict__
        if self.distance_cache:
            self_dict.update(
                self.distance_cache.get_or_compute(
                    donor_center=self_dict[cn.D_ALLOC_CENTER],
                    recipient_center=self_dict[cn.RECIPIENT_CENTER],
                    match_record=self
                )
            )
            return

        pat_dict = self.patient.__dict__
        d_alloc_country = self_dict[cn.D_ALLOC_COUNTRY]
        if (
            pat_dict[cn.RECIPIENT_COUNTRY] !=
            d_alloc_country
        ):
            self_dict[cn.GEOGRAPHY_MATCH] = cn.A
            self_dict[cn.MATCH_DISTANCE] = cn.INT
        else:
            # Determine match criteria & geography in case of
            # national allocation
            if (
                pat_dict[cn.RECIPIENT_CENTER] ==
                self_dict[cn.D_ALLOC_CENTER]
            ):
                self_dict[cn.GEOGRAPHY_MATCH] = cn.L
            elif self.allocation_regional:
                self_dict[cn.GEOGRAPHY_MATCH] = mgr.R
            elif self.allocation_national:
                self_dict[cn.GEOGRAPHY_MATCH] = cn.H

            # Determine the match distances
            if d_alloc_country in (
                self.sim_set['COUNTRIES_WITH_REGIONS']
            ):
                if (
                    pat_dict[cn.RECIPIENT_REGION] ==
                    self_dict[cn.D_ALLOC_REGION]
                ):
                    self_dict[cn.MATCH_DISTANCE] = cn.L
                else:
                    self_dict[cn.MATCH_DISTANCE] = cn.N
            else:
                self_dict[cn.MATCH_DISTANCE] = cn.L

        # Convert to information necessary for allocation
        self_dict[cn.MATCH_LOCAL] = (
            1 if self_dict[cn.MATCH_DISTANCE] == cn.L else 0
        )
        self_dict[cn.MATCH_NATIONAL] = (
            1 if self_dict[cn.MATCH_DISTANCE] == cn.N else 0
        )
        self_dict[cn.MATCH_INTERNATIONAL] = (
            1 if self_dict[cn.MATCH_DISTANCE] == cn.INT else 0
        )

    def _determine_allocation_national(self):
        """Determine whether allocation is regional"""
        self_dict = self.__dict__
        pat_dict = self.patient.__dict__

        if (
            self_dict[cn.D_ALLOC_COUNTRY] == pat_dict[cn.RECIPIENT_COUNTRY]
        ):
            return True
        else:
            return False

    def _determine_allocation_regional(self):
        """Determine whether allocation is regional"""
        self_dict = self.__dict__
        pat_dict = self.patient.__dict__

        if (
            self.allocation_national and
            (
                self_dict[cn.D_ALLOC_COUNTRY] == mgr.GERMANY or
                self_dict[cn.D_ALLOC_COUNTRY] == mgr.AUSTRIA or
                self_dict[cn.D_ALLOC_COUNTRY] == mgr.BELGIUM
            )
        ):
            if (
                pat_dict[cn.RECIPIENT_REGION] ==
                self_dict[cn.D_ALLOC_REGION]
            ):
                return True
            else:
                return False
        else:
            return False

    def _determine_extalloc_priority(self) -> None:
        """Prioritize local allocation in case of extended allocation."""

        ext_alloc_priorities = es.EXTALLOC_INTERNATIONAL_PRIORITY[
            self.__dict__[cn.DONOR_COUNTRY]
        ]
        max_int_alloc_priority = max(
            ext_alloc_priorities.values()
        )

        if self.__dict__[cn.GEOGRAPHY_MATCH] == cn.L:
            self.__dict__[cn.EXT_ALLOC_PRIORITY] = (
                max_int_alloc_priority + 3
            )
        elif self.allocation_regional:
            self.__dict__[cn.EXT_ALLOC_PRIORITY] = (
                max_int_alloc_priority + 2
            )
        elif self.allocation_national:
            self.__dict__[cn.EXT_ALLOC_PRIORITY] = (
                max_int_alloc_priority + 1
            )
        else:
            self.__dict__[cn.EXT_ALLOC_PRIORITY] = (
                ext_alloc_priorities[
                    self.__dict__[cn.RECIPIENT_COUNTRY]
                ]
            )

    def __repr__(self):
        return (
            f'{self.determine_mismatch_string()} offer to '
            f'{self.__dict__[cn.ID_RECIPIENT]}'
            f'({self.__dict__[cn.RECIPIENT_CENTER]}) '
            f'from {self.__dict__[cn.D_ALLOC_CENTER]} '
            f'with {self.total_match_points} points with '
            f'date first dial: {self.patient.date_first_dial}'
        )

    def __str__(self):
        """Match record"""
        if not self.store_score_components:
            points_str = (
                f'{str(self.total_match_points).rjust(4, " ")} '
                f'total match points'
            )
        else:
            if self.sc:
                points_str = (
                    f'{self.total_match_points}pts:\n\t' + ' '.join(
                        f"{v} {str(k).removeprefix('wfmr_xc_').upper()}".
                        ljust(12)
                        for k, v in self.sc.items()
                    )
                )
        try:
            if self.patient.date_first_dial:
                fdial = self.patient.date_first_dial.strftime("%Y-%m-%d")
        except Exception as e:
            fdial = 'none'
        if cn.MTCH_TIER in self.__dict__:
            tier_str = f' in tier {self.mtch_tier_str.ljust(2)}'
        else:
            tier_str = ''
        accr = (
            self.__dict__.get(cn.ACCEPTANCE_REASON)
            if self.__dict__.get(cn.ACCEPTANCE_REASON)
            else ""
        )
        return (
            f'{self.determine_mismatch_string()} offer{tier_str} to '
            f'{str(self.__dict__[cn.ID_RECIPIENT]).rjust(6, "0")} '
            f'({self.__dict__[cn.RECIPIENT_CENTER]}) '
            f'on {self.date_match.strftime("%Y-%m-%d")} '
            f'from {self.__dict__[cn.D_ALLOC_CENTER]} '
            f'with date first dial: {fdial} '
            f'with {points_str} ({type(self).__name__}, {accr})'
        )

    def __lt__(self, other):
        """Order by match tuple."""
        return self.match_tuple > other.match_tuple

    def __deepcopy__(self, memo):
        # Create a new instance of the class
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        # Copy all attributes, except for hla_system and bal_system
        for k, v in self.__dict__.items():
            if k in {
                'hla_system', 'bal_system',
                'center_travel_times', 'calc_points',
                'mmp_system', 'distance_cache'
            }:
                setattr(result, k, v)  # Shallow copy
            else:
                setattr(result, k, deepcopy(v, memo))  # Deep copy

        return result

    @property
    def match_tuple(self):
        if not self._match_tuple:
            self._match_tuple = self.return_match_tuple()
        return self._match_tuple

    @property
    def allocation_national(self):
        if self._alloc_nat is None:
            self._alloc_nat = self._determine_allocation_national()
        return self._alloc_nat

    @property
    def allocation_regional(self):
        if self._alloc_reg is None:
            self._alloc_reg = self._determine_allocation_regional()
        return self._alloc_reg

    @property
    def other_profile_compatible(self):
        if self._other_profile_compatible is not None:
            return self._other_profile_compatible
        if type(self.patient.profile) is Profile:
            self._other_profile_compatible = (
                self.patient.profile._check_acceptable(
                    self.donor
                )
            )
        else:
            self._other_profile_compatible = True
        self.__dict__[cn.PROFILE_COMPATIBLE] = self._other_profile_compatible
        return self._other_profile_compatible

    @property
    def no_unacceptable_antigens(self):
        if self._no_unacceptables is not None:
            return self._no_unacceptables
        else:
            if (ua := self.patient.unacceptable_antigens):
                self._no_unacceptables = (
                    ua.unacceptables.isdisjoint(
                        self.donor.all_antigens
                    )
                )
            else:
                self._no_unacceptables = True
        return self._no_unacceptables

    def _check_zipped_rules(
            self, attr_name, rule_dict: RuleTuple,
            default=None, value_required: Optional[bool] = False
    ) -> Optional[Any]:
        """Check rules and update attributes accordingly."""
        obj_dict = self.__dict__
        for key, conditions in rule_dict:
            for attr, value in conditions:
                if obj_dict[attr] != value:
                    break
            else:
                obj_dict[attr_name] = key
                if value_required:
                    return key
                else:
                    break
        else:
            obj_dict[attr_name] = default
            if value_required:
                return default

    @property
    def mtch_tier_str(self) -> str:
        return str(self.__dict__[cn.MTCH_TIER])


class MatchList:
    """Class which implements an MatchList
    ...

    Attributes   #noqa
    ----------
    donor: Donor
        Donor that is on offer
    match_date: datetime
        Date that match list is generated.
    match_list: List[MatchRecordETKAS]
        Match list, consisting of list of match records

    Methods
    -------

    """

    _ids = count(0)

    def __init__(
            self, patients: Generator[Patient, None, None], donor: Donor,
            match_date: datetime,
            hla_system: HLASystem,
            bal_system: BalanceSystem,
            calc_points: MatchPointFunction,
            sim_start_date: Optional[datetime],
            type_offer: int = 1,
            alloc_center: Optional[str] = None,
            record_class: Type[MatchRecord] = MatchRecord,
            sort: bool = True,
            store_score_components: bool = False,
            attr_order_match: Optional[Tuple[str]] = None,
            travel_time_dict: Optional[Dict[str, Any]] = None,
            distance_cache: Optional[MatchDistanceCache] = None
    ) -> None:

        self.__dict__[cn.MATCH_DATE] = match_date
        if sim_start_date is not None:
            self.match_time = (
                (match_date - sim_start_date) /
                timedelta(days=1)
            )
        else:
            self.match_time = 0

        self.id_mtr = next(self._ids)
        self.donor = donor
        self.sim_set = donor.sim_set
        if type(match_date) is not datetime:
            raise TypeError(
                f'match_date must be datetime, '
                f'not a {type(match_date)}'
            )

        alloc_center = (
            alloc_center if alloc_center is not None
            else donor.donor_center
        )
        alloc_region = self.sim_set.CENTERS_TO_REGIONS.get(
            alloc_center,
            None
        )
        alloc_country = self.sim_set.FL_CENTER_CODES[alloc_center[0]]
        self.__dict__[cn.D_ALLOC_COUNTRY] = alloc_country
        self.__dict__[cn.TYPE_OFFER_DETAILED] = type_offer

        if (alloc_country == mgr.GERMANY and alloc_region is None):
            raise Exception(
                f'No region defined for German center: {alloc_center}'
            )

        if travel_time_dict:
            self.center_travel_times = travel_time_dict[
                self.donor.__dict__[cn.DONOR_CENTER]
            ]
        else:
            self.center_travel_times = None

        # Info which is common for all match records
        # in this list.
        common_info = donor.needed_match_info
        common_info.update(
            {
                'match_time': self.match_time,
                cn.MATCH_DATE: match_date,
                cn.TYPE_OFFER_DETAILED: type_offer,
                cn.D_ALLOC_CENTER: alloc_center,
                cn.D_ALLOC_REGION: alloc_region,
                cn.D_ALLOC_COUNTRY: alloc_country,
                cn.ID_MTR: self.id_mtr,
                'attr_order_match': attr_order_match,
                'store_score_components': store_score_components,
                'center_travel_times': self.center_travel_times,
                'donor': donor,
            }
        )
        self.match_list = [
            record_class(
                sim_set=self.sim_set,
                patient=pat,
                common_info=common_info,
                calc_points=calc_points,
                hla_system=hla_system,
                bal_system=bal_system,
                distance_cache=distance_cache
            ) for pat in patients
        ]

        if sort:
            self.sorted = True
            self.match_list.sort()
        else:
            self.sorted = False

    def is_empty(self) -> bool:
        """Check if event is empty"""
        return len(self.match_list) == 0

    def return_match_list(
            self
    ) -> List[MatchRecord]:
        if self.sorted:
            return [m for m in self.match_list]
        else:
            raise Exception("Cannot return a match list that is not sorted!")

    def return_match_info(
        self
    ) -> List[Dict[str, Any]]:
        """Return match lists"""
        return [
            matchr.return_match_info() for matchr in self.match_list
        ]

    def initialize_extalloc_priorities(self, exp_class) -> None:
        """Initialize extended allocation priorities."""
        self.ext_alloc_priority = True
        for mr in self.match_list:
            if isinstance(mr, exp_class):
                mr._determine_extalloc_priority()

    def __str__(self) -> str:
        string = ''
        for evnt in sorted(self.match_list):
            string += str(evnt) + '\n'
        return string

    def __repr__(self):
        """Print the match list"""
        string = ''
        for evnt in sorted(self.match_list):
            string += str(evnt) + '\n'
        return string

    def __len__(self):
        return len(self.match_list)
