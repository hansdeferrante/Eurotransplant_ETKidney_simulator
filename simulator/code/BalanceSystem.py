#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

from typing import (
    List, Dict, Tuple, Optional, DefaultDict, Union, Any, Set, Deque
)
from datetime import timedelta
from copy import deepcopy
from collections import defaultdict, deque

import pandas as pd
from statistics import median
from simulator.code.utils.utils import DotDict, round_to_decimals
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
import simulator.magic_values.column_names as cn


class InternationalTransplantation:
    """International transplantation to be balanced
    ...

    Attributes   #noqa
    ----------
    transplant_type: str
        import or export
    transplant_value: str
        import (+1) or export (-1)
    txp_time: float
        time of transplantation
    txp_expiry_time: float
        time at which transplantation does not count any more

    Methods
    -------

    """

    def __init__(
        self,
        group: Tuple[str],
        level: str,
        party: str,
        party_center: str,
        transplant_type: str,
        txp_time: Optional[float] = None,
        txp_expiry_time: Optional[float] = None
    ):
        self.group = group
        self.level = level
        self.type = type
        self.party = party if party != mgr.LUXEMBOURG else mgr.BELGIUM
        if party_center in es.REG_BAL_CENTER_GROUPS:
            self.party_center = es.REG_BAL_CENTER_GROUPS[party_center]
        else:
            self.party_center = party_center
        self.transplant_type = transplant_type
        if transplant_type == mgr.EXPORT:
            self.transplant_value = -1
        elif transplant_type == mgr.IMPORT:
            self.transplant_value = +1
        else:
            raise Exception(
                f'Transplant type should be {mgr.EXPORT} or {mgr.IMPORT}'
                f'not {transplant_type}'
            )
        self.txp_time = txp_time
        self.txp_expiry_time = txp_expiry_time

    def __radd__(self, other):
        if isinstance(other, int):
            return self.transplant_value + other
        else:
            return self.transplant_value + other.transplant_value

    def __add__(self, other):
        return self.transplant_value + other.transplant_value

    def __str__(self):
        return (
            f'{self.group}, {self.transplant_type} for {self.party} '
            f'at {round_to_decimals(self.txp_time, 1)} which expires at '
            f'{round_to_decimals(self.txp_expiry_time, 1)}'
        )

    def __repr__(self):
        return (
            f'{self.group}, {self.transplant_type} for {self.party} '
            f'at {round_to_decimals(self.txp_time, 1)} which expires at '
            f'{round_to_decimals(self.txp_expiry_time, 1)}'
        )

    def __lt__(self, other):
        """Sort by expiry time."""
        return self.txp_expiry_time < other.txp_expiry_time


class BalanceSystem:
    """ Balancing system, potentially with a
        grouping var (e.g. blood type)
    ...

    Attributes   #noqa
    ----------

    Methods
    -------


    """

    def __init__(
        self,
        nations: Set[str],
        initial_national_balance: Dict[
            Tuple[str, ...],
            Dict[str, Deque[InternationalTransplantation]]
        ],
        group_vars: Optional[Set[str]] = None,
        update_balances: bool = True
    ) -> None:

        self._parties = nations
        self.group_vars = group_vars

        self.initial_national_balance = defaultdict(lambda: defaultdict(list))
        self.initial_national_balance.update(initial_national_balance)

        # Ensure that all parties are in the balance dictionary.
        for balance_dict in initial_national_balance.values():
            for party in self._parties:
                if party not in balance_dict:
                    balance_dict[party] = Deque()

        self.national_balances = deepcopy(initial_national_balance)

        self.balance_update_interval = 1
        self.time_last_update = -1000
        self.time_last_update_regbal = -1000
        self._last_nat_bal = None
        self._last_nat_bal_norm = None
        self._last_reg_bal = None
        self._last_reg_bal_norm = None
        self._first_expiry = None
        self.update_balances = update_balances

        self.empty_balances = {
            k: 0 for k in self.parties
        }

    def remove_expired_balance_items(self, current_time: float):
        """Removes expired balances"""
        if self.first_expiry:
            if current_time > self.first_expiry:
                for balances_by_country in self.national_balances.values():
                    for balances in balances_by_country.values():
                        if balances:
                            while (balances and balances[0].txp_expiry_time < current_time):
                                balances.popleft()
                self._first_expiry = None

    def normalize(self, d: Dict[str, int]):
        """ Normalize balances, i.e. substract largest import balance from
            each countries balance
        """
        max_dict = max(d.values())
        return {
            k: max_dict - v for k, v in d.items()
        }

    def return_national_balances(
        self,
        group_values: Tuple[str],
        normalize: bool = False,
        current_time: Optional[float] = 0
    ) -> Dict[str, int]:
        """Return national balance states"""

        # Re-calculate balances daily.
        if current_time:
            if (
                (current_time - self.time_last_update) <=
                self.balance_update_interval
            ):
                if normalize and self._last_nat_bal_norm:
                    return self._last_nat_bal_norm.get(
                        group_values,
                        self.empty_balances
                    )
                elif not normalize and self._last_nat_bal:
                    return self._last_nat_bal.get(
                        group_values,
                        self.empty_balances
                    )
            else:
                self.time_last_update = current_time
                if self.update_balances:
                    self.remove_expired_balance_items(
                        current_time=current_time
                    )

        # Using defaultdict and a single statement to flatten and
        # sum the values
        if self.national_balances:

            national_balances: Dict[str, Dict[str, int]] = (
                defaultdict(lambda: defaultdict(int))
            )

            # This single statement will accumulate sums for each country
            for grp, balances_by_countries in self.national_balances.items():
                for country, bal in balances_by_countries.items():
                    national_balances[grp][country] = sum(bal)

            national_balances[('median',)] = {
                cntry: round(
                    median(
                        national_balances[group][cntry]
                        for group in national_balances
                    )
                )
                for cntry in self.parties
            }
            self._last_nat_bal_norm = {
                k: self.normalize(v) for k, v in national_balances.items()
            }
            self._last_nat_bal = national_balances
        else:
            self._last_nat_bal_norm = defaultdict(
                lambda: {k: 0 for k in self.parties}
            )
            self._last_nat_bal = defaultdict(
                lambda: {k: 0 for k in self.parties}
            )

        if normalize:
            return self._last_nat_bal_norm.get(
                group_values,
                self.empty_balances
            )
        return self._last_nat_bal.get(
            group_values,
            self.empty_balances
        )

    def return_regional_balances(
        self,
        group_values: Tuple[str] = ('1',),
        normalize: bool = False,
        current_time: Optional[Union[float, int]] = None
    ) -> Dict[str, Dict[str, int]]:
        """Return regional balance states"""

        if current_time:
            if (current_time - self.time_last_update_regbal) <= 1:
                if normalize and self._last_reg_bal_norm:
                    return self._last_reg_bal_norm
                elif not normalize and self._last_reg_bal_norm:
                    return self._last_reg_bal_norm
            else:
                self.time_last_update_regbal = current_time

        regional_balances = defaultdict(lambda: defaultdict(int))
        for cntry in es.COUNTRIES_REGIONAL_BALANCES:

            for b in self.national_balances[group_values][cntry]:
                regional_balances[b.party][b.party_center] += (
                    b.transplant_value
                )

        if normalize:
            expected_balances = {
                cntry: (
                    sum(self.national_balances[group_values][cntry]) /
                    len(regional_balances[cntry])
                ) if len(regional_balances[cntry]) > 0 else 0
                for cntry in es.COUNTRIES_REGIONAL_BALANCES
            }
            regional_balances = {
                cntry: {
                    ctr: expected_balances[cntry] - balance
                    for ctr, balance in center_balances.items()
                } for cntry, center_balances in regional_balances.items()
            }

        for cntry, center_dict in regional_balances.items():
            for orig_center, group in es.REG_BAL_CENTER_GROUPS.items():
                if group in center_dict:
                    center_dict.update(
                        {orig_center: center_dict[group]}
                    )

        if normalize:
            self._last_reg_bal_norm = regional_balances
            return self._last_reg_bal_norm
        else:
            self._last_reg_bal = regional_balances
            return self._last_reg_bal

    def __str__(self) -> str:
        """Print the outstanding balances"""
        df_obl = self._current_bal_as_df(which_bal='national')
        return df_obl.to_string()

    def _current_bal_as_df(self, which_bal: str) -> pd.DataFrame:
        """Return current balances as pd.DataFrame"""
        match which_bal:
            case mgr.NATIONAL:
                which_bal = self.national_balances
                group_colnames = self.group_vars
                balkey = 'balance'
                df_obl = pd.DataFrame.from_records(
                    [
                        (group, country, sum(v))
                        for group in which_bal.keys()
                        for country, v in which_bal[group].items()
                    ],
                    columns=['group', 'party', balkey]
                )
                if group_colnames:
                    df_obl[list(group_colnames)] = pd.DataFrame(
                        df_obl['group'].tolist(), index=df_obl.index
                    )
                df_obl = df_obl.drop(columns=['group'])
            case mgr.REGIONAL:
                balkey = 'balance'
                if self.group_vars:
                    group_colnames = self.group_vars + (cn.D_COUNTRY,)
                else:
                    group_colnames = ('group', cn.D_COUNTRY,)
                which_bal = {
                    group + (country,): balances
                    for group in self.national_balances.keys()
                    for country, balances
                    in self.return_regional_balances(
                        group_values=group
                    ).items()
                }
                df_obl = pd.DataFrame.from_records(
                    [
                        (group, center, balance)
                        for group, centerdict in which_bal.items()
                        for center, balance in centerdict.items()
                    ],
                    columns=['group', cn.D_CENTER, balkey]
                )
                if group_colnames:
                    df_obl[list(group_colnames)] = pd.DataFrame(
                        df_obl['group'].tolist(), index=df_obl.index
                    )
                df_obl = df_obl.drop(columns=['group'])
            case _:
                raise Exception(
                    f'which_bal should be {mgr.REGIONAL} or {mgr.NATIONAL}'
                    f' not {which_bal}'
                )
        return df_obl

    def add_balance_from_txp(
            self, txp_time: float, expiry_time: float,
            rcrd: Dict[str, Any]
    ):
        """Add a balance to the balance system"""
        if self.first_expiry is None or expiry_time < self.first_expiry:
            self._first_expiry = None
        self._add_balance_from_record(
            rcrd=rcrd,
            txp_time=txp_time,
            expiry_time=expiry_time,
            target_dict=self.national_balances,
            group_vars=self.group_vars
        )

    @classmethod
    def _add_balance_from_record(
        cls,
        rcrd: Dict[str, Any],
        txp_time: float,
        expiry_time: float,
        target_dict: DefaultDict[
            Tuple[str, ...],
            DefaultDict[str, List[InternationalTransplantation]]
        ],
        group_vars: Optional[Set[str]] = None
    ) -> None:

        if group_vars:
            group = tuple(rcrd[gv] for gv in group_vars)
        else:
            group = ('1',)

        donor_party = (
            rcrd[cn.D_ALLOC_COUNTRY]
            if rcrd[cn.D_ALLOC_COUNTRY] != mgr.LUXEMBOURG
            else mgr.BELGIUM
        )
        rec_party = (
            rcrd[cn.RECIPIENT_COUNTRY]
            if rcrd[cn.RECIPIENT_COUNTRY] != mgr.LUXEMBOURG
            else mgr.BELGIUM
        )

        target_dict[group][donor_party].append(
            InternationalTransplantation(
                group=group,
                level=mgr.NATIONAL,
                party=donor_party,
                party_center=rcrd[cn.D_ALLOC_CENTER],
                transplant_type=mgr.EXPORT,
                txp_time=txp_time,
                txp_expiry_time=expiry_time
            )
        )
        target_dict[group][rec_party].append(
            InternationalTransplantation(
                group=group,
                level=mgr.NATIONAL,
                party=rec_party,
                party_center=rcrd[cn.RECIPIENT_CENTER],
                transplant_type=mgr.IMPORT,
                txp_time=txp_time,
                txp_expiry_time=expiry_time
            )
        )

    @classmethod
    def from_balance_df(
        cls, ss: DotDict,
        df_init_balances: pd.DataFrame,
        group_vars: Optional[Set[str]] = None,
        update_balances: bool = True
    ):
        # Add tstart and tstop columns for existing balances
        df_init_balances.loc[:, cn.TSTART] = (
            df_init_balances.loc[:, cn.D_DATE].values - ss.SIM_START_DATE
        ) / timedelta(days=1)
        df_init_balances.loc[:, cn.TSTOP] = (
            df_init_balances.loc[:, cn.TSTART] + ss.WINDOW_FOR_BALANCE
        )

        init_balances = defaultdict(lambda: defaultdict(deque))

        for rcrd in df_init_balances.to_dict(orient='records'):
            cls._add_balance_from_record(
                rcrd,
                target_dict=init_balances,
                group_vars=group_vars,
                txp_time=rcrd[cn.TSTART],
                expiry_time=rcrd[cn.TSTOP]
            )

        return cls(
            nations=es.ET_COUNTRIES,
            initial_national_balance=init_balances,
            group_vars=group_vars,
            update_balances=update_balances
        )

    @property
    def parties(self) -> Set[str]:
        """Parties involved in the Exception System"""
        return self._parties

    @property
    def first_expiry(self) -> float:
        """Return first expiry time"""
        if self._first_expiry is None:
            if self.national_balances:
                self._first_expiry = min(
                    bal_txps[0]
                    for balances in self.national_balances.values()
                    for bal_txps in balances.values()
                    if len(bal_txps) > 0
                ).txp_expiry_time
        return self._first_expiry
