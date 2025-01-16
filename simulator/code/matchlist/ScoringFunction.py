#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19th

@author: H.C. de Ferrante
"""

from typing import (
    Tuple, Union, Optional, TYPE_CHECKING,
    Any, Dict, Callable, Mapping
)
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.column_names as cn
from simulator.code.utils.utils import round_to_decimals, round_to_int
from collections import defaultdict
import numpy as np


def get_fun_name(x: Any):
    if type(x) is str:
        return x
    else:
        return x.__name__


class InteractionHandler:
    # This class implements an interaction handler.
    def __init__(
            self,
            interactions: Mapping[
                str,
                Tuple[
                    Union[str, Tuple[str, Any]],
                    Union[str, Tuple[str, Any]]
                ]
            ]
    ):
        self.interactions = interactions
        self.interaction_functions = self._precompute_interactions()

    def _precompute_interactions(
        self
    ) -> Dict[str, Callable[[Dict[str, Any]], Any]]:
        interaction_functions = {}
        for inter, intervars in self.interactions.items():
            if len(intervars) == 1:
                interaction_functions[inter] = (
                    lambda d, v=intervars[0]: d[v[0]] == v[1]
                )
            else:
                var1, var2 = intervars
                if isinstance(var1, str):
                    if isinstance(var2, str):
                        interaction_functions[inter] = (
                            lambda d, v1=var1, v2=var2: d[v1] * d[v2]
                        )
                    else:
                        interaction_functions[inter] = (
                            lambda d, v1=var1, v2=var2: (
                                d[v1] * (d[v2[0]] == v2[1])
                            )
                        )
                else:
                    if isinstance(var2, str):
                        interaction_functions[inter] = (
                            lambda d, v1=var1, v2=var2:
                            (d[v1[0]] == v1[1]) * d[v2]
                        )
                    else:
                        interaction_functions[inter] = (
                            lambda d, v1=var1, v2=var2:
                            (d[v1[0]] == v1[1]) * (d[v2[0]] == v2[1])
                        )
        return interaction_functions


class MatchPointFunction:
    """Class which implements a scoring function
    ...

    Attributes   #noqa
    ----------
    coef: dict[str, float]
        coefficients to use to calculate score
    intercept: float
        intercept for calculating score
    trafos: dict[str, str]
        transformations to apply to score component
    caps: dict[str, Tuple[float, float]]
        caps to apply to score component
    limits: Tuple[float, float]
        caps to apply to final score
    round: bool
        whether to round scores to nearest integer

    Methods
    -------
    calc_score() -> float
    """

    def __init__(
            self,
            coef: dict[str, float],
            points_comp_to_group: Dict[str, str],
            trafos: Optional[Dict[str, str]] = {},
            point_multiplier: Optional[Dict[str, Any]] = None,
            caps: Optional[Dict[str, Tuple[float, float]]] = {},
            limits: Optional[Tuple[float, float]] = None,
            clamp_defaults: Optional[Dict[str, int]] = None,
            round: bool = True
    ) -> None:
        self.coef = {k.lower(): fl for k, fl in coef.items()}
        self.variables = set(self.coef.keys())

        # Construct interactions and raw variables,
        # based on inputted coefficient dictionary.
        self.construct_interactions_and_raw_variables(
            self.coef
        )
        self.interaction_handler = InteractionHandler(
            self.coef_interactions
        )

        # Implementation of a gravity model
        if point_multiplier is not None:
            self.point_multiplier = {
                mn.lower(): es.TRAFOS[mn.lower()] for mn in point_multiplier
            }
        else:
            self.point_multiplier = None

        self.caps = caps
        if clamp_defaults:
            self.clamp_defaults = clamp_defaults

        # Transformations to apply. Pre-construct also the identify
        # transformation for untransformed variables.
        if trafos is None:
            self.trafos = {}
        else:
            self.trafos = {
                k.lower(): es.TRAFOS[trafo_name]
                for k, trafo_name in trafos.items()
            }
            for key in self.coef.keys():
                if key not in self.trafos:
                    self.trafos[key] = lambda x: x

        self._vars_to_construct = None
        self._vars_to_construct_initialized = False

        self.limits = limits
        self.round = round
        self.points_comp_to_group = points_comp_to_group

    def construct_interactions_and_raw_variables(
        self,
        coefs: Dict[str, float]
    ) -> None:
        # This function detects interactions in the input dictionary,
        # and initializes these coefficient interactions as a class attribute.
        # It also initializes all variables as a class attribute.
        # For instance, r_ped:mms_hla_dr will be detected as an
        # interaction between r_ped and mms_hla_dr. Both r_ped and mms_hla_dr
        # will be added to the all_raw_variables attribute.

        interactions = {
            key: key.split(':')
            for key in coefs.keys()
            if ':' in key or '-' in key
        }
        self.all_raw_variables = {
            var.split('-')[0] for var in self.coef.keys()
            if var not in interactions.keys()
        }.union(
            {
                var.split('-')[0]
                for var in sum(interactions.values(), [])
            }
        )
        self.coef_interactions = {
            key: list(
                val.split('-') if '-' in val else val for val in values
            ) for key, values in interactions.items()
        }

    def add_interactions(self, d: Dict[str, Any]) -> None:
        d[cn.INTERCEPT] = 1
        for inter, func in (
            self.interaction_handler.interaction_functions.items()
        ):
            try:
                d[inter] = func(d)
            except Exception as e:
                print("Could not add the following interaction:")
                print(inter, func)
                print(d)
                exit()

    def add_multipliers(self, d: Dict[str, Any]) -> None:
        if self.point_multiplier is not None:
            for var, trafo in self.point_multiplier.items():
                d[var] = trafo(d)

    def get_total_multiplier(self, d: Dict[str, Any]) -> float:
        if self.point_multiplier is None:
            return 1
        else:
            return np.prod(
                list(
                    d[variable]
                    for variable in self.point_multiplier.keys()
                )
            )

    def calc_score(
            self,
            match_record: dict[str, Any]
    ) -> float:
        """Calculate the score"""

        self.add_interactions(match_record)
        self.add_multipliers(match_record)

        # Calculate score. This is faster with a for-loop than
        # list comprehension.
        score = 0
        if self.trafos:
            for key, c in self.coef.items():
                score += self.trafos[key](match_record[key]) * c
        else:
            for key, c in self.coef.items():
                score += match_record[key] * c

        if self.point_multiplier is not None:
            total_multiplier = self.get_total_multiplier(match_record)
            score = score * total_multiplier

        if self.round:
            return round_to_int(score)
        return score

    def calc_patient_score(
        self,
        pd: dict[str, Any]
    ) -> float:

        # Calculate score. This is faster with a for-loop than
        # list comprehension.
        score = 0
        for key, coef in self.coef.items():
            if key in self.trafos and key in pd and pd[key]:
                score += self.trafos[key](pd[key]) * coef
            elif key in pd and pd[key]:
                score += pd[key] * coef
        if self.round:
            return round_to_int(score)

        return score

    def calc_score_components(
            self,
            match_record: dict[str, Any]
    ) -> Dict[str, float]:
        """Calculate the score"""

        sc = defaultdict(float)

        self.add_interactions(match_record)
        self.add_multipliers(match_record)

        # Calculate score
        for k, coef in self.coef.items():
            if k not in self.points_comp_to_group:
                raise Exception(f'{k} not in points components')
            if k in self.trafos:
                sc[self.points_comp_to_group[k]] += (
                    self.trafos[k](match_record[k]) * coef
                )
            else:
                sc[self.points_comp_to_group[k]] += match_record[k] * coef

        if self.point_multiplier is not None:
            total_multiplier = self.get_total_multiplier(match_record)
            return {
                k: round_to_int(v * total_multiplier) for k, v in sc.items()
            }
        return {
            k: round_to_int(v) for k, v in sc.items()
        }

    def __str__(self):
        fcoefs = list(
            f'{round_to_decimals(v, p=3)}*'
            f'{get_fun_name(self.trafos.get(k, "I"))}({k})'
            for k, v in self.coef.items()
        )
        if self.point_multiplier is None:
            return ' + '.join(fcoefs)
        else:
            return (
                f'{"*".join(self.point_multiplier.keys()).lower()}'
                f'*({" + ".join(fcoefs)})'
            )

    def set_vars_to_construct(self, in_dict: dict[str, Any]):
        self._vars_to_construct = tuple(
            var for var in self.all_raw_variables
            if var not in in_dict.keys() and var != cn.INTERCEPT
        )
        self._vars_to_construct_initialized = True

    @property
    def vars_to_construct(self):
        return self._vars_to_construct
