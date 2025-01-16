#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

from collections.abc import Mapping
from copy import deepcopy
from typing import Any, Tuple, Set, DefaultDict, FrozenSet, Callable
from math import isnan
import numpy as np
from re import sub
from collections import defaultdict

# Define a type hint for a single key-value pair.
KeyValuePair = Tuple[str, Any]
InnerTuple = Tuple[str, Tuple[KeyValuePair]]
RuleTuple = Tuple[InnerTuple, ...]


def nanOrNone(x):
    if x is None:
        return True
    if isnan(x):
        return True
    return False


def round_to_int(x: float):
    if isnan(x):
        return x
    return int(x + 0.5)


def round_to_decimals(x: float, p: int):
    if isnan(x):
        return x
    p = float(10**p)
    return int(x * p + 0.5) / p


def freeze_set_values(
    in_dict: DefaultDict[str, set]
) -> DefaultDict[str, FrozenSet]:
    out_dict = defaultdict(frozenset)
    for k, v in in_dict.items():
        out_dict[k] = frozenset(v)
    return out_dict


# Transformations
def identity(__x: Any):
    """Identitity transformation function"""
    return __x


def log_ceil(x):
    return np.log(np.ceil(x + 2))


def construct_piecewise_term(
    trafo: str, trafo_x: Callable = identity
) -> Callable:
    """Construct a piecewise term"""

    relation, constant = trafo.split('_')

    assert relation in ['under', 'over'], \
        f'Relation should be "under" or "over", not {relation}'

    constant = sub('(?<=[0-9])p(?=[0-9])', '.', constant)

    if 'log' in constant:
        constant = constant.replace('log', '')
    if 'pc' in constant:
        constant = float(constant.replace('pc', '')) / 100
    else:
        if 'm' in constant:
            constant = -float(''.join(c for c in constant if c.isdigit()))
        else:
            constant = float(constant)

    if relation == 'under':
        def fun(xvals):
            return max(0, trafo_x(constant) - trafo_x(xvals))
        return fun
    else:
        def fun(xvals):
            return max(0, trafo_x(xvals) - trafo_x(constant))
        return fun


class DotDict(dict):
    """Helper class which allows access with dot operator
    """

    def __init__(self, *args, **kwargs):
        super(DotDict, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for key, val in arg.items():
                    self[key] = val

        if kwargs:
            for key, val in kwargs.items():
                self[key] = val

    def __getattr__(self, attr) -> Any:
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(DotDict, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(DotDict, self).__delitem__(key)
        del self.__dict__[key]

    def __deepcopy__(self, memo=None):
        return DotDict(deepcopy(dict(self), memo=memo))
