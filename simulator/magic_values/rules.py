#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

Rules for allocation

@author: H.C. de Ferrante
"""


def le_91(_t):
    return _t <= 91


def le_183(_t):
    return _t <= 183


def le_275(_t):
    return _t <= 275


def le_365(_t):
    return _t <= 365


FRACTIONS_RETURN_WAITINGTIME = {
    le_91: 1,
    le_183: 0.75,
    le_275: 0.50,
    le_365: 0.25
}
