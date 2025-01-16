#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

import simulator.magic_values.etkidney_simulator_settings as es
from typing import Optional
from math import isnan


class StatusUpdate:
    """Class which implements all status updates.
    ...

    Attributes   #noqa
    ----------
    type_status: str
        the debtor of the obligation (i.e., the party who received the
         organ and has to return)
    synthetic: bool
        indicator whether the status update is synthetic
        (i.e. imputed from another patient)
    arrival_time: float
        time at which the status is updated
    status_detail: str
        detailed information regarding status
    status_value: str
        value for status update
    offset_arrival_to_simstart: float
        offset arrival to simulation start (days)
    before_sim_start: bool
        whether arrival is before simulation start
    """
    __slots__ = [
        'type_status', 'synthetic', 'status_detail',
        'status_detail2', 'status_value', 'arrival_time',
        'offset_arrival_to_simstart', 'before_sim_start'
    ]

    def __init__(
        self, type_status: str,
        arrival_time: float,
        sim_start_time: float,
        status_value: str = '',
        status_detail: str = '',
        status_detail2: str = ''
    ) -> None:

        self.type_status, self.synthetic = (type_status, False) \
            if type_status in es.STATUS_TYPES \
            else (type_status.lstrip('S'), True)

        assert self.type_status in es.STATUS_TYPES, \
            f'{self.type_status} is an invalid status type'

        self.status_value = status_value
        self.status_detail = status_detail
        self.status_detail2 = status_detail2

        self.arrival_time = arrival_time

        self.offset_arrival_to_simstart = sim_start_time
        self.before_sim_start = (
            self.offset_arrival_to_simstart + self.arrival_time <= 0
        )

    def __str__(self):
        return '{0} status update at time {1} ({2}: {3})'.format(
            self.type_status,
            self.arrival_time,
            self.status_detail,
            self.status_value
        )

    def __repr__(self):
        return '{0} status update at time {1} ({2}: {3})'.format(
            self.type_status,
            self.arrival_time,
            self.status_detail,
            self.status_value
        )

    def __lt__(self, other):
        return self.arrival_time < other.arrival_time


class ProfileUpdate(StatusUpdate):
    """Class which implements all status updates.
    ...

    Attributes   #noqa
    ----------
    type_status: str
        the debtor of the obligation (i.e., the party who received the
         organ and has to return)
    synthetic: bool
        indicator whether the status update is synthetic
        (i.e. imputed from another patient)
    arrival_time: float
        time at which the status is updated
    status_detail: str
        detailed information regarding status
    status_value: str
        value for status update
    offset_arrival_to_simstart: float
        offset arrival to simulation start (days)
    before_sim_start: bool
        whether arrival is before simulation start
    """
    __slots__ = [
        'type_status', 'synthetic', 'status_detail',
        'status_detail2', 'status_value', 'arrival_time',
        'offset_arrival_to_simstart', 'before_sim_start',
        'profile'
    ]

    def __init__(
            self, type_status: str,
            arrival_time: float,
            sim_start_time: float,
            profile: 'simulator.code.entities.Profile'
    ):
        StatusUpdate.__init__(
            self,
            type_status=type_status,
            arrival_time=arrival_time,
            sim_start_time=sim_start_time
        )
        self.profile = profile
