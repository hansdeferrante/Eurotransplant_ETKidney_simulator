#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31

@author: H.C. de Ferrante
"""

import simulator.magic_values.etkidney_simulator_settings as es
from typing import Optional


class Event:
    """Class which implements that an event occurs
    ...

    Attributes   #noqa
    ----------
    type_event: str
        type of event which occurs (patient or donor event)
    event_time: float
        time at which the event occurs (relative to the start date)
    identifier: int
        identifier of the subject pertaining to the event (donor or patient id)

    Methods
    -------
    update_event_time():
        updates event time to new event time.

    """
    __slots__ = ['type_event', 'event_time', 'identifier']

    def __init__(
        self,
        type_event: str,
        event_time: float,
        identifier: int
    ) -> None:

        assert type_event in es.EVENT_TYPES, \
            f'type_event must be one of {", ".join(es.EVENT_TYPES)}'
        self.type_event = type_event
        assert event_time is not None, \
            (
                f'event-time cannot be None for {type_event}'
                '-event for {identifier}'
            )
        self.event_time = event_time
        self.identifier = identifier

    def update_event_time(self, new_time: Optional[float]) -> None:
        """Update the event time"""
        assert new_time is not None, \
            (
                f'event-time cannot be None for {self.type_event}'
                f'-event for {self.identifier}'
            )
        assert new_time >= self.event_time, \
            (
                f'new event time needs to be after old event time '
                f'for {self} (new: {new_time})'
            )
        self.event_time = new_time

    def __str__(self):
        return '{0} event for entity {1} at time {2:.2f}'.format(
            self.type_event,
            self.identifier,
            self.event_time
        )

    def __lt__(self, other):
        return self.event_time < other.event_time
