#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 16-08-2022

@author: H.C. de Ferrante
"""

import heapq
from simulator.code.events.StatusUpdate import StatusUpdate, ProfileUpdate
import simulator.magic_values.magic_values_rules as mgr
from typing import List, Union, Tuple, Optional
import pandas as pd


class PatientStatusQueue:
    """Class which implements a heapqueue

    Attributes   #noqa
    ----------
    events : List[Union[StatusUpdate, ProfileUpdate]]
        list of status updates
    ...
    """
    def __init__(self, initial_events: List[StatusUpdate]):
        if initial_events:
            self.events = initial_events
            heapq.heapify(self.events)
        else:
            self.events = []

    def remove_status_types(
            self, remove_event_types: Tuple[str]
    ) -> None:
        """Remove certain status types"""
        self.events = [
            e for e in self.events if
            e.type_status not in remove_event_types
        ]
        heapq.heapify(self.events)

    def return_status_types(
        self, event_types: List[str]
    ) -> List[Union[StatusUpdate, ProfileUpdate]]:
        """Return whether status types are scheduled for the patient"""
        if self.events:
            return [
                e for e in self.events if
                e.type_status in event_types
            ]
        return []

    def return_time_to_exit(
        self,
        exit_statuses
    ) -> float:
        """Return time to exiting status"""
        urg_events = self.return_status_types([mgr.URG])
        if urg_events:
            for stat in urg_events:
                if stat.status_value in exit_statuses:
                    return stat.arrival_time
            else:
                return 0
        else:
            return 0

    def truncate_after(self, truncate_time) -> None:
        """Truncate event list to certain time"""
        self.events = [
            e for e in self.events if
            e.arrival_time <= truncate_time
        ]
        heapq.heapify(self.events)

    def add(self, event) -> None:
        """Add an event to the heapqueue"""
        heapq.heappush(self.events, event)

    def next(self) -> Union[StatusUpdate, ProfileUpdate]:
        """Pop event from heapqueue"""
        return heapq.heappop(self.events)

    def first(self) -> Union[StatusUpdate, ProfileUpdate]:
        """Return first event (do not remove it)"""
        return self.events[0]

    def is_empty(self) -> bool:
        """Check if event is empty"""
        return len(self.events) == 0

    def clear(self) -> None:
        """Clear event statuses queue"""
        self.events.clear()

    def __str__(self) -> str:
        string = ''
        sorted_events = sorted(self.events)
        for evnt in sorted_events:
            string += str(evnt) + '\n'
        return string
