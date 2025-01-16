#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 16-08-2022

@author: H.C. de Ferrante
"""
import heapq
from simulator.code.events.Event import Event
from typing import List


class EventQueue:
    """Class which implements a heapqueue for events

    Attributes   #noqa
    ----------
    events : List[Event]
        list of events

    ...
    """
    def __init__(self):
        self.events = []

    def add(self, event: Event) -> None:
        """Add an event to the heapqueue"""
        if event.event_time is not None:
            heapq.heappush(self.events, event)

    def next(self) -> Event:
        """Pop event from heapqueue"""
        return heapq.heappop(self.events)

    def first(self) -> Event:
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
