#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import unittest

import pandas as pd

sys.path.append("./")

import simulator.code.utils.read_input_files as rdr
import simulator.magic_values.etkidney_simulator_settings as es
from simulator.code.HLA.api import HLAStatsAPI
from simulator.code.HLA.typing import HLATyping


class TestHighResolutionHLAMismatch(unittest.TestCase):
    """Test for high-resolution donor+patient mismatch counting."""

    # Typings from Nizhny Novgorod (see https://github.com/JasonMendoza2008/HLATypingImputationBenchmarks/tree/main/Nizhny_Novgorod_Russia)
    RAW_HIGHRES_TYPINGS = [
        "A*23:01:01 A*23:01:01 B*49:01:01 B*44:03:01 C*07:01:01 C*04:01:01 DRB1*13:02:01 DRB1*07:01:01 DQB1*06:04:01 DQB1*02:02",
        "A*31:01:02:01 A*02:01:01:01 B*27:05:02 B*18:01:01 C*02:02:02 C*07:01:01 DRB1*12:01 DRB1*11:01:01 DQB1*03:01 DQB1*03:01",
        "A*02:01:01:01 A*02:01:01:01 B*51:01:01 B*73:01 C*15:02:01:01 C*15:05:01 DRB1*11:01:01 DRB1*04:05:01 DQB1*03:01 DQB1*02:02",
        "A*02:01:01:01 A*03:01:01:01 B*27:05:02 B*35:01:01 C*15:11 C*04:01:01 DRB1*11:01:01 DRB1*01:01:01 DQB1*03:01 DQB1*05:01",
        "A*02:01:01:01 A*25:01:01 B*27:05:02 B*18:01:01 C*01:02:01 C*12:03:01:01 DRB1*15:01:01 DRB1*03:01:01 DQB1*06:02:01 DQB1*02:01",
        "A*02:01:01:01 A*02:01:01:01 B*35:01:01:02 B*44:02:01 C*03:03:01 C*05:01:01:02 DRB1*11:03 DRB1*08:01:01 DQB1*03:01 DQB1*04:02:01",
        "A*01:01:01:01 A*03:01:01:01 B*08:01:01 B*07:02:01 C*07:01:01 C*07:02:01:03 DRB1*03:01:01 DRB1*15:01:01 DQB1*02:01 DQB1*06:02:01",
        "A*02:01:01:01 A*01:01:01:01 B*40:01:02 B*08:01:01 C*03:04:01:01 C*07:01:01 DRB1*11:01:01 DRB1*03:01:01 DQB1*03:01 DQB1*02:01",
        "A*24:02:01:01 A*24:02:01:01 B*51:01:01 B*39:06:02 C*16:02:01 C*07:02:01:01 DRB1*04:03:01 DRB1*08:01:01 DQB1*03:02 DQB1*04:02:01",
        "A*32:01:01 A*25:01:01 B*40:02:01 B*18:01:01 C*02:02:02 C*12:03:01:01 DRB1*11:01:01 DRB1*04:01:01 DQB1*03:01 DQB1*03:02",
        "A*33:03:01 A*68:01:02:02 B*58:01:01 B*51:01:01 C*03:02:02:01 C*14:02:01 DRB1*03:01:01 DRB1*11:01:01 DQB1*02:01 DQB1*03:01",
        "A*02:01:01:01 A*25:01:01 B*27:05:02 B*18:01:01 C*01:02:01 C*12:03:01:01 DRB1*01:01:01 DRB1*15:01:01 DQB1*05:01 DQB1*06:02:01"
    ]

    EXPECTED_MM_MATRIX = [
    [
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 0,
            "mms_hla_dr": 2,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 0,
            "mms_hla_dr": 2,
            "mm_total": -3
        }
    ],
    [
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 0,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -4
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 0,
            "mms_hla_dr": 2,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 0,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 1,
            "mm_total": -1
        }
    ],
    [
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -2
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -4
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 0,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 0,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        }
    ],
    [
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 0,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 0,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        }
    ],
    [
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 0,
            "mms_hla_dr": 2,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 1,
            "mm_total": -3
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 1,
            "mm_total": -1
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 1,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -5
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 1,
            "mmb_hla_b": 1,
            "mms_hla_dr": 2,
            "mm_total": -4
        },
        {
            "mmb_hla_a": 2,
            "mmb_hla_b": 2,
            "mms_hla_dr": 2,
            "mm_total": -6
        },
        {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0
        }
    ]
]

    @staticmethod
    def _collapse_to_two_field_typing(hla_typing: str) -> str:
        """Reduce each allele token to two-field resolution for match-table lookup."""
        collapsed = []
        for token in hla_typing.split():
            if "*" not in token:
                collapsed.append(token)
                continue
            locus, rest = token.split("*", 1)
            parts = rest.split(":")
            if len(parts) >= 2:
                collapsed.append(f"{locus}*{parts[0]}:{parts[1]}")
            else:
                collapsed.append(token)
        return " ".join(collapsed)

    @classmethod
    def setUpClass(cls) -> None:
        ss = rdr.read_sim_settings(
            os.path.join(es.DIR_SIM_SETTINGS, "sim_settings_test.yaml")
        )
        ss.needed_broad_mismatches = ("hla_b", "hla_a")
        ss.needed_split_mismatches = ("hla_dr",)
        cls.hla_stats_api = HLAStatsAPI(sim_set=ss)

        two_field_typings = [
            cls._collapse_to_two_field_typing(s)
            for s in cls.RAW_HIGHRES_TYPINGS
        ]
        cls.normalized_typings = rdr.fix_hla_string(
            pd.Series(two_field_typings)
        ).tolist()

    def test_count_mismatches_highres_donor_and_patient_matrix(self):
        diagonal_expected = {
            "mmb_hla_a": 0,
            "mmb_hla_b": 0,
            "mms_hla_dr": 0,
            "mm_total": 0,
        }

        for i_d, donor_hla in enumerate(self.normalized_typings):
            for i_p, patient_hla in enumerate(self.normalized_typings):
                d = HLATyping(
                    hla_string=donor_hla,
                    ontology=self.hla_stats_api.ontology,
                )
                p = HLATyping(
                    hla_string=patient_hla,
                    ontology=self.hla_stats_api.ontology,
                )

                mm = self.hla_stats_api.count_mismatches(d_hla=d, p_hla=p)

                self.assertIsNotNone(
                    mm,
                    msg=(
                        f"Mismatch result is None for donor idx {i_d}, "
                        f"patient idx {i_p}"
                    ),
                )
                self.assertEqual(
                    mm,
                    self.EXPECTED_MM_MATRIX[i_d][i_p],
                    msg=(
                        f"Unexpected mismatch result for donor idx {i_d}, "
                        f"patient idx {i_p}"
                    ),
                )

                if i_d == i_p:
                    self.assertEqual(
                        mm,
                        diagonal_expected,
                        msg=f"Diagonal mismatch should be zero at index {i_d}",
                    )


if __name__ == "__main__":
    unittest.main()
