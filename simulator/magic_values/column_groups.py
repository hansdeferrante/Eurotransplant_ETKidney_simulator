#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

blood group compatibility rules

@author: H.C. de Ferrante
"""

import sys

sys.path.append('./')
import simulator.magic_values.column_names as cn


MTCH_COLS = (
    cn.TOTAL_MATCH_POINTS, cn.MTCH_TIER
)


ACCEPTANCE_CODES = (
    cn.T1, cn.T3, cn.CR, cn.RR, cn.FP
)


HLA_MISMATCH_FREQS = (
    cn.HLA_MISMATCHFREQ_1ABDR,
    cn.HLA_MISMATCHFREQ_FH,
    cn.HLA_MISMATCHFREQ_1BDR,
    cn.HLA_MISMATCHFREQ_0DR
)
