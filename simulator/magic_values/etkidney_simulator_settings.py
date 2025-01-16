#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

Magic values for the simulator. Only magic values
which are to some extent immodifiable are included.

@author: H.C. de Ferrante
"""

from datetime import timedelta, datetime
from typing import Dict, TYPE_CHECKING, Union
import numpy as np
import simulator.magic_values.magic_values_rules as mgr
import simulator.magic_values.column_names as cn
import simulator.magic_values.column_groups as cg
import simulator.magic_values.magic_values_rules as mr
from simulator.code.utils.utils import (
    identity, log_ceil, log_ceil
)
from functools import reduce
from operator import or_
from math import exp
from itertools import product

if TYPE_CHECKING:
    from simulator.code.matchlist.MatchListETKASandESP import (
        MatchRecordESP, MatchRecordETKAS
    )

DAYS_PER_YEAR_FLOAT = 365.25
DAYS_PER_YEAR = timedelta(days=DAYS_PER_YEAR_FLOAT)

DEFAULT_ETKAS_TIER = 'B'

DEFAULT_ETKAS_ATTR_ORDER = (
    cn.MTCH_TIER,
    cn.TOTAL_MATCH_POINTS,
    cn.MM_TOTAL,
    cn.YEARS_ON_DIAL
)

DEFAULT_ESP_ATTR_ORDER = (
    cn.MTCH_TIER_REV,
    cn.PATIENT_IS_HU,
    cn.TOTAL_MATCH_POINTS
)


def reorder_product(*a):
    for tup in product(*a[::-1]):
        yield tup[::-1]


HLA_MQS = list(
    reorder_product(
        range(3),
        range(3),
        range(3)
    )
)

MISMATCH_STR_DEFINITION = (
    mr.MMB_HLA_A, mr.MMB_HLA_B, mr.MMS_HLA_DR
)

HLA_MQS_STR = list(
    ''.join((str(i) for i in item)) for item in HLA_MQS
)


PROFILE_HLA_MQS = {
    hla: f'profile_hla{mq_str}' for hla, mq_str in zip(HLA_MQS, HLA_MQS_STR)
}

# ET and NHS mismatch definitions
ET_MM_LEVELS = {
    (0, 0, 0): mr.FH,
    (1, 0, 0): mr.DR_plus,
    (0, 1, 0): mr.DR_plus,
    (2, 0, 0): mr.DR_plus,
    (1, 1, 0): mr.DR_plus,
    (0, 2, 0): mr.DR_plus,
    (2, 1, 0): mr.DR0,
    (1, 2, 0): mr.DR0,
    (2, 2, 0): mr.DR0,
    (0, 0, 1): mr.DR1,
    (1, 0, 1): mr.DR1,
    (2, 0, 1): mr.DR1,
    (0, 1, 1): mr.DR1,
    (0, 2, 1): mr.DR1,
    (1, 1, 1): mr.DR1,
    (1, 2, 1): mr.DR1,
    (2, 1, 1): mr.DR1,
    (2, 2, 1): mr.DR1,
    (0, 0, 2): mr.DR2,
    (1, 0, 2): mr.DR2,
    (2, 0, 2): mr.DR2,
    (0, 1, 2): mr.DR2,
    (0, 2, 2): mr.DR2,
    (1, 1, 2): mr.DR2,
    (1, 2, 2): mr.DR2,
    (2, 1, 2): mr.DR2,
    (2, 2, 2): mr.DR2
}

NHS_MM_LEVELS = {
    (0, 0, 0): '1',
    (1, 0, 0): '2',
    (0, 1, 0): '2',
    (1, 1, 0): '2',
    (2, 0, 0): '2',
    (2, 1, 0): '2',
    (0, 0, 1): '2',
    (1, 0, 1): '2',
    (2, 0, 1): '2',
    (0, 2, 0): '3',
    (1, 2, 0): '3',
    (2, 2, 0): '3',
    (0, 1, 1): '3',
    (1, 1, 1): '3',
    (2, 1, 1): '3',
    (0, 2, 1): '4',
    (1, 2, 1): '4',
    (2, 2, 1): '4',
    (0, 0, 2): '4',
    (1, 0, 2): '4',
    (2, 0, 2): '4',
    (0, 1, 2): '4',
    (1, 1, 2): '4',
    (2, 1, 2): '4',
    (0, 2, 2): '4',
    (1, 2, 2): '4',
    (2, 2, 2): '4'
}

YPB_MM_LEVELS = {
    (0, 0, 0): '1',
    (1, 0, 0): '2',
    (2, 0, 0): '2',
    (0, 1, 0): '3',
    (1, 1, 0): '3',
    (2, 1, 0): '3',
    (0, 0, 1): '3',
    (1, 0, 1): '3',
    (2, 0, 1): '3',
    (0, 2, 0): '4',
    (1, 2, 0): '4',
    (2, 2, 0): '4',
    (0, 1, 1): '4',
    (1, 1, 1): '4',
    (2, 1, 1): '4',
    (0, 0, 2): '4',
    (1, 0, 2): '4',
    (2, 0, 2): '4',
    (0, 2, 1): '5',
    (1, 2, 1): '5',
    (2, 2, 1): '5',
    (0, 1, 2): '5',
    (1, 1, 2): '5',
    (2, 1, 2): '5',
    (0, 2, 2): '6',
    (1, 2, 2): '6',
    (2, 2, 2): '7'
}

MATCH_DEF_TO_LEVELS = {
    mgr.NHS_MISMATCH_LEVEL: set(NHS_MM_LEVELS.values()),
    mgr.ET_MISMATCH_LEVEL: set(ET_MM_LEVELS.values()),
    mgr.YPB_MISMATCH_LEVEL: set(YPB_MM_LEVELS.values())
}


def age_difference_filter(__m: Dict):
    """Calculate age difference filter"""
    d_age = __m[cn.D_AGE]
    r_age = __m[cn.R_MATCH_AGE]

    if abs(d_age - r_age) <= 5:
        return 1
    elif r_age > d_age:
        return 1 / exp(0.10 * (r_age - d_age - 5)**0.95)
    else:
        return 1 / exp(0.02 * (d_age - r_age - 5)**0.75)


def age_difference_filter_wrapped(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]'
) -> float:
    """Calculate age difference filter"""
    return age_difference_filter(__m.__dict__)


def age_difference_filter_muted(__m: Dict):
    """Calculate age difference filter"""
    d_age = __m[cn.D_AGE]
    r_age = __m[cn.R_MATCH_AGE]

    if abs(d_age - r_age) <= 5:
        return 1
    elif r_age > d_age:
        return 1 / exp(0.015 * (r_age - d_age - 5)**0.95)
    else:
        return 1 / exp(0.01 * (d_age - r_age - 5)**0.60)


def age_decay_p5_wrapped(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]', ref=50, power=5
):
    """Calculate age difference filter"""
    __x = __m.__dict__[cn.R_MATCH_AGE]
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p6_wrapped(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]', ref=50, power=6
):
    """Calculate age difference filter"""
    __x = __m.__dict__[cn.R_MATCH_AGE]
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p7_wrapped(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]', ref=50, power=7
):
    """Calculate age difference filter"""
    __x = __m.__dict__[cn.R_MATCH_AGE]
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p5_r40_wrapped(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]',
    ref=40, power=5
):
    """Calculate age difference filter"""
    __x = __m.__dict__[cn.R_MATCH_AGE]
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p10_wrapped(
        __m: 'Union[MatchRecordETKAS, MatchRecordESP]',
        ref=50, power=10
):
    """Calculate age difference filter"""
    __x = __m.__dict__[cn.R_MATCH_AGE]
    return 1 / (1 + (__x / ref) ** power)


def construct_nhs_mm_level(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]'
) -> str:
    return NHS_MM_LEVELS[__m.mq]


def construct_ypb_mm_level(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]'
) -> str:
    return YPB_MM_LEVELS[__m.mq]


def construct_et_mm_level(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]'
) -> str:
    return ET_MM_LEVELS[__m.mq]


def construct_decreasing_age_filter_45_75(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]'
) -> float:
    r_age = __m.__dict__[cn.R_MATCH_AGE]
    if r_age <= 45:
        return 1
    elif r_age >= 75:
        return 0
    else:
        return 1 - (r_age - 45) / (75 - 45)


def construct_increasing_age_filter_45_100(
    __m: 'Union[MatchRecordETKAS, MatchRecordESP]'
) -> float:
    r_age = __m.__dict__[cn.R_MATCH_AGE]
    if r_age <= 45:
        return 0
    elif r_age >= 100:
        return 1
    else:
        return (r_age - 45) / (100 - 45)


def exponential_base2(__x: float, base=10):
    """Calculate age difference filter"""
    return (base**__x - 1) / (base - 1)


def exponential_base5(__x: float, base=10):
    """Calculate age difference filter"""
    return (base**__x - 1) / (base - 1)


def exponential_base20(__x: float, base=10):
    """Calculate age difference filter"""
    return (base**__x - 1) / (base - 1)


def exponential_base10(__x: float, base=10):
    """Calculate age difference filter"""
    return (base**__x - 1) / (base - 1)


def exponential_base50(__x: float, base=50):
    """Calculate age difference filter"""
    return (base**__x - 1) / (base - 1)


def exponential_base1000(__x: float, base=1000):
    """Calculate age difference filter"""
    return (base**__x - 1) / (base - 1)


def age_decay_p5(__x: float, ref=50, power=5):
    """Calculate age difference filter"""
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p6(__x: float, ref=50, power=6):
    """Calculate age difference filter"""
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p7(__x: float, ref=50, power=7):
    """Calculate age difference filter"""
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p5_r40(__x: float, ref=40, power=5):
    """Calculate age difference filter"""
    return 1 / (1 + (__x / ref) ** power)


def age_decay_p10(__x: float, ref=50, power=10):
    """Calculate age difference filter"""
    return 1 / (1 + (__x / ref) ** power)


TRAFOS = {
    cn.IDENTITY_TRAFO: identity,
    cn.LOG: np.log,
    cn.AGE_DIFFERENCE_FILTER: age_difference_filter,
    cn.AGE_DIFFERENCE_FILTER_MUTED: age_difference_filter_muted,
    cn.EXPONENTIAL_BASE2: exponential_base2,
    cn.EXPONENTIAL_BASE5: exponential_base5,
    cn.EXPONENTIAL_BASE20: exponential_base20,
    cn.EXPONENTIAL_BASE10: exponential_base10,
    cn.EXPONENTIAL_BASE50: exponential_base50,
    cn.EXPONENTIAL_BASE1000: exponential_base1000,
    cn.AGE_HLA_DECAY_P5: age_decay_p5,
    cn.AGE_HLA_DECAY_P6: age_decay_p6,
    cn.AGE_HLA_DECAY_P7: age_decay_p7,
    cn.AGE_HLA_DECAY_P5_R40: age_decay_p5_r40,
    cn.AGE_HLA_DECAY_P10: age_decay_p10
}

VARIABLE_CONSTRUCTION_FUNCTIONS = {
    mr.ET_MISMATCH_LEVEL: construct_et_mm_level,
    mr.NHS_MISMATCH_LEVEL: construct_nhs_mm_level,
    mr.YPB_MISMATCH_LEVEL: construct_ypb_mm_level,
    cn.LINEAR_AGE_FILTER_45_75: construct_decreasing_age_filter_45_75,
    cn.LINEAR_AGE_FILTER_45_100: construct_increasing_age_filter_45_100,
    cn.AGE_DIFFERENCE_FILTER: age_difference_filter_wrapped,
    cn.AGE_HLA_DECAY_P5: age_decay_p5_wrapped,
    cn.AGE_HLA_DECAY_P6: age_decay_p6_wrapped,
    cn.AGE_HLA_DECAY_P7: age_decay_p7_wrapped,
    cn.AGE_HLA_DECAY_P5_R40: age_decay_p5_r40_wrapped,
    cn.AGE_HLA_DECAY_P10: age_decay_p10_wrapped
}


# Directory with simulation settings
DIR_SIM_SETTINGS = 'simulator/sim_yamls/'
DIR_ACCEPTANCE_COEFS = 'simulator/magic_values/acceptance/'
DIR_POSTTXP_COEFS = 'simulator/magic_values/post_txp/'
DIR_TEST_LR = 'data/test/'

# Path to rescue probs
PATH_RESCUE_COX_BH = (
    'simulator/magic_values/acceptance/bh_rescue_triggered.csv'
)
PATH_RESCUE_COX_COEFS = (
    'simulator/magic_values/acceptance/coefs_rescue_triggered.csv'
)

# Paths to files with travel information
PATH_DRIVING_TIMES = (
    'simulator/magic_values/driving_distances/driving_travel_times.csv'
)

# Paths relevant for the acceptance module
ACCEPTANCE_PATHS = {
    k: DIR_ACCEPTANCE_COEFS + v for k, v in {
        'etkas_rd': 'coefs_recipient_driven_ETKAS.csv',
        'etkas_cd': 'coefs_center_acceptance_ETKAS.csv',
        'esp_rd': 'coefs_recipient_driven_ESP.csv',
        'esp_cd': 'coefs_center_acceptance_ESP.csv',
        'dkt': 'coefs_dkt.csv'
    }.items()
}
LR_TEST_FILES = {
    k: DIR_TEST_LR + v for k, v in {
        'etkas_rd': 'acceptance_pd_ETKAS.csv',
        'etkas_cd': 'acceptance_cd_ETKAS.csv',
        'esp_rd': 'acceptance_pd_ESP.csv',
        'esp_cd': 'acceptance_cd_ESP.csv',
        'dkt': 'test_cases_dkt.csv'
    }.items()
}

# Column mappings for post-transplant survival
POSTTXP_RELIST_STRATA_VARS = {
    cn.R_AGE_GROUP_POSTTXP: cn.R_MATCH_AGE,
    cn.TIME_SINCE_PREV_TXP_CAT: cn.SURV_OFFSET
}

# Paths relevant for post-transplant survival predictions
POSTTXP_RELISTPROB_PATHS = {
    None: DIR_POSTTXP_COEFS + 'prob_relist.csv'
}
POSTTXP_SURV_PATHS = {
    None: DIR_POSTTXP_COEFS + 'posttx_coefs.csv'
}
POSTTXP_SURV_TESTPATHS = {
    None: DIR_TEST_LR + 'posttx_testcases.csv'
}

OFFER_INHERIT_COLS = {
    'donor': [
        cn.D_AGE, cn.DONOR_BALANCE_AGE_GROUP, cn.ESP_DONOR,
        cn.DEATH_CAUSE_GROUP,
        cn.D_TUMOR_HISTORY, cn.D_MARGINAL_FREE_TEXT,
        cn.D_DCD, cn.D_DIABETES, cn.D_HYPERTENSION, cn.D_LAST_CREAT,
        cn.D_CARREST, cn.D_HBSAG, cn.D_HCVAB, cn.D_AGE,
        cn.D_HBCAB, cn.D_DIABETES, cn.D_HYPERTENSION,
        cn.D_LAST_CREAT, cn.D_TUMOR_HISTORY, cn.D_SMOKING,
        cn.DONOR_COUNTRY, cn.DONOR_REGION,
        cn.D_CMV, cn.D_URINE_PROTEIN, cn.D_AGE_GROUP
    ],
    'patient': [
        cn.PATIENT_SEX, cn.URGENCY_CODE, cn.IS_RETRANSPLANT,
        cn.PATIENT_COUNTRY, cn.RECIPIENT_COUNTRY,
        cn.RECIPIENT_REGION, cn.ID_REGISTRATION, cn.AM_STATUS,
        cn.KIDNEY_PROGRAM
    ]
}

# Settings for post-transplant module
POSTTXP_DISCRETE_MATCH_VARS = [cn.REREG_RETURN_DIAL_TIME, cn.RECIPIENT_COUNTRY]
POSTTXP_CONTINUOUS_MATCH_VARS = {
    cn.RETRANSPLANT: [
        cn.AGE_PREV_TXP, cn.TIME_SINCE_PREV_TXP,
        cn.TIME_LIST_TO_REALEVENT, cn.YEARS_ON_DIAL
    ],
    cn.OFFER: [
        cn.R_MATCH_AGE, cn.TIME_TO_REREG,
        cn.TIME_LIST_TO_REALEVENT, cn.YEARS_ON_DIAL
    ]
}
POSTTXP_MIN_MATCHES = 5
POSTTXP_MATCH_CALIPERS = [20.0, 2.0, 1.0, 3]

ACCEPTANCE_CODES = {cn.T1, cn.T3}


POSTTXP_TRANSFORMATIONS = [identity, log_ceil, log_ceil]

POSTTXP_COPY_VARS = [
    cn.PATIENT_COUNTRY, cn.RECIPIENT_CENTER, cn.RECIPIENT_COUNTRY,
    cn.RECIPIENT_REGION, cn.ID_RECIPIENT, cn.R_BLOODGROUP, cn.R_DOB,
    cn.PATIENT_SEX, 'profile', cn.PATIENT_FEMALE, 'hla_string', 'hla'
]

hlas_to_load = (
    mr.HLA_A, mr.HLA_B, mr.HLA_DR, mr.HLA_DQB
)

# Allowed values for various factors.
ALLOWED_BLOODGROUPS = set([mr.BG_A, mr.BG_B, mr.BG_AB, mr.BG_O])
ALLOWED_STATUSES = set([mr.NT, mr.T, mr.HU, mr.IMM, mr.HI])
TERMINAL_STATUSES = frozenset([mgr.R, mgr.D])
EXITING_STATUSES = frozenset([mgr.R, mgr.D, cn.FU])
ACTIVE_STATUSES = set([mr.T, mr.HU, mr.IMM, mr.HI])
ALL_STATUS_CODES = ALLOWED_STATUSES.union(EXITING_STATUSES)

# Countries & specification of country rules.
ET_COUNTRIES = set(
    [
        mgr.NETHERLANDS, mgr.BELGIUM, mgr.GERMANY, mgr.AUSTRIA,
        mgr.HUNGARY, mgr.SLOVENIA, mgr.CROATIA, mgr.LUXEMBOURG
    ]
)

COUNTRIES_REGIONAL_BALANCES = set([mgr.AUSTRIA])

# Priority in extended allocation.
EXTALLOC_INTERNATIONAL_PRIORITY = {
    mgr.AUSTRIA: {
        mgr.GERMANY: 0, mgr.CROATIA: 0, mgr.HUNGARY: 0, mgr.SLOVENIA: 0
    },
    mgr.BELGIUM: {mgr.GERMANY: 0, mgr.NETHERLANDS: 0},
    mgr.CROATIA: {
        mgr.GERMANY: 0, mgr.AUSTRIA: 1, mgr.HUNGARY: 1, mgr.SLOVENIA: 1
    },
    mgr.GERMANY: {mgr.GERMANY: 0},
    mgr.HUNGARY: {
        mgr.SLOVENIA: 1, mgr.GERMANY: 0, mgr.AUSTRIA: 1, mgr.CROATIA: 1
    },
    mgr.NETHERLANDS: {mgr.GERMANY: 1, mgr.BELGIUM: 0},
    mgr.SLOVENIA: {mgr.GERMANY: 0, mgr.AUSTRIA: 1}
}
EXTALLOC_INTERNATIONAL_PRIORITY = {
    cntry: {
        pc: priorities.get(pc, -1) for pc in ET_COUNTRIES
        if pc != cntry
    }
    for cntry, priorities in EXTALLOC_INTERNATIONAL_PRIORITY.items()
}

# Check ETKAS / ESP eligibility every 30 days
CHECK_ETKAS_ESP_ELIGIBILITY = 30
COUNTRIES_ESP_ETKAS_MUTUALLY_EXCLUSIVE = set([mgr.GERMANY])
ESP_CHOICE_ENFORCED = {
    mgr.GERMANY: datetime(year=2011, month=1, day=1)
}

# Regional balance center groups
REG_BAL_CENTER_GROUPS = {
    'AWDTP': mgr.CENTER_VIENNA,
    'AOLTP': mgr.CENTER_UPPERAUSTRIA
}

# DCD accepting countries
DCD_COUNTRIES = set([mgr.NETHERLANDS, mgr.AUSTRIA, mgr.BELGIUM])
DCD_ACCEPTING_COUNTRIES = DCD_COUNTRIES

# Save transplantation probabilities
WINDOWS_TRANSPLANT_PROBS = (
    28, 90, 365
)
PREFIX_SURVIVAL = 'psurv_posttxp'
COLS_TRANSPLANT_PROBS = list(
    f'{PREFIX_SURVIVAL}_{w}' for w in WINDOWS_TRANSPLANT_PROBS
)

HLA_DECAY_FUNS = {
    varfun for varfun in VARIABLE_CONSTRUCTION_FUNCTIONS.keys()
    if 'decay' in varfun
}

HLA_DECAY_INTERACTIONS = {
    f'{mq_def}-{mq_val}:{decay_fun}'
    for mq_def, mq_vals in MATCH_DEF_TO_LEVELS.items()
    for mq_val in mq_vals
    for decay_fun in HLA_DECAY_FUNS
}
HLA_MM_LEVELS_POINTS = {
    f'{mq_def}-{mq_val}'
    for mq_def, mq_vals in MATCH_DEF_TO_LEVELS.items()
    for mq_val in mq_vals
}

POINT_GROUPS = {
    cn.POINTS_HLA: (
        cn.INTERCEPT, mr.MMB_HLA_A, mr.MMB_HLA_B, mr.MMS_HLA_DR
    ),
    cn.POINTS_HLA_PED: (
        f'{var}:{cn.R_PED}' for var in (
            cn.INTERCEPT, mr.MMB_HLA_A, mr.MMB_HLA_B, mr.MMS_HLA_DR)
    ),
    cn.POINTS_HLA_AGE_MATCH: {
        f'nhs_mm_level-{value}:{cn.R_PED}' for value
        in set(
            NHS_MM_LEVELS.values()
        )
    }.union(
        HLA_DECAY_INTERACTIONS
    ).union(
        HLA_MM_LEVELS_POINTS
    ),
    cn.POINTS_AGE: {cn.LINEAR_AGE_FILTER_45_100, cn.AGE_DIFFERENCE_FILTER},
    'nat_bal': (cn.BALANCE_NAT, ),
    'reg_bal': (cn.BALANCE_REG,),
    cn.POINTS_WAIT: (cn.YEARS_ON_DIAL, cn.PREVIOUS_WT),
    cn.POINTS_DIST: (
        cn.MATCH_LOCAL, cn.MATCH_INTERNATIONAL, cn.MATCH_NATIONAL
    ),
    cn.POINTS_PED: (cn.R_PED,),
    cn.POINTS_MMP: (cn.ET_MMP, ),
    cn.POINTS_HAPLOTYPE_FREQUENCY: (
        cn.HLA_MISMATCHFREQ_1ABDR, cn.HLA_MISMATCHFREQ_FH,
        cn.HLA_MISMATCHFREQ_1BDR, cn.HLA_MISMATCHFREQ_0DR
    ),
    cn.POINTS_HU: (cn.PATIENT_IS_HU,),
    cn.POINTS_VPRA: (cn.VPRA,)
}

POINT_COMPONENTS_TO_GROUP = {
    iv: group
    for group, it in POINT_GROUPS.items() for iv in it
}

# Pre-specify columns for simulation results
OUTPUT_COLS_DISCARDS = (
    'reporting_date',
    cn.ID_DONOR, cn.D_WEIGHT, cn.D_DCD,
    'donor_age', cn.D_BLOODGROUP,
    cn.N_OFFERS, cn.N_PROFILE_TURNDOWNS,
    cn.TYPE_OFFER_DETAILED
)
OUTPUT_COLS_PATIENTS = (
    cn.ID_RECIPIENT, cn.ID_REGISTRATION, cn.R_PED,
    cn.RECIPIENT_CENTER, cn.LISTING_DATE,
    cn.EXIT_STATUS, cn.EXIT_DATE, cn.URGENCY_REASON,
    cn.FINAL_REC_URG_AT_TRANSPLANT,
    cn.LAST_NONNT_HU,
    cn.DATE_FIRST_DIAL,
    cn.TIME_SINCE_PREV_TXP,
    cn.TYPE_RETX,
    cn.RECIPIENT_COUNTRY,
    cn.ANY_ACTIVE,
    cn.R_DOB, cn.ANY_HU,
    cn.PATIENT_SEX,
    cn.INIT_URG_CODE,
    cn.R_BLOODGROUP,
    cn.DISEASE_SINCE,
    cn.DISEASE_GROUP,
    cn.URGENCY_CODE,
    cn.VPRA,
    cn.PRA,
    cn.VALID_PRA,
    cn.ET_MMP,
    cn.ET_HLA_MISMATCHFREQ,
    cn.HZ_HLA_A, cn.HZ_HLA_B, cn.HZ_HLA_DR,
    'hla',
    cn.AM_STATUS,
    cn.KIDNEY_PROGRAM
)
OUTPUT_COLS_SNAPSHOT = (
    cn.ID_RECIPIENT, cn.ID_REGISTRATION,
    cn.LISTING_DATE,
    cn.RECIPIENT_COUNTRY,
    cn.RECIPIENT_CENTER, cn.URGENCY_CODE,
    cn.AM_STATUS, cn.KIDNEY_PROGRAM,
    cn.R_BLOODGROUP,
    cn.R_DOB, cn.DATE_FIRST_DIAL,
    cn.TYPE_RETX, cn.TIME_SINCE_PREV_TXP,
    cn.VPRA
)
OUTPUT_COLS_EXITS = (
    (
        cn.ID_RECIPIENT, cn.ID_REGISTRATION, cn.URGENCY_CODE,
        cn.TYPE_RETX, cn.ID_DONOR, cn.DONOR_AGE,
        cn.D_DCD, cn.EXIT_STATUS, cn.URGENCY_REASON,
        cn.LISTING_DATE, cn.EXIT_DATE, cn.MATCH_CRITERIUM,
        cn.MATCH_DISTANCE, cn.GEOGRAPHY_MATCH,
        cn.MATCH_LOCAL, cn.MATCH_NATIONAL, cn.MATCH_INTERNATIONAL,
        cn.PATIENT_REGION, cn.DONOR_REGION,
        cn.MATCH_ABROAD, cn.RECIPIENT_CENTER,
        cn.RECIPIENT_REGION, cn.RECIPIENT_COUNTRY, cn.R_BLOODGROUP,
        cn.D_BLOODGROUP, cn.MATCH_DATE, cn.PATIENT_RANK,
        cn.RANK, cn.D_ALLOC_CENTER, cn.D_ALLOC_REGION,
        cn.D_ALLOC_COUNTRY, cn.R_MATCH_AGE, cn.R_PED,
        cn.PATIENT_IS_HU, cn.PATIENT_SEX,
        cn.ANY_HU, cn.ACCEPTANCE_REASON,
        cn.OFFERED_TO, cn.DISEASE_GROUP, cn.DISEASE_SINCE,
        cn.PROFILE_COMPATIBLE, cn.TYPE_OFFER_DETAILED,
        cn.PROB_ACCEPT_C, cn.PROB_ACCEPT_P, cn.DRAWN_PROB, cn.DRAWN_PROB_C,
        cn.PROB_DKT, cn.DKT,
        cn.TYPE_RECORD, cn.VPRA, cn.ET_MMP,
        cn.ET_HLA_MISMATCHFREQ, cn.HZ_HLA_A, cn.HZ_HLA_B, cn.HZ_HLA_DR,
        cn.EXT_ALLOC_PRIORITY, cn.EXT_ALLOC_TRIGGERED,
        cn.ALLOCATION_MECHANISM, cn.ALLOCATION_PROGRAM,
        cn.PRA, cn.VALID_PRA,
        cn.D_HYPERTENSION, cn.D_LAST_CREAT, cn.D_DIABETES,
        cn.DEATH_CAUSE_GROUP, cn.D_MALIGNANCY,
        cn.YEARS_ON_DIAL, cn.AM_STATUS,
        cn.KIDNEY_PROGRAM, cn.MTCH_TIER
    ) + tuple(cg.MTCH_COLS) + tuple(COLS_TRANSPLANT_PROBS) +
    tuple(MISMATCH_STR_DEFINITION) + tuple(POINT_GROUPS.keys())
)

OUTPUT_COLS_EXIT_CONSTRUCTED = (
    cn.ID_REREGISTRATION, cn.TIME_WAITED
)

MATCH_INFO_COLS = (
    cn.ID_MTR, cn.MATCH_DATE, cn.OFFERED_TO, cn.N_OFFERS, cn.ID_DONOR,
    cn.D_DCD, cn.D_COUNTRY, cn.D_ALLOC_COUNTRY,
    cn.D_ALLOC_REGION, cn.D_ALLOC_CENTER,
    cn.D_WEIGHT, cn.D_AGE,
    cn.TYPE_OFFER, cn.TYPE_OFFER_DETAILED,
    cn.RECIPIENT_CENTER, cn.RECIPIENT_REGION, cn.RECIPIENT_COUNTRY,
    cn.ID_RECIPIENT, cn.ID_REGISTRATION,
    cn.LISTING_DATE, cn.R_BLOODGROUP,
    cn.D_BLOODGROUP, cn.R_MATCH_AGE, cn.R_WEIGHT,
    cn.R_PED, cn.MATCH_CRITERIUM, cn.GEOGRAPHY_MATCH,
    cn.MATCH_DISTANCE, cn.DEATH_CAUSE_GROUP, cn.D_TUMOR_HISTORY,
    cn.PATIENT_SEX, cn.PATIENT_RANK, cn.RANK, cn.ACCEPTED,
    cn.ACCEPTANCE_REASON, cn.PROB_ACCEPT_C, cn.PROB_ACCEPT_P, cn.DRAWN_PROB,
    cn.DRAWN_PROB_C, cn.PROFILE_COMPATIBLE,
    cn.VPRA, cn.ET_MMP, cn.ET_HLA_MISMATCHFREQ,
    cn.HZ_HLA_A, cn.HZ_HLA_B, cn.HZ_HLA_DR,
    cn.PRA, cn.VALID_PRA, cn.TYPE_RECORD
) + tuple(cg.MTCH_COLS) + tuple(POINT_GROUPS.keys())


# Cut-off for transplantation, where candidate's
# reregistration time is returned
CUTOFF_REREG_RETURN_DIAL_TIME = 365

# Add type statuses
STATUS_TYPES = set(
    (
        mgr.URG, mgr.FU, mgr.PRF, mgr.DIAG,
        mgr.HLA, mgr.UNACC, mgr.AM, mgr.DIAL,
        mgr.PRA, mgr.PRG
    )
)
STATUS_TYPES_IGNORE_SYNTH = (
    mgr.PRF, mgr.DIAG, mgr.HLA, mgr.UNACC,
    mgr.AM, mgr.DIAL, mgr.PRG
)

# Event types
EVENT_TYPES = [cn.PAT, cn.DON, cn.BAL, cn.SNAPSHOT]

# Implemented acceptance module policies
PATIENT_ACCEPTANCE_POLICIES = ['Always', 'LR']
CENTER_ACCEPTANCE_POLICIES = ['Always', 'LR']


# Required donor and patient information for match records
MTCH_RCRD_DONOR_COLS = (
    cn.ID_DONOR, cn.DONOR_COUNTRY, cn.D_BLOODGROUP, cn.D_DCD,
    cn.D_AGE, cn.DONOR_BALANCE_AGE_GROUP
)

MTCH_RCRD_PAT_COLS = (
    cn.ID_RECIPIENT, cn.R_BLOODGROUP, cn.RECIPIENT_COUNTRY,
    cn.RECIPIENT_CENTER, cn.PATIENT_IS_HU, cn.PREVIOUS_WT
)

HLA_FORBIDDEN_CHARACTERS = ('\\*', ':', '-')


UNSPLITTABLES = {
    mr.DR3: {mr.DR17, mr.DR18}
}
UNSPLITTABLE_BROADS = set(UNSPLITTABLES.keys())
UNSPLITTABLE_SPLITS = reduce(or_, UNSPLITTABLES.values())

MATCH_TO_BROADS = {
    v: f'mmb_{v}' for v in mr.ALL_HLA_LOCI
}
MATCH_TO_SPLITS = {
    v: f'mms_{v}' for v in mr.ALL_HLA_LOCI
}

REFERENCE = 'reference'


def is_zero_mm(mm_dict):
    try:
        return (
            mm_dict['mmb_hla_a'] + mm_dict['mmb_hla_b'] +
            mm_dict['mms_hla_dr'] == 0
        )
    except Exception as e:
        return None


def is_at_most_1mm(mm_dict):
    try:
        return (
            mm_dict['mmb_hla_a'] + mm_dict['mmb_hla_b'] +
            mm_dict['mms_hla_dr'] <= 1
        )
    except Exception as e:
        return None


def is_at_most_1bdr(mm_dict):
    try:
        return mm_dict['mmb_hla_b'] + mm_dict['mms_hla_dr'] <= 1
    except Exception as e:
        return None


def is_no_dr_mm(mm_dict):
    try:
        return mm_dict['mms_hla_dr'] == 0
    except Exception as e:
        return None


HLA_FAVORABLE_MATCH_DEFINITIONS = {
    cn.HMP_FH: is_zero_mm,
    cn.HMP_1ABDR: is_at_most_1mm,
    cn.HMP_0DR: is_no_dr_mm,
    cn.HMP_1BDR: is_at_most_1bdr
}
