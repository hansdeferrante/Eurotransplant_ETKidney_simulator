#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

Magic values for the simulator

@author: H.C. de Ferrante
"""

import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es

DEFAULT_DATE_TIME_HMS = '%Y-%m-%d %H:%M:%S'


DTYPE_OFFERLIST = {
    cn.RECIPIENT_OFFERED: 'Int64',
    cn.PATIENT_DOB: 'object',
    cn.ID_MTR: 'Int64',
    cn.ID_DONOR: 'Int64',
    cn.WFTS_MATCH_DATE: 'object',
    cn.ID_TXP: 'Int64',
    cn.D_DCD: 'Int64',
    cn.D_WEIGHT: 'float64',
    cn.D_BLOODGROUP: 'object',
    cn.GRAFT_AGE: 'Int64',
    cn.GRAFT_BLOODGROUP: 'object',
    cn.TYPE_OFFER: 'Int64',
    cn.DONOR_COUNTRY: 'object',
    cn.DONOR_REGION: 'object',
    cn.DONOR_CENTER: 'object',
    cn.D_REGION: 'object',
    cn.D_ALLOC_CENTER: 'object',
    cn.R_MATCH_AGE: 'float64',
    cn.MATCH_CRITERIUM: 'object',
    cn.PATIENT_URGENCY: 'object',
    cn.R_BLOODGROUP: 'object',
    cn.REC_REQ: 'object',
    cn.REQ_INT: 'Int64',
    cn.PATIENT_COUNTRY: 'object',
    cn.PATIENT_REGION: 'object',
    cn.PATIENT_CENTER: 'object',
    cn.RECIPIENT_REGION: 'object',
    cn.R_WEIGHT: 'float64',
    cn.R_HEIGHT: 'float64',
    cn.RNK_OTHER: 'Int64',
    cn.RNK_ACCEPTED: 'Int64',
    cn.ACCEPTED: bool,
    cn.ACTION_CODE: 'object',
    cn.ACTION_CODE_DESC: 'object',
    cn.FINAL_REC_URG_AT_TRANSPLANT: 'object',
    cn.MATCH_COMMENT: 'object',
    cn.PATIENT_HLA_LITERAL: 'str',
    cn.D_HLA_FULL: 'str',
    cn.POINTS_WAIT: 'float64',
    cn.POINTS_DIST: 'Int64',
    cn.POINTS_PED: 'float64',
    cn.POINTS_BALANCE_NAT: 'Int64',
    cn.POINTS_BALANCE_REG: 'float64',
    cn.POINTS_DIST: 'Int64',
    cn.POINTS_HU: 'Int64',
    cn.POINTS_MMP: 'float64',
    cn.POINTS_MMP_SPLIT: 'float64',
    cn.TOTAL_MATCH_POINTS: 'float64',
    cn.UNACC_ANT: 'str',
    cn.MATCH_RANK_LAYER: 'str'
}


DTYPE_PATIENTLIST = {
    cn.ID_RECIPIENT: int,
    cn.ID_REGISTRATION: int,
    cn.TIME_REGISTRATION: 'object',
    cn.TIME_TO_DEREG: 'float64',
    cn.TIME_SINCE_PREV_TXP: 'Int64',
    cn.PREV_TXP_LIVING: 'Int64',
    cn.PREV_TXP_PED: 'Int64',
    cn.KTH_REG: int,
    cn.NTH_TRANSPLANT: int,
    cn.RECIPIENT_COUNTRY: str,
    cn.RECIPIENT_CENTER: str,
    cn.RECIPIENT_REGION: str,
    cn.R_BLOODGROUP: str,
    cn.R_WEIGHT: float,
    cn.R_HEIGHT: float,
    cn.R_DOB: 'object',
    cn.PATIENT_SEX: str,
    cn.R_AGE_LISTING: float,
    cn.PATIENT_HLA: str,
    cn.PREVIOUS_T: 'Int64',
    cn.WLKI_PROGRAMME: str,
    cn.DATE_FIRST_DIAL: 'object'
}

DTYPE_DONORLIST = {
    cn.D_DATE: object,
    cn.ID_DONOR: int,
    # cn.KTH_OFFER: int,
    # cn.TYPE_OFFER: str,
    # cn.TYPE_OFFER_DETAILED: int,
    cn.D_COUNTRY: str,
    cn.D_CENTER: str,
    cn.D_REGION: str,
    cn.D_ESP_SUBREGION: str,
    cn.D_AGE: int,
    # cn.D_WEIGHT: float,
    cn.D_BLOODGROUP: str,
    cn.D_HBSAG: 'Int64',
    cn.D_HCVAB: 'Int64',
    cn.D_HBCAB: 'Int64',
    cn.D_SEPSIS: 'Int64',
    cn.D_MENINGITIS: 'Int64',
    cn.D_MALIGNANCY: 'Int64',
    cn.D_DRUG_ABUSE: 'Int64',
    cn.D_DCD: 'Int64',
    cn.D_EUTHANASIA: 'Int64',
    cn.DONOR_HLA: 'str',
    cn.DEATH_CAUSE_GROUP: str,
    cn.D_TUMOR_HISTORY: 'Int64',
    cn.D_MARGINAL_FREE_TEXT: 'Int64',
    cn.D_SEX: str,
    cn.D_SMOKING: 'Int64',
    cn.N_KIDNEYS_AVAILABLE: 'Int64',
    cn.D_DIABETES: 'Int64',
    cn.D_HYPERTENSION: 'Int64',
    cn.D_SMOKING: 'Int64',
    cn.D_LAST_CREAT: 'float64',
    cn.D_URINE_PROTEIN: 'float64',
    cn.D_CARREST: 'Int64',
    cn.D_RESCUE: 'Int64',
    cn.D_CMV: 'Int64'
}

DTYPE_MATCH_POTENTIALS = {
    **{
        cn.PATIENT_HLA: str
    },
    **{
        hmp_name.lower(): float for hmp_name in
        es.HLA_FAVORABLE_MATCH_DEFINITIONS.keys()
    },
    **{
        hmp_name.lower().replace('hmp_', 'hlamismatchfreq_'): float
        for hmp_name
        in es.HLA_FAVORABLE_MATCH_DEFINITIONS.keys()
    }
}


DTYPE_DONORBALLIST = {
    cn.D_DATE: object,
    cn.D_ALLOC_COUNTRY: str,
    cn.D_ALLOC_CENTER: str,
    cn.D_AGE: int,
    cn.D_BLOODGROUP: str,
    cn.RECIPIENT_CENTER: str,
    cn.RECIPIENT_COUNTRY: str,
    cn.NON_ETKAS_ESP: int
}


DTYPE_DONOR_FILL_NAS = {
    cn.D_HBSAG: 0,
    cn.D_HCVAB: 0,
    cn.D_HBCAB: 0,
    cn.D_SEPSIS: 0,
    cn.D_MENINGITIS: 0,
    cn.D_MALIGNANCY: 0,
    cn.D_DRUG_ABUSE: 0,
    cn.D_DCD: 0,
    cn.D_EUTHANASIA: 0,
    cn.D_RESCUE: 0,
    cn.D_DIABETES: 0
}


DTYPE_PROFILES = {
    cn.ID_REGISTRATION: int,
    cn.TSTART: float,
    cn.PROFILE_MIN_DONOR_AGE: 'Int64',
    cn.PROFILE_MAX_DONOR_AGE: 'Int64',
    cn.PROFILE_DCD: 'Int64',
    cn.PROFILE_HBSAG: 'Int64',
    cn.PROFILE_HCVAB: 'Int64',
    cn.PROFILE_HBCAB: 'Int64',
    cn.PROFILE_SEPSIS: 'Int64',
    cn.PROFILE_MENINGITIS: 'Int64',
    cn.PROFILE_MALIGNANCY: 'Int64',
    cn.PROFILE_DRUG_ABUSE: 'Int64',
    cn.PROFILE_RESCUE: 'Int64',
    cn.PROFILE_EUTHANASIA: 'Int64',
    cn.PROFILE_ESP: 'Int64',
    cn.PROFILE_EXTENDED_ESP: 'Int64'
}
DTYPE_PROFILES.update(
    {
        v: 'Int64' for v in es.PROFILE_HLA_MQS.values()
    }
)

DTYPE_STATUSUPDATES = {
    cn.ID_REGISTRATION: int,
    cn.TSTART: float,
    cn.TYPE_UPDATE: 'category',
    cn.STATUS_VALUE: str,
    cn.STATUS_DETAIL: str,
    cn.STATUS_DETAIL2: str,
    cn.LISTING_DATE: object,
    cn.REMOVAL_REASON: str,
    cn.URGENCY_CODE: 'category'
}

DTYPE_DONORLIST = {
    k: v for
    k, v in DTYPE_DONORLIST.items()
}

DTYPE_OFFERLIST = {
    k.upper(): v for
    k, v in DTYPE_OFFERLIST.items()
}
