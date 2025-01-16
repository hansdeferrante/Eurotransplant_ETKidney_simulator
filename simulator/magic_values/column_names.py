#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

blood group compatibility rules

@author: H.C. de Ferrante
"""

import sys
sys.path.append('./')

# Column names for disease groups
DISEASE_GROUP = 'disease_group'
DISEASE_SINCE = 'disease_since'

# List of transplantations
TIME_REGISTRATION = 'time_registration'
DATE_TXP = 'date_txp'
GRAFT_DCD = 'graft_dcd'
TOTAL_MATCH_POINTS = 'total_match_points'
MATCH_COMMENT = 'match_comment'
RECIPIENT_CENTER = 'recipient_center'
RECIPIENT_COUNTRY = 'recipient_country'
R_DOB = 'r_dob'
R_WEIGHT = 'r_weight'
R_HEIGHT = 'r_height'
DONOR_WEIGHT = 'donor_weight'
DONOR_HEIGHT = 'donor_height'
DONOR_AGE = 'donor_age'
DONOR_BALANCE_AGE_GROUP = 'donor_age_group_bal'
DONOR_AGE_75P = 'donor_age_75p'
R_SEX = 'r_sex'
TXP_RESCUE = 'txp_rescue'
ALLOCATION_MECHANISM = 'allocation_mechanism'
ALLOCATION_PROGRAM = 'allocation_program'
MATCH_QUALITY = 'match_quality'
CURRENT_PATIENT_HLA = 'current_patient_hla'
PATIENT_HLA = 'r_hla'
PATIENT_HLA_LITERAL = 'patient_hla'
HZ_HLA_A = 'hz_hla_a'
HZ_HLA_B = 'hz_hla_b'
HZ_HLA_DR = 'hz_hla_dr'
DONOR_HLA = 'donor_hla'
DATE_FIRST_DIAL = 'date_first_dial'
PREVIOUS_WT = 'previous_wt_years'
YEARS_ON_DIAL = 'years_on_dial'
ON_DIAL = 'on_dialysis'
MM_TOTAL = 'mm_total'

# Match list information
RECIPIENT_OFFERED = 'patient_offered'
GEOGRAPHY_MATCH = 'match_geography'
ESP_MATCH_LOCATION = 'esp_match_location'
MATCH_ABROAD = 'match_abroad'
MATCH_DISTANCE = 'distance'
MATCH_LOCAL = 'match_l'
MATCH_NATIONAL = 'match_n'
MATCH_INTERNATIONAL = 'match_i'
R_PED = 'r_pediatric'
D_PED = 'd_pediatric'
ID_MTR = 'id_mtr'
MATCH_DATE = 'date_match'
WFTS_MATCH_DATE = 'wfts_match_date'
ID_TXP = 'id_txp'
N_OFFERS = 'n_offers'
N_PROFILE_TURNDOWNS = 'n_profile_turndowns'

# Match geographies for ESP program
MATCH_LSR = 'match_lsr'
MATCH_SUBREGIONAL = 'match_subregional'

MATCH_RANK_LAYER = 'patient_match_rank_layer'
POINTS_WAIT = 'wfmr_xc_wait'
POINTS_DIST = 'wfmr_xc_dist'
POINTS_PED = 'wfmr_xc_bonus_paed'
POINTS_BALANCE_NAT = 'wfmr_xc_balance'
POINTS_BALANCE_REG = 'wfmr_xc_balance_reg'
NON_ETKAS_ESP = 'nonetkasesp'
POINTS_DIST = 'wfmr_xc_dist'
POINTS_HU = 'wfmr_xc_bonus_hu'
POINTS_MMP_SPLIT = 'wfmr_xc_mmp_split'
POINTS_MMP = 'wfmr_xc_mmp'
POINTS_VPRA = 'points_vpra'
POINTS_HAPLOTYPE_FREQUENCY = 'points_haplotype_freq'

BALANCE_REG = 'bal_reg'
BALANCE_NAT = 'bal_nat'

ANY_UNACC = 'any_unacceptables'
VPRA = 'vpra'
VPRA_PERCENT = 'vpra_percent'
PRA = 'pra'
PRG = 'prg'
VALID_PRA = 'valid_pra'

# Mismatch probabilities & mismatch frequencies.
# Columns starting with ET_ are calculated under
# the assumption that HLAs are independently distributed
# according to HLA match frequency tables. Those
# without prefix are calculated assuming that BGs,
# vPRAs and HLAs are independent (but not HLA loci).
# Those with count suffixes are counted in the donor pool.
ET_MMP = 'et_mmp'
ET_MMP_BROAD = 'et_mmp_broad'
ET_MMP_SPLIT = 'et_mmp_split'
ET_HLA_MISMATCHFREQ = 'et_hla_mismatchfreq'
HLA_MISMATCHFREQ = 'hla_mismatchfreq'

HLA_MISMATCHFREQ_1ABDR = 'hlamismatchfreq_1abdr'
HLA_MISMATCHFREQ_FH = 'hlamismatchfreq_fh'
HLA_MISMATCHFREQ_1BDR = 'hlamismatchfreq_1bdr'
HLA_MISMATCHFREQ_0DR = 'hlamismatchfreq_0dr'

WLKI_PROGRAMME = 'wlki_programme'
KIDNEY_PROGRAM = 'kidney_program'
UNACC_ANT = 'unacceptable_antigens'

D_FULLY_HOMOZYGOUS = 'd_fully_homozygous'
R_HOMOZYGOSITY_LEVEL = 'r_homozygosity_level'
ZERO_MISMATCH = 'zero_mismatch'

ID_DONOR = 'id_donor'
N_KIDNEYS_AVAILABLE = 'n_kidneys_available'
D_DCD = 'graft_dcd'
D_WEIGHT = 'd_weight'
D_AGE = 'donor_age'
D_AGE_GROUP = 'donor_age_group'

ESP_DONOR = 'esp_donor'
GRAFT_AGE = 'graft_age'
D_HLA_FULL = 'd_hla_full'
D_BLOODGROUP = 'd_bloodgroup'
GRAFT_BLOODGROUP = 'graft_bloodgroup'
TYPE_OFFER = 'type_offer'
TYPE_OFFER_DETAILED = 'type_offer_detailed'
D_COUNTRY = 'd_country'
DONOR_COUNTRY = 'donor_country'
DONOR_CENTER = 'donor_center'
DONOR_REGION = 'donor_region'
D_REGION = 'd_region'
D_ESP_SUBREGION = 'd_center_subregion_esp'
D_CENTER = 'd_center'
D_ALLOC_COUNTRY = 'donor_alloc_country'
D_ALLOC_REGION = 'donor_alloc_region'
D_ALLOC_CENTER = 'donor_alloc_center'
R_DOB = 'r_dob'
PATIENT_DOB = 'patient_dob'
R_MATCH_AGE = 'r_match_age'
PATIENT_ESP_AGED = 'patient_esp_aged'
R_AGE_GROUP = 'r_age_group'
R_AGE_GROUP_POSTTXP = 'r_age_group_posttxp'
MATCH_CRITERIUM = 'match_criterium'
EXT_ALLOC_PRIORITY = 'rescue_priority'
EXT_ALLOC_TRIGGERED = 'ext_alloc_triggered'
TYPE_RECORD = 'type_record'
URGENCY_CODE = 'urgency_code'
PATIENT_URGENCY = 'patient_urgency'
INIT_URG_CODE = 'init_urg_code'
ANY_ACTIVE = 'any_active'
ANY_HU = 'any_hu'
R_BLOODGROUP = 'r_bloodgroup'
REC_REQ = 'rec_req'
REQ_INT = 'req_int'
RECIPIENT_COUNTRY = 'recipient_country'
RECIPIENT_REGION = 'recipient_region'
RECIPIENT_CENTER = 'recipient_center'
PATIENT_COUNTRY = 'patient_country'
PATIENT_REGION = 'patient_region'
PATIENT_CENTER = 'patient_center'
R_WEIGHT = 'r_weight'
R_HEIGHT = 'r_height'
RNK_OTHER = 'rnk_other'
RNK_ACCEPTED = 'rnk_accepted'
ACCEPTED = 'accepted'
ACCEPTANCE_REASON = 'acceptance_reason'
K_PREVIOUS_CENTER_REJECTIONS = 'k_previous_center_rejections'
PROB_ACCEPT_P = 'prob_accept_p'
PROB_DKT = 'prob_dkt'
DKT = 'dkt'
DRAWN_PROB = 'drawn_prob_p'
DRAWN_PROB_C = 'drawn_prob_c'
PROB_ACCEPT_C = 'prob_accept_c'
CBH_RESCUE = 'cum_bh'
N_OFFERS_TILL_RESCUE = 'n_offers_till_rescue'
ACTION_CODE = 'action_code'
ACTION_CODE_DESC = 'action_code_desc'
FINAL_REC_URG_AT_TRANSPLANT = 'final_rec_urg_at_transplant'
LAST_NONNT_HU = 'last_nonnt_hu'
MATCH_COMMENT = 'match_comment'
MTCH_TIER = 'mtch_tier'
MTCH_TIER_REV = 'mtch_tier_rev'
OFFERED_TO = 'offered_to'

# Additional information needed for post-txp survival model
YEAR_TXP_RT2014 = 'year_txp_rt2014'


# Recipient information
ID_RECIPIENT = 'id_recipient'
ID_REGISTRATION = 'id_registration'
ID_REREGISTRATION = 'id_reregistration'
TIME_REGISTRATION = 'inc_date_time'
KTH_REG = 'kth_registration'
NTH_TRANSPLANT = 'n_previous_transplants'
LISTING_DATE = 'inc_date_time'
TIME_TO_DEREG = 'time_to_dereg'
TIME_WAITED = 'time_waited'
TIME_SINCE_PREV_TXP = 'time_since_prev_txp'
TIME_SINCE_PREV_TXP_CAT = 'time_to_fail_group'
SURV_OFFSET = 'surv_offset'
PREV_TXP_LIVING = 'prev_txp_living'
PREV_TXP_PED = 'prev_txp_ped'
REREG_RETURN_DIAL_TIME = 'rereg_return_dial_time'
PATIENT_SEX = 'patient_sex'
PATIENT_FEMALE = 'patient_female'
AM_STATUS = 'am_status'

# Retransplantation information
TYPE_RETX = 'type_retx'
NO_RETRANSPLANT = 'retx_none'
RETRANSPLANT = 'retx_real'
IS_RETRANSPLANT = 'retransplant'
DATE_RETRANSPLANT = 'date_retx'
RETRANSPLANT_DURING_SIM = 'retx_simperiod'
PREV_TX_DATE = 'prev_tx_date'
PATIENT_RELISTING = 'rereg'
PATIENT_FAILURE = 'patient_failure'
PATIENT_RELISTING_DATE = 'patient_relist_date'
PATIENT_RELISTING_DATE_RAW = 'patient_relist_date_raw'
PATIENT_FAILURE_DATE = 'patient_failure_date'
PATIENT_FAILURE_DATE_RAW = 'patient_failure_date_raw'
TIME_LIST_TO_REALEVENT = 'time_list_to_realevent'
TIME_TO_REREG = 'time_to_rereg'
EXIT_STATUS_REREG = 'exit_status_rereg'
RETRANSPLANTED = 'retx'
EVENT = 'event'
TIME_TO_EVENT = 'time_to_event'
TIME_TO_PATIENT_FAILURE = 'time_to_patient_failure'
TIME_TO_RETX = 'time_to_retx'
TIME_TO_CENS = 'time_to_cens'
SIM_START_TIME = 'sim_start_time'

# Previously accrued waiting time
PREVIOUS_T = 'previous_t'

# Retransplantation information for post-transplant module
AGE_PREV_TXP = 'age_prev_txp'
OFFER = 'offer'

# Donor information
D_DATE = 'd_date'
ID_DONOR = 'id_donor'
D_CMV = 'd_cmv'
D_HBSAG = 'graft_hbsag'
D_HCVAB = 'graft_hcvab'
D_HBCAB = 'graft_hbcab'
D_SEPSIS = 'graft_sepsis'
D_MENINGITIS = 'graft_meningitis'
D_MALIGNANCY = 'donor_malignancy'
D_DRUG_ABUSE = 'donor_drug_abuse'
D_MARGINAL_FREE_TEXT = 'donor_marginal_free_text'
D_TUMOR_HISTORY = 'donor_tumor_history'
D_DCD = 'graft_dcd'
D_EUTHANASIA = 'graft_euthanasia'
DEATH_CAUSE_GROUP = 'death_cause_group'
D_SEX = 'd_sex'
D_SMOKING = 'donor_smoking'
D_DIABETES = 'donor_diabetes'
D_HYPERTENSION = 'donor_hypertension'
D_LAST_CREAT = 'donor_last_creat'
D_URINE_PROTEIN = 'd_urine_protein'
D_RESCUE = 'd_rescue'
D_CARREST = 'donor_carrest'
RESCUE = 'rescue'
INTERNATIONAL_RESCUE = 'international_rescue'
NONLOCAL_ESP_DONOR = 'nonlocal_esp_donor'
NONLOCAL_DCD_DONOR = 'nonlocal_dcd_donor'

# Status update information
TYPE_UPDATE = 'type_update'
EXIT_DATE = 'exit_date'
EXIT_STATUS = 'exit_status'
TSTART = 'tstart'
TSTOP = 'tstop'
URGENCY_REASON = 'urgency_reason'
REMOVAL_REASON = 'removal_reason'
R_DOB = 'r_dob'
R_AGE_LISTING = 'age_at_listing'
STATUS_VALUE = 'variable_value'
STATUS_DETAIL = 'variable_detail'
STATUS_DETAIL2 = 'variable_detail2'

# Status update types
URG = 'URG'
FU = 'FU'
R = 'R'
D = 'D'
FU = 'FU'
HU = 'HU'
T = 'T'
NT = 'NT'

T1 = 'T1'
T3 = 'T3'
FP = 'FP'
RR = 'RR'
CR = 'CR'

# Match geography
A = 'A'     # abroad
H = 'H'     # home
R = 'R'     # regional
S = 'S'     # subregional
L = 'L'     # local
N = 'N'     # national
INT = 'I'   # international

# Event types
PAT = 'patient'
BAL = 'balance'
SNAPSHOT = 'snapshot'
DON = 'don'

# Indicator columns for blood group rules.
PATIENT_IS_HU = 'patient_is_hu'

# For calculation of scores
INTERCEPT = 'intercept'
TRAFOS = 'trafos'

# Transformations
IDENTITY_TRAFO = 'i'
LOG = 'log'
AGE_DIFFERENCE_FILTER = 'age_difference_filter'
AGE_DIFFERENCE_FILTER_MUTED = 'age_difference_filter_muted'
EXPONENTIAL_BASE2 = 'exponential_base2'
EXPONENTIAL_BASE5 = 'exponential_base5'
EXPONENTIAL_BASE20 = 'exponential_base20'
EXPONENTIAL_BASE10 = 'exponential_base10'
EXPONENTIAL_BASE50 = 'exponential_base50'
EXPONENTIAL_BASE1000 = 'exponential_base1000'
LINEAR_AGE_FILTER_45_75 = 'decreasing_age_filter_45_75'
LINEAR_AGE_FILTER_45_100 = 'increasing_age_filter_45_100'

AGE_HLA_DECAY_P5 = 'age_hla_decay_p5'
AGE_HLA_DECAY_P6 = 'age_hla_decay_p6'
AGE_HLA_DECAY_P7 = 'age_hla_decay_p7'
AGE_HLA_DECAY_P5_R40 = 'age_hla_decay_p5_r40'
AGE_HLA_DECAY_P10 = 'age_hla_decay_p10'

PATIENT_RANK = 'recipient_rank'
RANK = 'rnk_other'

# Blood groups
BG_A = 'A'
BG_AB = 'AB'
BG_O = 'O'
BG_B = 'B'

# Profile variables
PROFILE_COMPATIBLE = 'prof_compatible'
PROFILE_MIN_DONOR_AGE = 'profile_min_donor_age'
PROFILE_MAX_DONOR_AGE = 'profile_max_donor_age'
PROFILE_DCD = 'profile_dcd'
PROFILE_HBSAG = 'profile_hbsag'
PROFILE_HCVAB = 'profile_hcvab'
PROFILE_HBCAB = 'profile_hbcab'
PROFILE_SEPSIS = 'profile_sepsis'
PROFILE_MENINGITIS = 'profile_meningitis'
PROFILE_MALIGNANCY = 'profile_malignancy'
PROFILE_DRUG_ABUSE = 'profile_drug_abuse'
PROFILE_RESCUE = 'profile_rescue'
PROFILE_EUTHANASIA = 'profile_euthanasia'
PROFILE_ESP = 'profile_esp'
PROFILE_EXTENDED_ESP = 'profile_extended_esp'

ALLELE = 'allele'
SPLIT = 'split'
BROAD = 'broad'

TOTAL_MATCH_POINTS = 'total_match_points'
POINTS_HLA = 'wfmr_xc_mism'
POINTS_HLA_PED = 'wfmr_xc_mism_ped'
POINTS_HLA_AGE_MATCH = 'hla_age'
POINTS_AGE = 'points_age'

POINTS_ETKAS = 'points_etkas'
MULTIPLIER_ETKAS = 'multiplier_etkas'
TRAFOS_ETKAS = 'trafos_etkas'
POINTS_ESP = 'points_esp'

# Travel info
FROM_CENTER = 'from_center'
TO_CENTER = 'to_center'
DRIVING_TIME = 'driving_time'

# HLA match potentials
HMP_FH = 'hmp_FH'
HMP_1BDR = 'hmp_1BDR'
HMP_1ABDR = 'hmp_1ABDR'
HMP_0DR = 'hmp_0DR'
