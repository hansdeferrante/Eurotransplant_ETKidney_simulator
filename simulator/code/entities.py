#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

import warnings
from typing import (
    Dict, Tuple, Optional, Callable, Any
)
from datetime import timedelta, datetime
from copy import deepcopy
from warnings import warn

import pandas as pd
import numpy as np
from simulator.code.utils.utils import (
    round_to_int, round_to_decimals, DotDict,
    nanOrNone
)
import simulator.magic_values.column_names as cn

import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
from simulator.code.events.StatusUpdate import StatusUpdate, ProfileUpdate
from simulator.code.events.PatientStatusQueue import PatientStatusQueue
from simulator.code.utils.SimResults import SimResults
from simulator.code.HLA.HLASystem import HLAProfile, Unacceptables

from simulator.magic_values.inputfile_settings import DEFAULT_DATE_TIME_HMS


class Donor:
    """
    Class which implements a donor.

    Attributes   #noqa
    ----------
    id_donor : int
        ID for the donor.
    donor_country : str
        Name of the country of donation.
    donor_region : str
        Name of the donor region.
    donor_center : str
        Name of the donating center (5-letter code).
    donor_age : int
        Donor age.
    bloodgroup : str
        Blood group of the donor.
    reporting_date : datetime
        Time of entry on the kidney waitlist.
    graft_dcd : bool
        Donor after circulatory death status.
    donor_weight : float
        Weight of the donor.
    donor_height : float
        Height of the donor.
    sepsis : bool
        Sepsis status.
    meningitis : bool
        Meningitis status.
    malignancy : bool
        Malignancy status.
    drug_abuse : bool
        Drug abuse status.
    euthanasia : bool
        Euthanasia status.
    rescue : bool
        Rescue status.
    death_cause_group : str
        Grouping of causes of death.
    hla_string : str
        Human leukocyte antigen string.
    hla : Optional['HLAProfile']
        HLA profile object, if available.
    n_kidneys_available : int
        Number of kidneys available from the donor.
    hypertension : bool
        Hypertension status.
    diabetes : int
        Diabetes status.
    cardiac_arrest : bool
        Cardiac arrest status.
    last_creat : float
        Last creatinine level of the donor.
    urine_protein : int
        Level of urine protein.
    smoker : bool
        Smoking status.
    hbsag : bool
        Hepatitis B surface antigen status.
    hcvab : bool
        Hepatitis C antibody status.
    hbcab : bool
        Hepatitis B core antibody status.
    cmv : bool
        Cytomegalovirus status.
    donor_marginal_free_text : Optional[int]
        Free text for marginal donor information.
    tumor_history : Optional[int]
        Tumor history of the donor.
    height : Optional[float]
        Height of the donor.
    weight : Optional[float]
        Weight of the donor.
    _needed_match_info : dict
        Cached needed match info.
    _offer_inherit_cols : dict
        Cached columns inherited for an offer.

    Methods
    -------
    arrival_at(sim_date: datetime):
        Returns the time at which donor arrives, as time from simulation date.
    get_pediatric(ped_fun: Callable) -> bool:
        Checks whether the donor is pediatric.
    """

    def __init__(
        self, sim_set: DotDict, id_donor: int, n_kidneys_available: int,
        bloodgroup: str, donor_country: str,
        donor_center: str, donor_region: str,
        donor_subregion_esp: str,
        reporting_date: datetime, donor_dcd: int,
        hla: str,
        hypertension: bool,
        diabetes: int,
        cardiac_arrest: bool,
        last_creat: float,
        urine_protein: int,
        smoker: bool,
        age: int, hbsag: bool,
        hcvab: bool, hbcab: bool,
        cmv: bool, sepsis: bool, meningitis: bool,
        malignancy: bool, drug_abuse: bool,
        euthanasia: bool,
        rescue: bool, death_cause_group: str,
        hla_system: 'entities.HLASystem' = None,
        donor_marginal_free_text: Optional[int] = 0,
        tumor_history: Optional[int] = 0,
        height: Optional[float] = None,
        weight: Optional[float] = None
    ) -> None:

        self.sim_set = sim_set
        self.id_donor = id_donor
        self.reporting_date = reporting_date
        self.n_kidneys_available = n_kidneys_available

        # Geographic information
        assert donor_country in es.ET_COUNTRIES, \
            f'donor country should be one of:\n\t' \
            f'{", ".join(es.ET_COUNTRIES)}, not {donor_country}'
        self.donor_country = donor_country
        self.donor_region = (
            donor_region if donor_region in sim_set.ALLOWED_REGIONS
            else sim_set.CENTERS_TO_REGIONS[donor_center]
        )
        self.donor_center = donor_center
        self.donor_subregion_esp = (
            donor_subregion_esp if isinstance(donor_subregion_esp, str)
            else donor_center
        )
        self.__dict__[cn.D_PED] = None

        # Donor blood group
        assert bloodgroup in es.ALLOWED_BLOODGROUPS, \
            f'blood group should be one of:\n\t' \
            f'{", ".join(es.ALLOWED_BLOODGROUPS)}'
        self.d_bloodgroup = str(bloodgroup)

        self.graft_dcd = bool(donor_dcd)

        # Profile info
        self.age_esp_eligible = self.sim_set.DONOR_AGE_ESP_ELIGIBLE.get(
            self.donor_country,
            65
        )
        if age is not None:
            self.__dict__[cn.D_AGE] = age
            self.__dict__[cn.ESP_DONOR] = str(
                int(age >= self.age_esp_eligible)
            )
            if age < 18:
                self.__dict__[cn.D_AGE_GROUP] = 'lt18'
                self.__dict__[cn.DONOR_BALANCE_AGE_GROUP] = mgr.LT18
            elif age < 35:
                self.__dict__[cn.D_AGE_GROUP] = '18_to_34'
                self.__dict__[cn.DONOR_BALANCE_AGE_GROUP] = mgr.YOUNGADULT
            elif age < 50:
                if age < 45:
                    self.__dict__[cn.D_AGE_GROUP] = '35_to_44'
                else:
                    self.__dict__[cn.D_AGE_GROUP] = '45_to_54'
                self.__dict__[cn.DONOR_BALANCE_AGE_GROUP] = mgr.YOUNGADULT
            elif age < 65:
                if age < 55:
                    self.__dict__[cn.D_AGE_GROUP] = '45_to_54'
                else:
                    self.__dict__[cn.D_AGE_GROUP] = '55_to_64'
                self.__dict__[cn.DONOR_BALANCE_AGE_GROUP] = mgr.OLDADULT
            else:
                self.__dict__[cn.D_AGE_GROUP] = '65p'
                self.__dict__[cn.DONOR_BALANCE_AGE_GROUP] = mgr.ELDERLY
        self.__dict__[cn.D_HBSAG] = hbsag if hbsag is not None and  \
            hbsag is not np.nan else False
        self.__dict__[cn.D_HCVAB] = hcvab if hcvab is not None and \
            hcvab is not np.nan else False
        self.__dict__[cn.D_HBCAB] = hbcab if hbcab is not None and \
            hbcab is not np.nan else False
        self.__dict__[cn.D_CMV] = cmv if cmv is not None and  \
            cmv is not np.nan else False
        self.sepsis = sepsis if sepsis is not None and \
            sepsis is not np.nan else False
        self.meningitis = meningitis if meningitis is not None and \
            meningitis is not np.nan else False
        self.malignancy = malignancy if malignancy is not None and \
            malignancy is not np.nan else False
        self.drug_abuse = drug_abuse if drug_abuse is not None and \
            drug_abuse is not np.nan else False
        self.euthanasia = euthanasia if euthanasia is not None and \
            euthanasia is not np.nan else False
        self.rescue = rescue if rescue is not None and \
            rescue is not np.nan else False
        self.death_cause_group = death_cause_group \
            if death_cause_group is not None and \
            death_cause_group is not np.nan else False

        self.donor_weight = weight
        self.donor_height = height

        # Extra information needed for predicting acceptance
        self.__dict__[cn.D_MARGINAL_FREE_TEXT] = donor_marginal_free_text \
            if not nanOrNone(donor_marginal_free_text) else 0
        self.__dict__[cn.D_TUMOR_HISTORY] = tumor_history \
            if not nanOrNone(tumor_history) else 0
        self.__dict__[cn.D_DIABETES] = diabetes \
            if not nanOrNone(diabetes) else 0
        self.__dict__[cn.D_SMOKING] = smoker \
            if not nanOrNone(smoker) else 0
        self.__dict__[cn.D_HYPERTENSION] = hypertension \
            if not nanOrNone(hypertension) else 0
        self.__dict__[cn.D_LAST_CREAT] = last_creat \
            if not nanOrNone(last_creat) else 0
        self.__dict__[cn.D_URINE_PROTEIN] = int(urine_protein) \
            if not nanOrNone(urine_protein) else 0
        self.__dict__[cn.D_CARREST] = cardiac_arrest \
            if not nanOrNone(cardiac_arrest) else 0

        self._needed_match_info = None
        self._offer_inherit_cols = None

        self.hla_string = hla
        if hla_system:
            self.hla = HLAProfile(hla, hla_system)
            self.all_antigens = frozenset(hla_system.return_all_antigens(hla))

    def arrival_at(self, sim_date: datetime) -> float:
        """Retrieve calendar time arrival time."""
        return (
            self.reporting_date - sim_date
        ) / timedelta(days=1)

    def get_pediatric(self, ped_fun: Callable) -> bool:
        """Check whether donor is pediatric"""
        if self.__dict__[cn.D_PED] is None:
            self.__dict__[cn.D_PED] = ped_fun(self)
        return self.__dict__[cn.D_PED]

    def __str__(self) -> str:
        return (
            f'Donor {self.id_donor}, reported on '
            f'{self.reporting_date} in {self.donor_center}'
            f', dcd: {self.graft_dcd}, bg: {self.d_bloodgroup} '
            f'and age {self.__dict__[cn.D_AGE]} with '
            f'{self.n_kidneys_available} ESP/ETKAS kidneys'
        )

    def __repr__(self) -> str:
        return (
            f'Donor {self.id_donor}, reported on '
            f'{self.reporting_date} in {self.donor_center}'
            f', dcd: {self.graft_dcd}, bg: {self.d_bloodgroup}'
            f' and age {self.__dict__[cn.D_AGE]}'
        )

    @classmethod
    def from_dummy_donor(cls, **kwargs):
        # Provide default arguments
        default_args = {
            'id_donor': 1255,
            'donor_country': 'Belgium',
            'donor_region': 'BLGTP',
            'donor_center': 'BLGTP',
            'donor_subregion_esp': 'Bel_1',
            'bloodgroup': 'O',
            'reporting_date': pd.Timestamp('2000-01-01'),
            'weight': 50,
            'donor_dcd': False,
            'hla': None,
            'hla_system': None,
            'hypertension': False,
            'diabetes': False,
            'cardiac_arrest': False,
            'last_creat': 1.5,
            'smoker': False,
            'age': 45,
            'hbsag': False,
            'hcvab': False,
            'hbcab': False,
            'sepsis': False,
            'meningitis': False,
            'malignancy': False,
            'drug_abuse': False,
            'euthanasia': False,
            'rescue': False,
            'death_cause_group': 'Anoxia',
            'n_kidneys_available': 2,
            'urine_protein': 1,
            'cmv': False
        }
        # Update default arguments with provided kwargs
        default_args.update(kwargs)
        return cls(**default_args)

    @property
    def needed_match_info(self) -> Dict[str, Any]:
        if not self._needed_match_info:
            self._needed_match_info = {
                k: self.__dict__[k]
                for k
                in es.MTCH_RCRD_DONOR_COLS
            }
        return self._needed_match_info

    @property
    def offer_inherit_cols(self) -> Dict[str, Any]:
        if not self._offer_inherit_cols:
            self._offer_inherit_cols = {
                k: self.__dict__[k] for k in es.OFFER_INHERIT_COLS['donor']
            }
        return self._offer_inherit_cols

    def __deepcopy__(self, memo):
        # Create a new instance of the class
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        # Copy all attributes, except for hla_system and bal_system
        for k, v in self.__dict__.items():
            if k in {'hla', 'all_antigens', 'sim_set'}:
                setattr(result, k, v)  # Shallow copy of HLA
            else:
                setattr(result, k, deepcopy(v, memo))  # Deep copy

        return result


class Patient:
    """Class which implements a patient

    Attributes   #noqa
    ----------
    id_recipient: int
        integer for the recipient ID
    r_dob: datetime
        patient DOB
    age_at_listing: int
        age in days at listing
    recipient_country: str
        name of country of listing (not nationality!)
    recipient_region: str
        name of region where candidate's from
    recipient_center: str
        name of listing center (5-letter code).
    bloodgroup: str
        recipient blood group
    hla_string: str
        string of the HLA
    hla: HLAProfile
        instance of class HLAProfile
    listing_date: datetime
        Time of entry on the kidney waitlist
    urgency_code: str
        Status of listing (HU/NT/T/HI/I)
    urgency_reason: Optional[str]
        urgency reason
    hu_since: float
        time since registration since when patient is HU
    future_statuses: Optional[PatientStatusQueue]
        status queue (ordered list of status updates)
    exit_status: Optional[str]
        exit status
    exit_date: Optional[datetime]
        exit time
    sim_set: DotDict,
        simulation settings under which pat. is created
    profile: Optional[Profile]
        patient's donor profile
    patient_sex: str
        sex of the patient,
    time_since_prev_txp: float
        Time elapsed since previous TXP
    id_registration: int
        registration ID [optional]
    seed: int
        Optional seed for random number generator
    time_to_dereg: float
        time to deregistration (in days from registration)
    rng_acceptance: Generator
        Random number generator for acceptance probabilities

    Methods
    -------
    get_age_at_t(time: float)
        Return age at time t
    get_next_update_time()
        Return time at which next update is issued
    get_acceptance_prob():
        Return an acceptance probability
    """

    def __init__(
        self, id_recipient: int, recipient_country: str,
        recipient_center: str, recipient_region: str, bloodgroup: str,
        listing_date: datetime, urgency_code: str,
        hla: str,
        r_dob: datetime,
        sim_set: DotDict,
        sex: Optional[str] = None,
        type_retx: str = cn.NO_RETRANSPLANT,
        time_since_prev_txp: Optional[float] = None,
        prev_txp_ped: Optional[bool] = None,
        prev_txp_living: Optional[bool] = None,
        id_reg: Optional[int] = None,
        date_first_dial: Optional[datetime] = None,
        seed: Optional[int] = None,
        hla_system: Optional['entities.HLASystem'] = None,
        mmp_system: Optional['entities.MMPSystem'] = None,
        previous_wt: Optional[float] = None,
        kidney_program: Optional[str] = None
    ) -> None:

        self.id_recipient = int(id_recipient)
        self.id_registration = id_reg if id_reg else None
        self.sim_set = sim_set
        self._seed_acceptance_behavior(seed=seed)

        # Load patient bloodgroup
        assert bloodgroup in es.ALLOWED_BLOODGROUPS, \
            f'blood group should be one of:\n\t' \
            f'{", ".join(es.ALLOWED_BLOODGROUPS)}'
        self.r_bloodgroup = str(bloodgroup)

        # Check listing status
        assert urgency_code in es.ALLOWED_STATUSES, \
            f'listing status should be one of:\n\t' \
            f'{", ".join(es.ALLOWED_STATUSES)}, not {urgency_code}'
        self.urgency_code = str(urgency_code)
        self.__dict__[cn.PATIENT_IS_HU] = str(urgency_code) == cn.HU
        self.__dict__[cn.ANY_HU] = self.__dict__[cn.PATIENT_IS_HU]
        self.__dict__[cn.AM_STATUS] = False
        self.__dict__[cn.ANY_ACTIVE] = False
        self.urgency_reason = None

        # Set other information relating to listing statuses
        # Status updates
        self.future_statuses = None
        self.__dict__[cn.EXIT_STATUS] = None
        self.__dict__[cn.EXIT_DATE] = None

        # Set patient center, region, and country
        self.__dict__[cn.RECIPIENT_CENTER] = str(recipient_center)
        self.__dict__[cn.RECIPIENT_REGION] = (
            str(recipient_region)
            if recipient_region in sim_set.ALLOWED_REGIONS
            else sim_set.CENTERS_TO_REGIONS[recipient_center]
        )
        self.__dict__[cn.RECIPIENT_COUNTRY] = str(recipient_country)
        self.__dict__[cn.PATIENT_COUNTRY] = str(recipient_country)

        # Set date of birth, listing date, and dial date
        if not isinstance(r_dob, datetime):
            raise TypeError(
                f'r_dob must be datetime, not a {type(r_dob)}'
            )
        self.r_dob = r_dob
        self.initialize_etkas_esp_eligibility()
        self.set_listing_date(listing_date=listing_date)
        self.set_dial_date(date_first_dial)

        # Initialize patient information
        self.initialized = False
        self.profile = None
        self.__dict__[cn.PATIENT_SEX] = sex
        self.__dict__[cn.PATIENT_FEMALE] = int(sex == 'Female')
        self.disease_group = None

        # Set transplant history & previous waiting time
        self.__dict__[cn.TIME_SINCE_PREV_TXP] = (
            None if pd.isnull(time_since_prev_txp) else time_since_prev_txp
        )
        if time_since_prev_txp is None or pd.isnull(time_since_prev_txp):
            self.prev_txp_ped = 0
            self.prev_txp_living = 0
            self.__dict__[cn.AGE_PREV_TXP] = None
            self.__dict__[cn.TIME_SINCE_PREV_TXP]
        else:
            self.prev_txp_ped = prev_txp_ped
            self.prev_txp_living = prev_txp_living
            self.__dict__[cn.TIME_SINCE_PREV_TXP] = time_since_prev_txp
            self.__dict__[cn.AGE_PREV_TXP] = (
                (
                    self.age_days_at_listing
                ) - time_since_prev_txp
            ) / es.DAYS_PER_YEAR_FLOAT
        assert type_retx in (
            cn.NO_RETRANSPLANT, cn.RETRANSPLANT, cn.RETRANSPLANT_DURING_SIM
        )
        self.__dict__[cn.TYPE_RETX] = type_retx
        self.__dict__[cn.IS_RETRANSPLANT] = (
            0 if type_retx and type_retx == cn.NO_RETRANSPLANT
            else 1 if type_retx
            else None
        )

        # Add dummy indicating whether candidate is eligible
        # for returned dial time
        self.__dict__[cn.REREG_RETURN_DIAL_TIME] = (
            (
                self.__dict__[cn.TIME_SINCE_PREV_TXP] <=
                es.CUTOFF_REREG_RETURN_DIAL_TIME
            ) if self.__dict__[cn.TIME_SINCE_PREV_TXP]
            else False
        )
        self.pediatric_retransplant = (
            self.__dict__[cn.IS_RETRANSPLANT] and
            self.__dict__[cn.AGE_PREV_TXP] and
            self.__dict__[cn.AGE_PREV_TXP] < 18
        ) or (
            self.prev_txp_ped and self.__dict__[cn.TIME_SINCE_PREV_TXP] <= 91
        )
        if previous_wt:
            self.__dict__[cn.PREVIOUS_WT] = previous_wt / \
                es.DAYS_PER_YEAR_FLOAT
        else:
            self.__dict__[cn.PREVIOUS_WT] = 0

        # Geographic information
        assert recipient_country in es.ET_COUNTRIES, \
            f'listing country should be one of:\n\t' \
            f'{", ".join(es.ET_COUNTRIES)}, not {recipient_country}'

        # Initialize HLA information
        self.hla_system = hla_system
        self.mmp_system = mmp_system
        for fb in es.HLA_FORBIDDEN_CHARACTERS:
            if type(hla) is str and fb in hla:
                warnings.warn(
                    f'Forbidden character {fb} in patient '
                    f'HLA-input string!\n({hla})'
                )
                break

        if type(hla) is str:
            self.set_match_hlas(hla_string=hla)
        else:
            self.set_match_hlas(hla_string='')

        # Participation in ETKAS or ESP (needed for Germany, if ESP)
        self.set_chosen_program(kidney_program)

        # Initialize empty unacceptables
        self.unacceptable_antigens = Unacceptables(hla_system)

        # Items we want to store, for rapid access through functions
        # and properties
        self._vpra = None
        self.valid_pra = 0
        self.__dict__[cn.R_PED] = None
        self._needed_match_info = None
        self._offer_inherit_cols = None
        self._active = None
        self._pediatric = None
        self._time_dial_start = None
        self._et_mmp = None
        self._am = None
        self._bg = None
        self._dcd_country = None
        self._eligibility_last_assessed = -1e10
        self._time_dial_at_listing = None

    def set_diag(self, upd: StatusUpdate) -> None:
        if self.__dict__[cn.DISEASE_GROUP] != upd.status_detail:
            self.__dict__[cn.DISEASE_GROUP] = upd.status_detail
            self.__dict__[cn.DISEASE_SINCE] = upd.arrival_time

    def set_dial_date(self, date_first_dial: Optional[datetime]) -> None:
        # Set date first dialysis. Also remove the cached pediatric status
        # of the candidate; the candidate could previously become pediatric
        # if dialysis starts at age 16, despite age at listing of 16.
        self.date_first_dial = date_first_dial
        self._time_dial_start = None
        if date_first_dial:
            self.age_first_dial = (
                (date_first_dial - self.r_dob) / es.DAYS_PER_YEAR
            )
            self._pediatric = None
        else:
            self.age_first_dial = None
            self._pediatric = None

    def set_chosen_program(self, kidney_program: Optional[str]) -> None:
        # Process an update to ETKAS / ESP choice for candidate.
        self.__dict__[cn.KIDNEY_PROGRAM] = (
            kidney_program if isinstance(kidney_program, str) else mgr.ETKAS
        )
        self._esp_eligible = None
        self._etkas_eligible = None
        self._eligibility_last_assessed = -1e10

    def set_listing_date(self, listing_date: datetime) -> None:
        assert isinstance(listing_date, datetime), \
            f'listing date should be datetime, not {type(listing_date)}'

        self.__dict__[cn.LISTING_DATE] = listing_date
        self.age_days_at_listing = round_to_int(
            (listing_date - self.r_dob) / timedelta(days=1)
        )

        self.listing_offset = (
            listing_date - self.sim_set.SIM_START_DATE
        ) / timedelta(days=1)

    def initialize_etkas_esp_eligibility(self) -> None:
        # Calculate date at which patient becomes ESP eligible,
        # as well as simulation time when patient becomes ESP
        # eligible
        self.age_esp_eligible = self.sim_set.CAND_AGE_ESP_ELIGIBLE.get(
            self.__dict__[cn.RECIPIENT_COUNTRY]
        )
        self.date_esp_eligible = (
            self.r_dob + es.DAYS_PER_YEAR *
            self.age_esp_eligible
        )
        self.sim_time_esp_eligible: float = (
            self.date_esp_eligible - self.sim_set.SIM_START_DATE
        ) / timedelta(days=1)
        self.always_etkas_aged: bool = (
            self.date_esp_eligible > self.sim_set.SIM_END_DATE
        )
        self.country_always_etkas_eligible: bool = (
            not self.__dict__[cn.RECIPIENT_COUNTRY]
            in es.COUNTRIES_ESP_ETKAS_MUTUALLY_EXCLUSIVE
        )
        self.program_choice_required_after: float = (
            self.sim_set['times_esp_eligible'].get(
                self.__dict__[cn.RECIPIENT_COUNTRY],
                1e10
            )
        )

    def set_match_hlas(self, hla_string: str) -> None:
        self.hla_string = hla_string
        if self.hla_system:
            self.hla = HLAProfile(hla_string, hla_system=self.hla_system)
            self.__dict__.update(
                self.hla.homozygosity_per_locus
            )

    def set_urgency_code(
        self, code: str, reason: Optional[str] = None
    ) -> None:
        """Update urgency code with string"""
        if code not in es.ALL_STATUS_CODES:
            raise Exception(
                f'listing status should be one of:\n\t'
                f'{", ".join(es.ALL_STATUS_CODES)},\n not {code}'
            )
        self.__dict__[cn.URGENCY_CODE] = code
        self.__dict__[cn.PATIENT_IS_HU] = code == cn.HU
        self._active = None
        if reason:
            self.__dict__[cn.URGENCY_REASON] = reason

    def set_transplanted(
            self, tx_date: datetime,
            donor: Optional[Donor] = None,
            sim_results: Optional[SimResults] = None,
            match_record: Optional[
                'simulator.code.matchlist.MatchList.MatchRecord'
            ] = None
    ) -> None:
        """Set patient to transplanted"""
        if self.future_statuses:
            self.future_statuses.clear()
        self.__dict__[cn.EXIT_DATE] = tx_date
        self.__dict__[cn.EXIT_STATUS] = cn.FU
        self.__dict__[cn.FINAL_REC_URG_AT_TRANSPLANT] = (
            self.__dict__[cn.URGENCY_CODE]
        )
        self.set_urgency_code(cn.FU)
        if sim_results:
            sim_results.save_transplantation(
                pat=self, matchr=match_record
            )

    def get_acceptance_prob(self) -> float:
        """Simulate an acceptance probability"""
        prob = self.rng_acceptance.random()
        self.__dict__[cn.DRAWN_PROB] = round_to_decimals(prob, 3)
        return prob

    def get_dial_time_sim_start(self) -> Optional[float]:
        """Get dialysis start time, relative to simulation time"""
        if self._time_dial_start:
            return self._time_dial_start
        if self.date_first_dial:
            self._time_dial_start = (
                self.date_first_dial - self.sim_set.SIM_START_DATE
            ) / timedelta(days=1)
            return self._time_dial_start
        else:
            return None

    def get_dial_time_at_listing(self) -> Optional[float]:
        """Get dialysis start time, relative to simulation time"""
        if self._time_dial_at_listing:
            return self._time_dial_at_listing
        if self.date_first_dial is not None:
            self._time_dial_at_listing = (
                self.date_first_dial - self.__dict__[cn.TIME_REGISTRATION]
            ) / timedelta(days=1)
            return self._time_dial_at_listing
        else:
            return 0

    def get_pediatric(self, ped_fun: Callable, match_age: float) -> bool:
        """Check whether recipient is pediatric. Only check for
            candidates who were pediatric at past matches.
        """
        if self._pediatric is True or self._pediatric is None:
            self._pediatric = ped_fun(
                candidate_country=self.__dict__[cn.PATIENT_COUNTRY],
                age_first_dial=self.age_first_dial,
                age_at_listing=(
                    self.age_days_at_listing / es.DAYS_PER_YEAR_FLOAT
                ),
                match_age=match_age,
                prev_txp_ped=self.pediatric_retransplant,
                time_since_prev_txp=self.__dict__[cn.TIME_SINCE_PREV_TXP]
            )
            self.__dict__[cn.R_PED] = self._pediatric
            return self._pediatric
        if self._pediatric is False:
            return self._pediatric

    def get_et_mmp(self) -> float:
        """Calculates the ET mismatch probability, under the assumption that
            HLAs, blood types, and unacceptables are independent. This is how
            ET calculates the mismatch probability.
        """
        if self._et_mmp:
            return self._et_mmp
        elif self.hla_system and self.mmp_system:
            et_don_mmprobs, et_hla_mmfreqs = (
                self.mmp_system.calculate_broad_split_mmps(
                    p=self
                )
            )
            self._et_mmp = et_don_mmprobs.get(
                cn.ET_MMP_SPLIT,
                et_don_mmprobs.get(cn.ET_MMP_BROAD, 0)
            )
            self.__dict__[cn.ET_HLA_MISMATCHFREQ] = et_hla_mmfreqs.get(
                cn.ET_MMP_SPLIT, et_hla_mmfreqs.get(cn.ET_MMP_BROAD, 0)
            )

            return self._et_mmp
        else:
            raise Exception(
                'Cannot calculate ET-MMP, '
                'without HLA / MMP system.')

    def get_hla_mismatchfreqs(self) -> Dict[str, float]:
        """Get the HLA mismatch frequency, i.e. the probability
            no suitable donor among next 1,000 based on counting
            (and only based on HLA). For this, we first calculate
            the empirical hla match probability, and then calculate
            the hla match frequency among 1,000 donors based on this.

        """
        if self.__dict__.get(cn.HLA_MISMATCHFREQ):
            return self.__dict__[cn.HLA_MISMATCHFREQ]
        elif self.hla_system and self.mmp_system:
            if self.hla_string not in self.mmp_system.match_potentials:
                print(
                    f'Warning: re-calculating match potential '
                    f'for:\n\t{self.hla_string}\n of patient '
                    f'{self.id_recipient}.\n These were not initialized.\n'
                )
                hmps = self.mmp_system.calculate_hla_match_probability(
                    self.hla,
                    hla_match_pot_definitions=(
                        es.HLA_FAVORABLE_MATCH_DEFINITIONS
                    )
                )
                hla_matchfreqs = {
                    k.lower().replace(
                        'hmp_',
                        'hlamismatchfreq_'
                    ): round(100 * (1 - v)**1000, 5)
                    for k, v in hmps.items()
                }
                self.mmp_system.match_potentials[self.hla_string] = {
                    **hmps,
                    **hla_matchfreqs
                }
            if self.mmp_system.match_potentials[self.hla_string]:

                self.__dict__[cn.HLA_MISMATCHFREQ] = (
                    self.mmp_system.match_potentials[self.hla_string]
                )
                return self.__dict__[cn.HLA_MISMATCHFREQ]
            else:
                raise Exception(
                    'No match potentials available for {self.hla_string}'
                )
        else:
            raise Exception(
                'Cannot calculate mismatch frequency, without an MMP system.'
            )

    def update_etkas_esp_eligibility(self, s: float) -> None:
        """Update a candidate's ETKAS/ESP eligibility.

        1. Candidates aged less than 65 for the entire
           simulation are ETKAS eligible
        2. Candidates in countries who do not enforce
            ETKAS/ESP choice are ETKAS eligible, and ESP eligible
            if aged 65+ or extended ESP accepting
        3. In countries which enforce ETKAS/ESP choice, <65 aged candidates are
           always ETKAS-eligible, 65+ candidates are ETKAS-eligible if
            they chose ETKAS.
        """
        if self.always_etkas_aged:
            self._etkas_eligible = True
            if self.profile:
                self._esp_eligible = self.profile.extended_esp
            else:
                self._esp_eligible = False
        elif (
            self.country_always_etkas_eligible or
            s <= self.program_choice_required_after
        ):
            self._etkas_eligible = True
        else:
            self._etkas_eligible = (
                s < self.sim_time_esp_eligible or
                (
                    s >= self.sim_time_esp_eligible and
                    self.__dict__[cn.KIDNEY_PROGRAM] == mgr.ETKAS
                )
            )
        self._esp_eligible = (
            s >= self.sim_time_esp_eligible or
            (self.profile and self.profile.extended_esp)
        )
        self._eligibility_last_assessed = s

    def get_etkas_eligible(self, s: float) -> bool:
        """Return whether the patient is ETKAS eligible."""
        if self.always_etkas_aged or self.country_always_etkas_eligible:
            return True
        elif (
                s and
                (
                    s - self._eligibility_last_assessed
                ) >= es.CHECK_ETKAS_ESP_ELIGIBILITY
        ):
            self.update_etkas_esp_eligibility(s=s)
        return self._etkas_eligible

    def get_esp_eligible(self, s: float) -> bool:
        """Return whether the patient is ESP eligible"""
        if (
                s and
                (
                    s - self._eligibility_last_assessed
                ) >= es.CHECK_ETKAS_ESP_ELIGIBILITY
        ):
            self.update_etkas_esp_eligibility(s=s)
        return self._esp_eligible

    def get_next_update_time(self) -> Optional[float]:
        """
            Return time at which next update is issued,
            offset from calendar time
        """
        if (
            self.future_statuses is not None and
            not self.future_statuses.is_empty()
        ):
            return (
                self.listing_offset +
                self.future_statuses.first().arrival_time
            )

    def _seed_acceptance_behavior(self, seed: Optional[int]) -> None:
        """Use common random numbers to pre-set patient acceptance behavior"""
        if seed is None:
            self.seed = 1
        elif self.id_registration is None:
            warn(
                f'No registration ID is set. '
                f'Setting seed to the registration ID.'
            )
            self.seed = self.seed
        else:
            self.seed = seed * self.id_registration
        self.rng_acceptance = np.random.default_rng(self.seed)

    def reset_matchrecord_info(self) -> None:
        self._needed_match_info = None
        self._offer_inherit_cols = None

    def preload_status_updates(
        self,
        fut_stat: PatientStatusQueue
    ) -> None:
        """Initialize status updates for patient"""
        self.initialized = True
        self.__dict__[cn.TIME_LIST_TO_REALEVENT] = (
            fut_stat.return_time_to_exit(
                exit_statuses=es.EXITING_STATUSES
            )
        )
        self.future_statuses = fut_stat

    def do_patient_update(
        self,
        sim_results: Optional[SimResults] = None
    ):
        """Update patient with first next status"""
        if (
            self.future_statuses is not None and
            not self.future_statuses.is_empty()
        ):
            stat_update = self.future_statuses.next()
            upd_time = stat_update.arrival_time

            if not self.__dict__[cn.ANY_ACTIVE]:
                if not stat_update.before_sim_start:
                    if self.__dict__[cn.URGENCY_CODE] != cn.NT:
                        self.__dict__[cn.ANY_ACTIVE] = True

            # Now reset cached match record information;
            # we no longer need them.
            # Also reset other information, inherited to match records
            self.reset_matchrecord_info()

            # Update the patient with the update status
            if self.sim_set.USE_REAL_FU and stat_update.synthetic:
                pass
            elif stat_update.type_status == mgr.URG:
                if stat_update.status_value == cn.FU:
                    if (
                        self.sim_set.USE_REAL_FU or
                        stat_update.before_sim_start
                    ):
                        self.set_transplanted(
                            self.__dict__[cn.LISTING_DATE] + timedelta(
                                days=stat_update.arrival_time
                            )
                        )
                    elif self.am:
                        self.set_transplanted(
                            self.__dict__[cn.LISTING_DATE] + timedelta(
                                days=stat_update.arrival_time
                            )
                        )
                elif stat_update.status_value in es.TERMINAL_STATUSES:
                    self.future_statuses.clear()
                    self.__dict__[cn.FINAL_REC_URG_AT_TRANSPLANT] = (
                        self.urgency_code
                    )
                    self.set_urgency_code(
                        stat_update.status_value,
                        stat_update.status_detail
                    )
                    self.__dict__[cn.EXIT_STATUS] = stat_update.status_value
                    self.__dict__[cn.EXIT_DATE] = (
                        self.__dict__[cn.LISTING_DATE] +
                        timedelta(days=stat_update.arrival_time)
                    )
                    if sim_results:
                        sim_results.save_exit(self)
                else:
                    self._update_urgency(stat_update)
            elif stat_update.type_status == mgr.DIAG:
                self._update_diag(stat_update)
            elif stat_update.type_status == mgr.HLA:
                self._update_hla(stat_update)
            elif stat_update.type_status == mgr.AM:
                self._update_am(stat_update)
            elif stat_update.type_status == mgr.DIAL:
                self._update_dialysis_status(stat_update)
            elif stat_update.type_status == mgr.UNACC:
                self._update_unacceptables(stat_update)
            elif stat_update.type_status == mgr.PRA:
                self._update_pra(stat_update)
            elif stat_update.type_status == mgr.PRG:
                self._update_program(stat_update)
            elif isinstance(stat_update, ProfileUpdate):
                self.profile = stat_update.profile
                self._other_profile_compatible = None
                self._etkas_eligible = None
                self._esp_eligible = None
            else:
                print(
                    f'Stopped early for recipient {self.id_recipient} '
                    f'due to unknown status ({stat_update.type_status}):'
                    f'\n{stat_update.type_status}'
                )
                exit()

    def trigger_historic_updates(
        self
    ) -> None:
        """Trigger all updates which occured before sim start"""
        if self.future_statuses is not None:
            while (
                not self.future_statuses.is_empty() and
                (
                    self.future_statuses.first().before_sim_start or
                    (self.future_statuses.first().arrival_time < 0)
                )
            ):
                self.do_patient_update()
            if self.__dict__[cn.URGENCY_CODE] != 'NT':
                self.__dict__[cn.ANY_ACTIVE] = True
            self.__dict__[cn.INIT_URG_CODE] = (
                self.__dict__[cn.URGENCY_CODE]
            )

    def schedule_death(self, fail_date: datetime) -> None:
        """
            Schedule a death update at the fail date. Needed for constructing
            synthetic re-registrations.
        """
        death_event = StatusUpdate(
            type_status=mgr.URG,
            arrival_time=(
                fail_date - self.__dict__[cn.LISTING_DATE]
            ) / timedelta(days=1),
            status_detail='',
            status_value=mgr.D,
            sim_start_time=(
                self.__dict__[cn.LISTING_DATE] - self.sim_set.SIM_START_DATE
            ) / timedelta(days=1)
        )
        self.future_statuses.add(
            death_event
        )
        self.future_statuses.truncate_after(
            truncate_time=(
                fail_date - self.__dict__[cn.LISTING_DATE]
            ) / timedelta(days=1)
        )

    def _update_am(self, upd: StatusUpdate) -> None:
        self._etkas_eligible = None
        self._esp_eligible = None
        if upd.status_value == 'A':
            self.__dict__[cn.AM_STATUS] = True
            self._am = True
        else:
            self.__dict__[cn.AM_STATUS] = False
            self._am = False

    def _update_program(self, upd: StatusUpdate) -> None:
        self._etkas_eligible = None
        self._esp_eligible = None
        if upd.status_value in {mgr.ETKAS, mgr.ESP}:
            self.set_chosen_program(upd.status_value)
        else:
            raise Exception(
                f'{upd.status_value} is not a valid program choice'
            )

    def _update_diag(self, upd: StatusUpdate) -> None:
        """Update downmarked status"""
        self.set_diag(upd)

    def _update_dialysis_status(self, upd: StatusUpdate) -> None:
        dt = datetime.strptime(upd.status_value, DEFAULT_DATE_TIME_HMS)
        self.set_dial_date(dt)

    def _update_hla(self, upd: StatusUpdate) -> None:
        self.set_match_hlas(upd.status_value)
        self._et_mmp = None

    def _update_unacceptables(self, upd: StatusUpdate) -> None:
        if type(upd.status_value) is str:
            self.unacceptable_antigens = Unacceptables(
                self.hla_system,
                unacc_string=upd.status_value
            )
            self._vpra = float(upd.status_detail)
        else:
            self.unacceptable_antigens = Unacceptables(self.hla_system)
            self._vpra = None
        self._et_mmp = None

    def _update_urgency(self, upd: StatusUpdate) -> None:
        """Update urgency code"""
        assert upd.status_value in es.ALLOWED_STATUSES, \
            f'listing status should be one of:\n\t' \
            f'{", ".join(es.ALLOWED_STATUSES)}, not {upd.status_value}'
        self.set_urgency_code(
            upd.status_value,
            upd.status_detail
        )

        # Track other information
        if self.urgency_code == mgr.T:
            self.__dict__[cn.LAST_NONNT_HU] = False
            if not self.__dict__[cn.ANY_ACTIVE]:
                self.__dict__[cn.ANY_ACTIVE] = True
        elif self.urgency_code == cn.HU:
            self.__dict__[cn.LAST_NONNT_HU] = True
            if not self.__dict__[cn.ANY_ACTIVE]:
                self.__dict__[cn.ANY_ACTIVE] = True
            self.__dict__[cn.ANY_HU] = True

    def _update_pra(self, upd: StatusUpdate):
        self.valid_pra = float(upd.status_value)
        if isinstance(upd.status_detail, str):
            self.__dict__[cn.PRA] = float(upd.status_detail)
        else:
            self.__dict__[cn.PRA] = 0

    @property
    def am(self):
        if self._am is None:
            self._am = self.__dict__[cn.AM_STATUS]
        return self._am

    @property
    def vpra(self) -> float:
        if self._vpra is None:
            self._vpra = self.hla_system.calculate_vpra_from_string(
                self.unacceptable_antigens.unacc_string
            )
        return self._vpra

    @property
    def dcd_country(self) -> bool:
        if self._dcd_country is None:
            self._dcd_country = (
                self.__dict__[cn.RECIPIENT_COUNTRY]
                in es.DCD_ACCEPTING_COUNTRIES
            )
        return self._dcd_country

    @property
    def bloodgroup(self) -> str:
        if self._bg is None:
            self._bg = self.__dict__[cn.R_BLOODGROUP]
        return self._bg

    @property
    def needed_match_info(self) -> Dict[str, Any]:
        if not self._needed_match_info:
            self._needed_match_info = {
                k: self.__dict__[k]
                for k
                in es.MTCH_RCRD_PAT_COLS
            }
        return self._needed_match_info

    @property
    def offer_inherit_cols(self) -> Dict[str, Any]:
        if not self._offer_inherit_cols:
            self._offer_inherit_cols = {
                k: self.__dict__[k] for k in es.OFFER_INHERIT_COLS['patient']
            }
        return self._offer_inherit_cols

    @property
    def active(self) -> bool:
        if self._active is None:
            self._active = (
                self.urgency_code in es.ACTIVE_STATUSES
            )
        return self._active

    def is_initialized(self) -> bool:
        """Whether future statuses have been loaded for the patient"""
        return self.initialized

    def __deepcopy__(self, memo) -> None:
        # Create a new instance of the class
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        # Copy all attributes, except for hla_system and bal_system
        for k, v in self.__dict__.items():
            if k in {
                'hla_system', 'bal_system',
                'center_travel_times', 'calc_points',
                'mmp_system', 'distance_cache',
                'hla', 'unacceptable_antigens',
                'sim_set'
            }:
                setattr(result, k, v)  # Shallow copy
            else:
                setattr(result, k, deepcopy(v, memo))  # Deep copy

        return result

    def __str__(self) -> str:
        if self.__dict__[cn.EXIT_STATUS] is not None:
            return (
                f'Patient {self.id_recipient}, exited on '
                f'{self.__dict__[cn.EXIT_DATE]} with status '
                f' {self.__dict__[cn.EXIT_STATUS]} '
                f'(delisting reason: {self.urgency_reason})'
                f' with bloodgroup {self.r_bloodgroup}'
            )

        return (
            f'Patient {self.id_recipient}, listed on '
            f'{self.__dict__[cn.LISTING_DATE].strftime("%d-%m-%Y")} '
            f'at center {self.__dict__[cn.RECIPIENT_CENTER]} with '
            f'current status {self.urgency_code}'
            f' with bloodgroup {self.r_bloodgroup}'
        )

    def __repr__(self) -> str:
        return (
            f'Patient {self.id_recipient}, listed on '
            f'{self.__dict__[cn.LISTING_DATE].strftime("%d-%m-%Y")}: '
            f'at center {self.__dict__[cn.RECIPIENT_CENTER]}'
            f' with bloodgroup {self.r_bloodgroup}\n'
        )


class Profile:
    """Class which implements an obligation

    Attributes   #noqa
    ----------
    min_age : int
        Minimum acceptable age.
    max_age : int
        Maximum acceptable age.
    hbsag : bool
        Hepatitis B surface antigen status.
    hcvab : bool
        Hepatitis C antibody status.
    hbcab : bool
        Hepatitis B core antibody status.
    sepsis : bool
        Sepsis status.
    meningitis : bool
        Meningitis status.
    malignancy : bool
        Malignancy status.
    drug_abuse : bool
        Drug abuse status.
    rescue : bool
        Rescue status.
    euthanasia : bool
        Euthanasia status.
    dcd : bool
        Donor after circulatory death status.
    match_qualities : Dict[Tuple[int, int, int], bool]
        Dictionary representing match qualities.
    esp : bool
        ESP status.
    extended_esp : bool
        Extended ESP status.

    Methods
    -------
    _check_acceptable(don: Donor, verbose: bool = False) -> bool
        Check whether donor is acceptable according to profile.
    _check_hla_acceptable(mq: Optional[Tuple[int, ...]] = None,
        verbose: bool = False) -> bool
        Check whether HLA is acceptable according to profile.
    """
    __slots__ = [
        'min_age', 'max_age', 'hbsag',
        'hcvab', 'hbcab', 'sepsis', 'meningitis', 'malignancy',
        'drug_abuse', 'rescue', 'euthanasia', 'dcd', 'match_qualities',
        'esp', 'extended_esp'
    ]

    def __init__(
        self, min_age: int, max_age: int,
        hbsag: bool, hcvab: bool, hbcab: bool,
        sepsis: bool, meningitis: bool,
        malignancy: bool, drug_abuse: bool,
        rescue: bool, euthanasia: bool, dcd: bool,
        match_qualities: Dict[Tuple[int, int, int], bool],
        esp: bool, extended_esp: bool
    ) -> None:

        self.min_age = min_age
        self.max_age = max_age
        self.hbsag = hbsag
        self.hcvab = hcvab
        self.hbcab = hbcab
        self.sepsis = sepsis
        self.meningitis = meningitis
        self.malignancy = malignancy
        self.drug_abuse = drug_abuse
        self.rescue = rescue
        self.euthanasia = euthanasia
        self.dcd = dcd
        self.match_qualities = match_qualities
        self.esp = esp
        self.extended_esp = extended_esp

    def _check_acceptable(
        self, don: Donor, verbose=False
    ) -> bool:
        """Check whether donor is acceptable to patient."""

        don_dict = don.__dict__
        if don_dict[cn.D_DCD] > self.dcd:
            if verbose:
                print('DCD-incompatible')
            return False
        if don_dict[cn.D_AGE] < self.min_age:
            if verbose:
                print(f'{don_dict[cn.D_AGE]} >= {self.min_age}')
            return False
        if don_dict[cn.D_AGE] > self.max_age:
            if verbose:
                print(f'Donor age {don_dict[cn.D_AGE]} > {self.max_age}')
            return False
        if don.__dict__[cn.D_HBSAG] > self.hbsag:
            if verbose:
                print('HBsag incompatible')
            return False
        if don.__dict__[cn.D_HCVAB] > self.hcvab:
            if verbose:
                print('HCVab incompatible')
            return False
        if don.__dict__[cn.D_HBCAB] > self.hbcab:
            if verbose:
                print('HBCab incompatible')
            return False
        if don.sepsis > self.sepsis:
            if verbose:
                print('sepsis incompatible')
            return False
        if don.meningitis > self.meningitis:
            if verbose:
                print('meningitis incompatible')
            return False
        if don.malignancy > self.malignancy:
            if verbose:
                print('malignancy incompatible')
            return False
        if don.rescue > self.rescue:
            if verbose:
                print('rescue incompatible')
            return False
        if don.euthanasia > self.euthanasia:
            if verbose:
                print('Euthanasia incompatible')
            return False
        esp_donor = don_dict[cn.D_AGE] >= don.age_esp_eligible
        if esp_donor > self.esp:
            if verbose:
                print('Check whether ESP is accepted')
            return False
        if esp_donor and self.extended_esp == 0:
            if verbose:
                print('Check whether extended ESP donor is accepted')
            return False

        return True

    def _check_hla_acceptable(
        self, mq: Optional[Tuple[int, int, int]] = None, verbose=False
    ) -> bool:
        if mq is not None and self.match_qualities[mq] == 0:
            if verbose:
                print('HLA-incompatible')
            return False
        elif mq is not None:
            return True
        return False

    def __str__(self) -> str:
        return str(
            {slot: getattr(self, slot) for slot in self.__slots__}
        )
