#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 17:33:44 2022

Scripts to read in input files.

@author: H.C. de Ferrante
"""

from datetime import datetime, timedelta
from typing import Optional, List, Dict, Any, Tuple, Callable
from collections import defaultdict
from numpy import ndarray, nan, array, isnan
import pandas as pd

import yaml

from simulator.code.utils.utils import DotDict
import simulator.magic_values.inputfile_settings as dtypes
import simulator.magic_values.column_names as cn
import simulator.magic_values.column_groups as cg
import simulator.magic_values.magic_values_rules as mgr
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
from simulator.code.matchlist.ScoringFunction import MatchPointFunction


def _read_with_datetime_cols(
        input_path: str,
        dtps: Dict[str, Any],
        casecols: bool = False,
        usecols: Optional[List[str]] = None,
        datecols: Optional[List[str]] = None,
        **kwargs
) -> pd.DataFrame:
    """Read in pd.DataFrame with datecols as datetime"""

    if usecols is None:
        usecols = list(dtps.keys())
    if casecols:
        usecols = [c.upper() for c in usecols]
        dtps = {k.upper(): v for k, v in dtps.items()}
        if datecols:
            datecols = [c.upper() for c in datecols]

    data_ = pd.read_csv(
        input_path,
        dtype=dtps,
        usecols=lambda x: x in usecols,
        **kwargs
    )
    if datecols:
        for col in datecols:
            col_format = (
                "%Y-%m-%d"
                if all(data_[col].str.len() == 10)
                else "%Y-%m-%d %H:%M:%S"
            )
            data_[col] = pd.to_datetime(
                data_[col],
                format=col_format
            )
            data_[col] = pd.Series(
                array(data_[col].values.astype('datetime64[ns]')),
                dtype='object'
            )

    assert isinstance(data_, pd.DataFrame), \
        f'Expected DataFrame, not {type(data_)}'

    if casecols:
        data_.columns = data_.columns.str.lower()

    return data_


def read_hla_match_table(input_path: str):
    d_hla_ = pd.read_csv(input_path)
    return d_hla_


def fix_hla_string(str_col: pd.Series) -> pd.Series:
    """Harmonize a HLA string column being read in"""

    for forbidden_character in es.HLA_FORBIDDEN_CHARACTERS:
        str_col = str_col.replace(forbidden_character, '', regex=True)
    str_col = str_col.str.upper()
    str_col = str_col.replace(
        to_replace={
            '(?<=\\s)RB': 'DRB',
            '(?<=\\s)QB': 'DQB',
            'CW(?=([A-Z]|[0-9]){4})': 'C',
            '/[0-9]+': '',
            '\\s\\s': ' '
        },
        regex=True
    )
    return str_col


def read_patients(
    input_path: str,
    datecols: Optional[List[str]] = None,
    usecols: Optional[List[str]] = None,
    **kwargs
) -> pd.DataFrame:
    """"Read in patient file."""

    if datecols is None:
        datecols = [cn.TIME_REGISTRATION, cn.R_DOB, cn.DATE_FIRST_DIAL]
    if usecols is None:
        usecols = list(dtypes.DTYPE_PATIENTLIST.keys())

    data_ = _read_with_datetime_cols(
        input_path=input_path,
        dtps=dtypes.DTYPE_PATIENTLIST,
        usecols=usecols,
        datecols=datecols,
        **kwargs
    )
    assert isinstance(data_, pd.DataFrame), \
        f'Expected DataFrame, not {type(data_)}'

    # Sort recipients by id & time of registration.
    data_: pd.DataFrame = data_.sort_values(
        by=[cn.ID_RECIPIENT, cn.TIME_REGISTRATION]
    )
    data_.reset_index(inplace=True)

    # Add re-transplantation information
    idx = data_[cn.TIME_SINCE_PREV_TXP].notna()

    data_.loc[idx, cn.PREV_TX_DATE] = (
        pd.to_datetime(data_.loc[idx, cn.TIME_REGISTRATION]) -
        pd.to_timedelta(data_.loc[idx, cn.TIME_SINCE_PREV_TXP], unit='days')
    )

    return data_


def read_rescue_baseline_hazards(
    input_path: str
) -> Tuple[Dict[str, ndarray], Tuple[str]]:
    """Read in rescue parameters for Cox PH model"""
    rescue_probs = pd.read_csv(
        input_path,
        delimiter=','
    )
    if 'strata' in rescue_probs.columns:
        rescue_probs = {
            n: x.iloc[:, 0:2] for n, x in rescue_probs.groupby('strata')
        }
        rescue_prob_dict = defaultdict(dict)
        strata_vars = None
        for k, v in rescue_probs.items():
            if ',' in k:
                k1, k2 = k.split(',')
                bh_var1, bh_var1_level = k1.split('=')
                bh_var2, bh_var2_level = k2.split('=')
                rescue_prob_dict[bh_var1_level].update(
                    {
                        bh_var2_level: {
                            cn.N_OFFERS_TILL_RESCUE:
                            v['offers_before_rescue'].to_numpy(),
                            cn.CBH_RESCUE: v['hazard'].to_numpy()
                        }
                    }
                )
                if not strata_vars or len(strata_vars) == 1:
                    strata_vars = (bh_var1.strip(), bh_var2.strip())
            else:
                bh_var, bh_level = k.split('=')

                rescue_prob_dict[bh_level] = {
                    cn.N_OFFERS_TILL_RESCUE:
                    v['offers_before_rescue'].to_numpy(),
                    cn.CBH_RESCUE: v['hazard'].to_numpy()
                }
                if not strata_vars:
                    strata_vars = (bh_var, )

        return rescue_prob_dict, strata_vars
    else:
        rescue_probs = {
            nan: {
                cn.N_OFFERS_TILL_RESCUE:
                rescue_probs.loc[:, 'offers_before_rescue'].to_numpy(),
                cn.CBH_RESCUE: rescue_probs.loc[:, 'hazard'].to_numpy()
            }
        }
    return rescue_probs, None


def read_donors(
        input_path: str,
        datecols: Optional[List[str]] = None,
        usecols: Optional[List[str]] = None,
        **kwargs
) -> pd.DataFrame:
    """"Read in patient file."""

    if usecols is None:
        usecols = list(dtypes.DTYPE_DONORLIST.keys())
    if datecols is None:
        datecols = [cn.D_DATE]

    data_: pd.DataFrame = _read_with_datetime_cols(
        input_path=input_path,
        dtps=dtypes.DTYPE_DONORLIST,
        casecols=True,
        usecols=usecols,
        datecols=datecols,
        **kwargs
    )

    for col in (cn.DONOR_HLA, ):
        data_[col] = fix_hla_string(data_[col])

    # Remove donors not from ET
    data_ = data_[data_[cn.D_COUNTRY].isin(es.ET_COUNTRIES)]

    data_ = data_.fillna(
        value=dtypes.DTYPE_DONOR_FILL_NAS
    )

    return data_


def read_donor_balances(
        input_path: str,
        datecols: Optional[List[str]] = None,
        usecols: Optional[List[str]] = None,
        **kwargs
) -> pd.DataFrame:
    """"Read in patient file."""

    if usecols is None:
        usecols = list(dtypes.DTYPE_DONORBALLIST.keys())
    if datecols is None:
        datecols = [cn.D_DATE]

    data_: pd.DataFrame = _read_with_datetime_cols(
        input_path=input_path,
        dtps=dtypes.DTYPE_DONORBALLIST,
        casecols=True,
        usecols=usecols,
        datecols=datecols,
        **kwargs
    )

    # Remove donors not from ET
    data_ = data_[data_[cn.D_ALLOC_COUNTRY].isin(es.ET_COUNTRIES)]
    data_ = data_.fillna(
        value=dtypes.DTYPE_DONOR_FILL_NAS
    )
    return data_


def read_historic_donor_balances(
        input_path: str,
        sim_start_date: datetime,
        max_window_length: 99999,
        datecols: Optional[List[str]] = None,
        usecols: Optional[List[str]] = None,
        **kwargs
):
    data_ = read_donor_balances(
        input_path=input_path,
        usecols=usecols,
        datecols=datecols,
        **kwargs
    )

    # Filter to donors which occured prior to simulation start date
    data_ = data_[data_[cn.D_DATE] <= sim_start_date]
    data_ = data_[
        data_[cn.D_DATE] >= (
            sim_start_date - timedelta(days=max_window_length)
        )
    ]

    # Add donor age group
    data_[cn.DONOR_BALANCE_AGE_GROUP] = pd.cut(
        data_['donor_age'],
        bins=[-float('inf'), 17, 49, 64, float('inf')],
        labels=[
            mgr.LT18, mgr.YOUNGADULT, mgr.OLDADULT, mgr.ELDERLY
        ],
        right=True  # Include the right edge in intervals
    )

    return data_


def read_nonetkasesp_balances(
        input_path: str,
        sim_start_date: datetime,
        sim_end_date: datetime,
        datecols: Optional[List[str]] = None,
        usecols: Optional[List[str]] = None,
        **kwargs
) -> pd.DataFrame:
    """"Read in patient file."""

    data_ = read_donor_balances(
        input_path=input_path,
        usecols=usecols,
        datecols=datecols,
        **kwargs
    )

    # Filter to donors which occured during the sim period
    data_ = data_[data_[cn.D_DATE] >= sim_start_date]
    data_ = data_[data_[cn.D_DATE] <= sim_end_date]

    # Add donor age group
    data_[cn.DONOR_BALANCE_AGE_GROUP] = pd.cut(
        data_['donor_age'],
        bins=[-float('inf'), 17, 49, 64, float('inf')],
        labels=[
            mgr.LT18, mgr.YOUNGADULT, mgr.OLDADULT, mgr.ELDERLY
        ],
        right=True  # Include the right edge in intervals
    )

    # Filter to non-ETKAS/ESP donors.
    data_ = data_[data_[cn.NON_ETKAS_ESP] == 1]

    data_ = data_.fillna(
        value=dtypes.DTYPE_DONOR_FILL_NAS
    )

    return data_


def read_donor_pool(
        input_path: str,
        datecols: Optional[List[str]] = None,
        usecols: Optional[List[str]] = None,
        **kwargs
) -> pd.DataFrame:
    """"Read in patient file."""

    if usecols is None:
        usecols = list(dtypes.DTYPE_DONORLIST.keys())

    data_: pd.DataFrame = _read_with_datetime_cols(
        input_path=input_path,
        dtps=dtypes.DTYPE_DONORLIST,
        casecols=True,
        usecols=usecols,
        datecols=datecols,
        **kwargs
    )

    data_[cn.DONOR_HLA] = fix_hla_string(
        data_[cn.DONOR_HLA]
    )

    data_ = data_.fillna(
        value=dtypes.DTYPE_DONOR_FILL_NAS
    )

    return data_


def read_hla_match_potentials(
    input_path: str,
    datecols: Optional[List[str]] = None,
    usecols: Optional[List[str]] = None,
    **kwargs
) -> pd.DataFrame:

    if usecols is None:
        usecols = list(dtypes.DTYPE_MATCH_POTENTIALS.keys())
    data_ = pd.read_csv(input_path)
    data_: pd.DataFrame = _read_with_datetime_cols(
        input_path=input_path,
        dtps=dtypes.DTYPE_MATCH_POTENTIALS,
        casecols=False,
        usecols=usecols,
        datecols=datecols,
        **kwargs
    )
    return data_


def read_travel_times(
    path_drive_time: str = es.PATH_DRIVING_TIMES
) -> Dict[str, Dict[str, float]]:
    """Read in expected travelling times by car between centers"""

    travel_info = pd.read_csv(
        path_drive_time
    )

    # Select the right travel time.
    travel_dict = (
        travel_info.groupby(cn.FROM_CENTER)[
            [cn.FROM_CENTER, cn.TO_CENTER, cn.DRIVING_TIME]
        ].apply(lambda x: x.set_index(cn.TO_CENTER).to_dict(orient='index'))
        .to_dict()
    )

    return travel_dict


def read_profiles(
    input_path: str,
    usecols: Optional[List[str]] = None,
    **kwargs
) -> pd.DataFrame:
    """"Read in patient file."""

    if usecols is None:
        usecols = list(dtypes.DTYPE_PROFILES.keys())

    data_ = _read_with_datetime_cols(
        input_path,
        dtps=dtypes.DTYPE_PROFILES,
        casecols=False,
        **kwargs
    )

    assert isinstance(data_, pd.DataFrame), \
        f'Expected DataFrame, not {type(data_)}'
    data_.columns = data_.columns.str.lower()

    return data_


def read_status_updates(
        input_path: str, sim_set: DotDict,
        start_date_col: Optional[str] = None,
        end_date_col: Optional[str] = None,
        **kwargs) -> pd.DataFrame:
    """Read in status updates."""

    d_s = _read_with_datetime_cols(
        input_path=input_path,
        dtps=dtypes.DTYPE_STATUSUPDATES,
        usecols=list(dtypes.DTYPE_STATUSUPDATES.keys()),
        datecols=[cn.LISTING_DATE],
        **kwargs
    )

    assert isinstance(d_s, pd.DataFrame), \
        f'Expected DataFrame, not {type(d_s)}'
    d_s.columns = d_s.columns.str.lower()

    end_date: datetime = (
        sim_set.SIM_END_DATE if not end_date_col
        else sim_set.__dict__[end_date_col]
    )
    d_s = d_s.loc[
        d_s[cn.LISTING_DATE] <= end_date,
        :
    ]

    # Set variable detail to removal reason for patients removed.
    d_s.loc[
        d_s[cn.STATUS_VALUE] == 'R',
        cn.STATUS_DETAIL
    ] = d_s.loc[
        d_s[cn.STATUS_VALUE] == 'R',
        cn.REMOVAL_REASON
    ]

    # Fix hla strings
    for type_upd in (mgr.HLA, mgr.UNACC):
        d_s.loc[d_s[cn.TYPE_UPDATE] == type_upd, cn.STATUS_VALUE] = (
            fix_hla_string(
                d_s.loc[
                    d_s[cn.TYPE_UPDATE] == type_upd,
                    cn.STATUS_VALUE]
            )
        )

    # Remove HLA updates without an A, B, and DR
    mask_invalid_hla = (
        (d_s.loc[:, cn.TYPE_UPDATE] == mgr.HLA) &
        (
            ~d_s.loc[:, cn.STATUS_VALUE].str.
            contains('A.*B.*DR.*', regex=True, na=True)
        )
    )
    print(f'Removed {sum(mask_invalid_hla)} invalid HLA status updates')
    d_s = d_s.loc[~mask_invalid_hla, :]

    # Remove multiple updates occuring at the same time-stamp.
    # Their order will not be maintained by the heapqueue
    d_s = d_s.drop_duplicates(
        subset=[cn.ID_REGISTRATION, cn.TSTART, cn.TYPE_UPDATE],
        keep='last'
    )

    return d_s


def pediatric_donor_function_factory(age_thresh: int) -> Callable:
    def check_etkas_ped_don(donor: 'entities.Donor') -> bool:
        """Checks whether a donor is pediatric (i.e. under age 18)"""
        if not isnan(donor.__dict__[cn.D_AGE]):
            return donor.__dict__[cn.D_AGE] < age_thresh
        return False
    return check_etkas_ped_don


def pediatric_patient_function_factory(
    age_thresh: int, dial_req: bool
) -> Callable:
    def check_etkas_ped_rec(
        candidate_country: str,
        age_at_listing: float,
        match_age: float,
        age_first_dial: Optional[float],
        prev_txp_ped: Optional[bool],
        time_since_prev_txp: Optional[float]
    ) -> bool:
        """Checks whether a patient is pediatric. Can either be:
            -   pediatric age at match (<18)
            -   pediatric age at listing (non-DE)
            -   on dialysis within 91 days of failed pediatric
                transplantation (non-DE)

            Until March 2023, the patient had to have started dialysis
            before age 18. After March 2023, being listed before age
            18 is enough.
        """
        if not isnan(match_age) and match_age < age_thresh:
            return True
        if candidate_country != mgr.GERMANY:
            if not dial_req and age_at_listing and age_at_listing < age_thresh:
                return True
            elif (
                dial_req and age_first_dial and
                age_first_dial < age_thresh
            ):
                return True
            elif (
                prev_txp_ped and time_since_prev_txp and
                time_since_prev_txp <= age_thresh
            ):
                if (
                    age_at_listing - time_since_prev_txp < age_thresh and
                    age_first_dial
                ):
                    if age_first_dial <= match_age:
                        return True
                    else:
                        return False
        return False
    return check_etkas_ped_rec


def read_sim_settings(
        ss_path: str,
        date_settings: Optional[List[str]] = None
) -> DotDict:
    """Read in simulation settings"""
    with open(ss_path, "r", encoding='utf-8') as file:
        sim_set: Dict[str, Any] = yaml.load(file, Loader=yaml.FullLoader)

    if date_settings is None:
        date_settings = [
            'SIM_END_DATE', 'SIM_START_DATE',
            'LOAD_RETXS_TO', 'LOAD_RETXS_FROM'
        ]

    # Fix dates
    min_time = datetime.min.time()
    for k in date_settings:
        sim_set[k] = datetime.combine(
            sim_set[k], min_time
        )

    sim_set['calc_etkas_score'] = MatchPointFunction(
        coef=sim_set[cn.POINTS_ETKAS.upper()],
        points_comp_to_group=es.POINT_COMPONENTS_TO_GROUP,
        point_multiplier=sim_set.get(cn.MULTIPLIER_ETKAS.upper()),
        trafos=sim_set.get(cn.TRAFOS_ETKAS.upper())
    )

    sim_set['calc_esp_score'] = MatchPointFunction(
        coef=sim_set[cn.POINTS_ESP.upper()],
        points_comp_to_group=es.POINT_COMPONENTS_TO_GROUP
    )

    # Read in geographic definitions from yml file.
    with open(
        sim_set['PATH_GEOGRAPHIC_DEFS'], "r",
        encoding='utf-8'
    ) as file:
        for k, v in yaml.load(file, Loader=yaml.FullLoader).items():
            sim_set[k] = v
    sim_set['ALLOWED_REGIONS'] = set(
        sim_set['CENTERS_TO_REGIONS'].values()
    )

    # Read in ESP match tiers from yml file.
    with open(
        sim_set['PATH_ESP_TIERS'], "r",
        encoding='utf-8'
    ) as file:
        try:
            sim_set['ESP_MATCH_TIERS'] = {
                country: tuple((k[0], tuple(v)) for k, v in innerdict.items())
                for country, innerdict in
                yaml.load(file, Loader=yaml.FullLoader).items()
            }
        except Exception as e:
            print(f"Error loading YAML: {e}")

    # Read in ESP match tiers from yml file.
    with open(
        sim_set['PATH_ETKAS_TIERS'], "r",
        encoding='utf-8'
    ) as file:
        try:
            sim_set['ETKAS_MATCH_TIERS'] = tuple(
                (
                    k[0], tuple(k[1])
                ) for k in yaml.load(file, Loader=yaml.FullLoader).items()
            )
        except Exception as e:
            print(f"Error loading YAML: {e}")

    sim_set['determine_mismatchfreqs'] = (
        (
            len(sim_set['calc_etkas_score'].variables.
                intersection(cg.HLA_MISMATCH_FREQS))) |
        (
            len(sim_set['calc_esp_score'].variables.
                intersection(cg.HLA_MISMATCH_FREQS))
        )
    )

    sim_set['times_esp_eligible'] = {
        k: (v - sim_set['SIM_START_DATE']) / timedelta(days=1)
        for k, v in es.ESP_CHOICE_ENFORCED.items()
    }

    # Default age at which donors/candidates are ESP
    if 'DONOR_AGE_ESP_ELIGIBLE' not in sim_set:
        sim_set['DONOR_AGE_ESP_ELIGIBLE'] = {
            cntry: 65 for cntry in es.ET_COUNTRIES
        }
    if 'CAND_AGE_ESP_ELIGIBLE' not in sim_set:
        sim_set['CAND_AGE_ESP_ELIGIBLE'] = {
            cntry: 65 for cntry in es.ET_COUNTRIES
        }

    # Pediatric definitions from March 2021
    sim_set['check_d_pediatric'] = pediatric_donor_function_factory(
        sim_set['PEDIATRIC_DONOR_AGE']
    )
    sim_set['check_r_pediatric'] = pediatric_patient_function_factory(
        sim_set['PEDIATRIC_CANDIDATE_AGE'],
        sim_set.get('PEDIATRIC_REQUIRES_DIALYSIS', False)
    )

    return DotDict(sim_set)
