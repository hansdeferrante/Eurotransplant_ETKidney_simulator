#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 17:33:44 2022

Scripts to read in input files.

@author: H.C. de Ferrante
"""

from typing import List, Dict, Dict, Tuple, Any, Iterator, Hashable, Optional
from datetime import timedelta
from itertools import groupby
import pandas as pd
import numpy as np
import warnings

import simulator.magic_values.etkidney_simulator_settings as es
from simulator.code.utils.utils import DotDict
import simulator.magic_values.column_names as cn
import simulator.magic_values.magic_values_rules as mgr
import simulator.code.utils.read_input_files as rdr
from simulator.code.events.StatusUpdate import StatusUpdate, ProfileUpdate
from simulator.code.events.PatientStatusQueue import PatientStatusQueue
from simulator.code.entities import (
    Donor, Patient, Profile
)
from simulator.code.BalanceSystem import BalanceSystem
from simulator.code.HLA.HLASystem import HLASystem
from simulator.code.HLA.MMPSystem import MMPSystem


def _read_patients_rich(
        sim_set: DotDict,
        start_date_col: Optional[str] = None,
        end_date_col: Optional[str] = None,
        nonnacols: Tuple[str, ...] = (cn.TIME_TO_DEREG, cn.PATIENT_HLA),
        **kwargs
) -> pd.DataFrame:
    """
        Read patients in with whether it is a retransplantation
        between the start and end date times.
    """

    # Read in patient data
    d_patients: pd.DataFrame = rdr.read_patients(
        sim_set.PATH_PATIENTS, **kwargs
    )
    d_patients.loc[:, cn.PATIENT_HLA] = rdr.fix_hla_string(
        d_patients.loc[:, cn.PATIENT_HLA]
    )

    # Read only patient data in before end of simulation.
    end_date = (
        sim_set.SIM_END_DATE if end_date_col is None
        else sim_set.__dict__[end_date_col]
    )
    d_patients = d_patients.loc[
        d_patients[cn.LISTING_DATE] <= end_date,
        :
    ]

    # Add whether patient is re-transplanted.
    d_patients[cn.RETRANSPLANT] = d_patients[cn.NTH_TRANSPLANT] > 0
    idx = d_patients.loc[:, cn.RETRANSPLANT]
    d_patients.loc[idx, cn.RETRANSPLANT_DURING_SIM] = (
        d_patients.loc[idx, cn.PREV_TX_DATE] >= sim_set.SIM_START_DATE
    )
    d_patients.loc[~idx, cn.RETRANSPLANT_DURING_SIM] = False

    # Add column for type re-transplantation
    d_patients[cn.TYPE_RETX] = np.select(
        condlist=[
            d_patients[cn.RETRANSPLANT_DURING_SIM].astype(bool).to_numpy(),
            d_patients[cn.RETRANSPLANT].to_numpy()
        ],
        choicelist=[
            cn.RETRANSPLANT_DURING_SIM,
            cn.RETRANSPLANT
        ],
        default=cn.NO_RETRANSPLANT
    )

    # Remove patients deregistered before simulation start
    start_date = (
        sim_set.SIM_START_DATE if start_date_col is None
        else sim_set.__dict__[start_date_col]
    )

    # Filter to selected period
    d_patients = d_patients.loc[
        (
            d_patients.apply(
                lambda x:
                x[cn.LISTING_DATE] + timedelta(days=x[cn.TIME_TO_DEREG]),
                axis=1
            ) >= start_date
        ),
        :
    ]

    # Remove patients with missing required info
    # Drop patients with NAs for columns not allowed to be nan
    for nonnacol in nonnacols:
        patients_with_nas = d_patients.loc[:, nonnacol].isna().sum()
        if patients_with_nas:
            print(
                f'Dropping {patients_with_nas} patients '
                f'with unknown {nonnacol}'
            )
            d_patients = d_patients.loc[
                ~ d_patients[nonnacol].isna(),
                :
            ]

    return d_patients


def _rcrd_to_patient(
    sim_set: DotDict, hla_system: HLASystem,
    mmp_system: MMPSystem, rcrd: dict[Hashable, Any]
) -> Patient:
    """Convert a record to a Patient object"""
    return Patient(
        id_recipient=rcrd[cn.ID_RECIPIENT],
        recipient_country=rcrd[cn.RECIPIENT_COUNTRY],
        recipient_center=rcrd[cn.RECIPIENT_CENTER],
        recipient_region=rcrd[cn.RECIPIENT_REGION],
        bloodgroup=rcrd[cn.R_BLOODGROUP],
        listing_date=rcrd[cn.TIME_REGISTRATION],
        urgency_code='NT',
        r_dob=rcrd[cn.R_DOB],
        sim_set=sim_set,
        sex=rcrd[cn.PATIENT_SEX],
        type_retx=rcrd[cn.TYPE_RETX],
        time_since_prev_txp=rcrd[cn.TIME_SINCE_PREV_TXP],
        prev_txp_living=rcrd[cn.PREV_TXP_LIVING],
        prev_txp_ped=rcrd[cn.PREV_TXP_PED],
        id_reg=rcrd[cn.ID_REGISTRATION],
        seed=sim_set.SEED,
        hla=rcrd[cn.PATIENT_HLA],
        hla_system=hla_system,
        mmp_system=mmp_system,
        previous_wt=rcrd[cn.PREVIOUS_T],
        kidney_program=rcrd[cn.WLKI_PROGRAMME],
        date_first_dial=rcrd[cn.DATE_FIRST_DIAL]
    )


def load_patients(
    sim_set: DotDict, hla_system: HLASystem,
    mmp_system: MMPSystem, **kwargs
) -> dict[int, Patient]:
    """Load list of patients"""

    # Load patients. The start and end date
    # are simulation start and end by default.
    d_patients = _read_patients_rich(
        sim_set=sim_set,
        **kwargs
    )

    # If we simulate re-transplantations, do not load patient registrations
    # that are future re-transplantations for patients not yet transplanted.
    patient_dict = {
        rcrd[cn.ID_REGISTRATION]: _rcrd_to_patient(
            sim_set=sim_set, rcrd=rcrd, hla_system=hla_system,
            mmp_system=mmp_system
        ) for rcrd in d_patients.to_dict(orient='records')
        if (
            not sim_set.SIM_RETX or not (
                sim_set.SIM_RETX and
                rcrd[cn.TYPE_RETX] == cn.RETRANSPLANT_DURING_SIM and
                (
                    pd.isna(rcrd[cn.PREV_TXP_LIVING]) or
                    rcrd[cn.PREV_TXP_LIVING] == 0
                )
            )
        )
    }

    return patient_dict


def load_retransplantations(
    sim_set: DotDict, hla_system: HLASystem, mmp_system: MMPSystem, **kwargs
) -> dict[int, Patient]:
    """Load list of patients"""

    # Load patients. The start and end date
    # are simulation start and end by default.
    d_patients = _read_patients_rich(
        sim_set=sim_set
    )

    # Do not include patients listed outside the transplantation registration
    # window; gives bias.
    d_patients = d_patients.loc[
        (
            d_patients.loc[:, cn.LISTING_DATE] >=
            sim_set.__dict__['LOAD_RETXS_FROM']
        ) & (
            d_patients.loc[:, cn.LISTING_DATE] <=
            sim_set.__dict__['LOAD_RETXS_TO']
        ),
        :
    ]

    # Construct a list of patients from patient data
    patient_dict = {
        rcrd[cn.ID_REGISTRATION]: _rcrd_to_patient(
            sim_set=sim_set, hla_system=hla_system,
            rcrd=rcrd, mmp_system=mmp_system
        ) for rcrd in d_patients.to_dict(orient='records')
        if rcrd[cn.TYPE_RETX] != cn.NO_RETRANSPLANT
    }

    return patient_dict


def load_balances(
        sim_set: DotDict,
        update_balances: bool = True
) -> BalanceSystem:
    """Load list of donors"""

    # Read in to be balanced transplantations.
    # Select only transplantations inserted before the
    # simulation start date which did not end yet.
    d_bal = rdr.read_historic_donor_balances(
        sim_set.PATH_BALANCES,
        sim_start_date=sim_set.SIM_START_DATE,
        max_window_length=sim_set.WINDOW_FOR_BALANCE
    )

    if sim_set.REMOVE_BALANCES_BEFORE_STARTTIME:
        d_bal = d_bal.loc[
            d_bal.d_date >= pd.Timestamp(sim_set.STARTTIME_FOR_BALANCE)
        ]

    # Add tstart and tstop columns for existing balances
    balance_system = BalanceSystem.from_balance_df(
        ss=sim_set,
        df_init_balances=d_bal,
        group_vars=sim_set.BALANCE_GROUP_VARS,
        update_balances=update_balances
    )

    return balance_system


def load_nonetkasesp_balances(sim_set: DotDict) -> Dict[int, Dict[str, Any]]:
    """Load non-ETKAS/ESP balances"""

    d_bal = rdr.read_nonetkasesp_balances(
        sim_set.PATH_BALANCES,
        sim_start_date=sim_set.SIM_START_DATE,
        sim_end_date=sim_set.SIM_END_DATE
    )

    bal_dict = {
        i: rcrd for i, rcrd in enumerate(d_bal.to_dict('records'))
    }

    return bal_dict


def load_donors(
        sim_set: DotDict,
        hla_system: HLASystem,
        **kwargs
) -> dict[int, Donor]:
    """Load list of donors"""

    # Read in donor data data
    d_don = rdr.read_donors(sim_set.PATH_DONORS, **kwargs)
    d_don = d_don.loc[
        (d_don[cn.D_DATE] >= sim_set.SIM_START_DATE) &
        (d_don[cn.D_DATE] <= sim_set.SIM_END_DATE),
        :
    ]

    # Construct a list of patients from patient data
    donor_dict = {
        rcrd[cn.ID_DONOR]: Donor(
            sim_set=sim_set,
            id_donor=rcrd[cn.ID_DONOR],
            n_kidneys_available=rcrd[cn.N_KIDNEYS_AVAILABLE],
            bloodgroup=rcrd[cn.D_BLOODGROUP],
            donor_country=rcrd[cn.D_COUNTRY],
            donor_region=rcrd[cn.D_REGION],
            donor_subregion_esp=rcrd[cn.D_ESP_SUBREGION],
            donor_center=rcrd[cn.D_CENTER],
            reporting_date=rcrd[cn.D_DATE],
            donor_dcd=rcrd[cn.D_DCD],
            age=rcrd[cn.D_AGE],
            cmv=rcrd[cn.D_CMV],
            hbsag=rcrd[cn.D_HBSAG],
            hcvab=rcrd[cn.D_HCVAB],
            hbcab=rcrd[cn.D_HBCAB],
            sepsis=rcrd[cn.D_SEPSIS],
            meningitis=rcrd[cn.D_MENINGITIS],
            malignancy=rcrd[cn.D_MALIGNANCY],
            drug_abuse=rcrd[cn.D_DRUG_ABUSE],
            euthanasia=rcrd[cn.D_EUTHANASIA],
            tumor_history=rcrd[cn.D_TUMOR_HISTORY],
            donor_marginal_free_text=rcrd[cn.D_MARGINAL_FREE_TEXT],
            death_cause_group=rcrd[cn.DEATH_CAUSE_GROUP],
            diabetes=rcrd[cn.D_DIABETES],
            hla_system=hla_system,
            hla=rcrd[cn.DONOR_HLA],
            hypertension=rcrd[cn.D_HYPERTENSION],
            last_creat=rcrd[cn.D_LAST_CREAT],
            urine_protein=rcrd[cn.D_URINE_PROTEIN],
            smoker=rcrd[cn.D_SMOKING],
            cardiac_arrest=rcrd[cn.D_CARREST],
            rescue=(
                False if sim_set.get('SIMULATE_RESCUE', False)
                else bool(rcrd[cn.D_RESCUE])
            )
        ) for rcrd in d_don.to_dict(orient='records')
    }

    return donor_dict


def preload_profiles(
        patients: Dict[int, Patient],
        sim_set: DotDict, **kwargs
) -> None:
    """Add profile information for patients"""

    # Read in donor data data
    d_profiles = rdr.read_profiles(sim_set.PATH_PROFILES, **kwargs)
    # For patient in patients
    for rcrd in d_profiles.to_dict(orient='records'):
        if (
            rcrd[cn.ID_REGISTRATION] in patients.keys() and
            patients[rcrd[cn.ID_REGISTRATION]].future_statuses is not None
        ):
            mqs = {
                k: rcrd[v] for k, v in es.PROFILE_HLA_MQS.items()
            }
            patients[rcrd[cn.ID_REGISTRATION]]. \
                future_statuses.add(
                    ProfileUpdate(
                        type_status=mgr.PRF,
                        arrival_time=rcrd[cn.TSTART],
                        profile=Profile(
                            min_age=rcrd[cn.PROFILE_MIN_DONOR_AGE],
                            max_age=rcrd[cn.PROFILE_MAX_DONOR_AGE],
                            hbsag=rcrd[cn.PROFILE_HBSAG],
                            hcvab=rcrd[cn.PROFILE_HCVAB],
                            hbcab=rcrd[cn.PROFILE_HBCAB],
                            sepsis=rcrd[cn.PROFILE_SEPSIS],
                            meningitis=rcrd[cn.PROFILE_MENINGITIS],
                            malignancy=rcrd[cn.PROFILE_MALIGNANCY],
                            drug_abuse=rcrd[cn.PROFILE_DRUG_ABUSE],
                            rescue=rcrd[cn.PROFILE_RESCUE],
                            euthanasia=rcrd[cn.PROFILE_EUTHANASIA],
                            dcd=rcrd[cn.PROFILE_DCD],
                            match_qualities=mqs,
                            esp=rcrd[cn.PROFILE_ESP],
                            extended_esp=rcrd[cn.PROFILE_EXTENDED_ESP]
                        ),
                        sim_start_time=(
                            patients[
                                rcrd[cn.ID_REGISTRATION]
                            ].__dict__[cn.LISTING_DATE] -
                            sim_set.SIM_START_DATE
                        ) / timedelta(days=1)
                    )
            )


def to_dict_fast(df):
    cols = list(df)
    col_arr_map = {col: df[col].astype(object).to_numpy() for col in cols}
    records = []
    for i in range(len(df)):
        record = {col: col_arr_map[col][i] for col in cols}
        records.append(record)
    return records


def preload_status_updates(
        patients: Dict[int, Patient],
        sim_set: DotDict,
        end_date_col: Optional[str] = None,
        **kwargs
) -> None:
    """ Read in status updates for selected patients,
        and preload them in the patients.
    """

    # Read in program updates for patients
    d_program_updates = (
        rdr.read_status_updates(
            sim_set.PATH_PROGRAM_UPDATES,
            end_date_col=end_date_col,
            sim_set=sim_set,
            casecols=True,
            **kwargs
        ).sort_values(
            by=[cn.ID_REGISTRATION, cn.TSTART, cn.STATUS_DETAIL],
            ascending=(True, True, False)
        )
    )
    d_program_updates = d_program_updates.loc[
        d_program_updates[cn.ID_REGISTRATION].isin(
            patients.keys()
        )
    ]

    # Read in status updates for patients
    d_status_updates = (
        rdr.read_status_updates(
            sim_set.PATH_STATUS_UPDATES,
            end_date_col=end_date_col,
            sim_set=sim_set,
            **kwargs
        ).sort_values(
            by=[cn.ID_REGISTRATION, cn.TSTART, cn.STATUS_DETAIL],
            ascending=(True, True, False)
        )
    )
    d_status_updates = d_status_updates.loc[
        d_status_updates[cn.ID_REGISTRATION].isin(
            patients.keys()
        )
    ]

    # Bind program updates to the status updates.
    d_status_updates = pd.concat(
        [
            d_program_updates.loc[
                d_program_updates.id_registration.isin(
                    d_status_updates.id_registration
                )
            ],
            d_status_updates
        ]
    )

    # Sort by status time, and make sure exits are
    # sorted last.
    d_status_updates['status_is_exit'] = np.where(
        d_status_updates.variable_value.isin(es.TERMINAL_STATUSES) &
        d_status_updates.type_update.isin({mgr.URG, f'S{mgr.URG}'}),
        1,
        0
    )

    d_status_updates = d_status_updates.sort_values(
        [cn.ID_REGISTRATION, cn.TSTART, 'status_is_exit'], na_position='first'
    )

    d_status_updates['already_exited'] = (
        d_status_updates.groupby(cn.ID_REGISTRATION)['status_is_exit']
        .cumsum()
        .groupby(d_status_updates[cn.ID_REGISTRATION])
        .shift(fill_value=0)
    )
    d_status_updates = (
        d_status_updates[d_status_updates['already_exited'] == 0]
        .drop(columns=['status_is_exit', 'already_exited'])
    )

    # Remove status updates which are real FU's
    d_status_updates = d_status_updates.loc[
        d_status_updates.loc[:, cn.TYPE_UPDATE] != 'SFU',
        :
    ]

    # Check whether imputation procedure was run correctly
    # (i.e. almost all patients have a terminal status)
    d_last_statuses = d_status_updates.groupby(
        cn.ID_REGISTRATION
    ).last().loc[:, cn.URGENCY_CODE]
    fb_last_statuses = (
        d_last_statuses.loc[~ d_last_statuses.isin([cn.R, cn.D])].index
    )
    if len(fb_last_statuses):
        d_status_updates = d_status_updates.loc[
            ~ d_status_updates.loc[:, cn.ID_REGISTRATION].isin(
                fb_last_statuses
            ),
            :
        ]
    if len(fb_last_statuses) > 5 and len(fb_last_statuses) < 25:
        warnings.warn(
            f'Removed {len(fb_last_statuses)} registrations '
            f'without R/D as last status update'
        )
    elif len(fb_last_statuses) > 5:
        raise ValueError(
            f'There are {len(fb_last_statuses)} registrations without '
            f'terminal R/D status. All registrations must end with R/D.'
        )

    d_status_updates[cn.SIM_START_TIME] = (
        d_status_updates.loc[:, cn.LISTING_DATE].values -
        sim_set.SIM_START_DATE
    ) / timedelta(days=1)

    for id_reg, status_upd in groupby(
        to_dict_fast(d_status_updates),
        lambda x: x[cn.ID_REGISTRATION]
    ):
        if int(id_reg) in patients:
            patients[int(id_reg)].preload_status_updates(
                fut_stat=dstat_to_queue(
                    rcrds=status_upd,
                    sim_set=sim_set
                )
            )


def dstat_to_queue(
        rcrds: Iterator[Dict[Hashable, Any]],
        sim_set: DotDict,
) -> PatientStatusQueue:
    """Read in statuses and convert to a heap queue instance"""

    status_queue = PatientStatusQueue(
        initial_events=[
            StatusUpdate(
                type_status=rcrd[cn.TYPE_UPDATE],
                arrival_time=rcrd[cn.TSTART],
                status_detail=rcrd[cn.STATUS_DETAIL],
                status_detail2=rcrd[cn.STATUS_DETAIL2],
                status_value=rcrd[cn.STATUS_VALUE],
                sim_start_time=rcrd[cn.SIM_START_TIME]
            ) for rcrd in rcrds
        ]
    )

    return status_queue
