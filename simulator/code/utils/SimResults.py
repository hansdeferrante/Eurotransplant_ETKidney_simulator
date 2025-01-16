#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

from functools import reduce
from datetime import datetime, timedelta
from typing import Tuple, Optional, Dict, List
import numpy as np
import os
import pandas as pd
import typing
import gzip
import shutil
import pathlib

from simulator.code.utils.utils import round_to_decimals
from simulator.code.utils.utils import DotDict
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.column_names as cn
from simulator.magic_values.inputfile_settings import DEFAULT_DATE_TIME_HMS
import simulator.magic_values.magic_values_rules as mgr

if typing.TYPE_CHECKING:
    from simulator.code import entities
    from simulator.code.matchlist import MatchList


class SimResults:
    """Class which implements a heapqueue

    Attributes   #noqa
    ----------
    transplantations: List[Dict[str, Any]]
        List of transplantations that occured
    exits: List[Dict[str, Any]]
        List of patient exits
    discards: List[Dict[str, Any]]
        List of discarded donors (not accepted by any)
    path_txp: str
        path for transplantations
    path_exits: str
        path for waitlist exits
    path_dsc: str
        path for discarded organs
    path_ml: str
        path for matchlists.
    ...
    """

    def __init__(
            self,
            cols_to_save_exit: Tuple[str, ...],
            cols_to_save_discard: Tuple[str, ...],
            cols_to_save_patients: Tuple[str, ...],
            cols_to_save_snapshot: Tuple[str, ...],
            sim_set: DotDict
    ):
        self.cols_to_save_exit = cols_to_save_exit
        self.cols_to_save_txp = cols_to_save_exit
        self.cols_to_save_discard = cols_to_save_discard
        self.cols_to_save_patients = cols_to_save_patients
        self.cols_to_save_snapshot = cols_to_save_snapshot
        self.transplantations = []
        self.posttransplant = []
        self.exits = []
        self.discards = []
        self.match_lists = []
        self.snapshots_waitlist = []

        self.path_txp: str = (
            sim_set.RESULTS_FOLDER + sim_set.PATH_TRANSPLANTATIONS
        )
        self.path_exits: str = sim_set.RESULTS_FOLDER + sim_set.PATH_EXITS
        self.path_dsc: str = sim_set.RESULTS_FOLDER + sim_set.PATH_DISCARDS
        self.path_pat: str = (
            sim_set.RESULTS_FOLDER + sim_set.PATH_FINAL_PATIENT_STATUS
        )
        self.path_ml: str = (
            sim_set.RESULTS_FOLDER + sim_set.PATH_MATCH_LISTS
        )
        self.sim_start_date = sim_set.SIM_START_DATE
        self.path_snapshot: str = (
            sim_set.RESULTS_FOLDER + sim_set.PATH_WAITLIST_SNAPSHOTS
        ) if sim_set.PATH_WAITLIST_SNAPSHOTS else None

        # Additional mismatches on loci to save.
        if sim_set.MATCHES_AT_TRANSPLANT_BROAD:
            self.broad_loci_to_save = set(sim_set.MATCHES_AT_TRANSPLANT_BROAD)
            self.cols_to_save_txp = self.cols_to_save_txp + tuple(
                'mmb_' + locus for locus in self.broad_loci_to_save
            )
        else:
            self.broad_loci_to_save = None
        if sim_set.MATCHES_AT_TRANSPLANT_SPLIT:
            self.split_loci_to_save = set(sim_set.MATCHES_AT_TRANSPLANT_SPLIT)
            self.cols_to_save_txp = self.cols_to_save_txp + tuple(
                'mms_' + locus for locus in self.split_loci_to_save
            )
        else:
            self.split_loci_to_save = None
        self.save_ml: bool = sim_set.SAVE_MATCH_LISTS

    def save_transplantation(
        self,
        pat: 'entities.Patient',
        matchr: 'MatchList.MatchRecord'
    ) -> None:
        """Save a transplantation to the transplantation list."""

        # Prepare potentially missing columns
        matchr.patient.get_pediatric(
            ped_fun=pat.sim_set.check_r_pediatric,
            match_age=matchr.__dict__[cn.R_MATCH_AGE]
        )

        # Combined dictionary
        result_dict = matchr.return_match_info(
            cols=self.cols_to_save_txp
        )

        result_dict[cn.TIME_WAITED] = round_to_decimals(
            (
                pat.__dict__[cn.EXIT_DATE] - pat.__dict__[cn.LISTING_DATE]
            ) / timedelta(days=1),
            2
        )

        result_dict[cn.ALLOCATION_PROGRAM] = (
            mgr.ESP if (
                matchr.__class__.__name__ == 'MatchRecordESP'
            ) else mgr.ETKAS
            if matchr.__class__.__name__ == 'MatchRecordETKAS'
            else None
        )
        result_dict[cn.EXT_ALLOC_TRIGGERED] = (
            1 if matchr.__dict__.get(cn.EXT_ALLOC_PRIORITY) is not None
            else 0
        )

        if self.broad_loci_to_save and pat.hla_system:
            result_dict.update(
                {
                    'mmb_' + k: v
                    for k, v in pat.hla_system._count_broad_mismatches(
                        d_hla=matchr.donor.hla,
                        p_hla=pat.hla,
                        loci=set(self.broad_loci_to_save)
                    ).items()
                }
            )
        if self.split_loci_to_save and pat.hla_system:
            result_dict.update(
                {
                    'mms_' + k: v
                    for k, v in pat.hla_system._count_split_mismatches(
                        d_hla=matchr.donor.hla,
                        p_hla=pat.hla,
                        loci=set(self.split_loci_to_save)
                    ).items()
                }
            )

        if result_dict[cn.ALLOCATION_PROGRAM] == mgr.ESP:
            result_dict[cn.ALLOCATION_MECHANISM] = mgr.ESP
        elif result_dict[cn.ALLOCATION_PROGRAM] == mgr.ETKAS:
            if result_dict[cn.EXT_ALLOC_TRIGGERED] == 1:
                if result_dict[cn.ACCEPTANCE_REASON] == cn.T3:
                    result_dict[cn.ALLOCATION_MECHANISM] = (
                        f'{mgr.ETKAS}_RESCUE'
                    )
                else:
                    result_dict[cn.ALLOCATION_MECHANISM] = f'{mgr.ETKAS}_EA'
            else:
                if result_dict[cn.URGENCY_CODE] == mgr.HU:
                    result_dict[cn.ALLOCATION_MECHANISM] = f'{mgr.ETKAS}_HU'
                else:
                    result_dict[cn.ALLOCATION_MECHANISM] = (
                        f'{mgr.ETKAS}_REGULAR'
                    )
        else:
            print(
                f'Warning: unknown allocation mechanism? '
                f'({result_dict[cn.ALLOCATION_PROGRAM]} program)'
            )

        self.transplantations.append(
            result_dict
        )

    def save_posttransplant(
        self,
        date_failure: Optional[datetime],
        date_relist: Optional[datetime],
        cens_date: datetime,
        rereg_id: Optional[int],
        matchr_: 'MatchList.MatchRecord'
    ) -> None:
        """Save relevant post-transplant information"""

        result_dict = {
            cn.ID_DONOR: matchr_.donor.__dict__[cn.ID_DONOR],
            cn.ID_RECIPIENT: matchr_.patient.__dict__[cn.ID_RECIPIENT],
            cn.ID_REREGISTRATION: rereg_id,
            cn.PATIENT_FAILURE_DATE: (
                date_failure if (
                    date_failure and
                    date_failure < cens_date and
                    date_relist is None
                ) else None
            ),
            cn.PATIENT_RELISTING_DATE: (
                date_relist if date_relist and date_relist < cens_date
                else None
            ),
            cn.PATIENT_FAILURE_DATE_RAW: date_failure,
            cn.PATIENT_RELISTING_DATE_RAW: date_relist,
            cn.PATIENT_RELISTING: (
                1 if date_relist and date_relist < cens_date
                else 0 if date_relist
                else 0
            )
        }

        self.posttransplant.append(result_dict)

    def save_exit(self, pat: 'entities.Patient') -> None:
        """Save an exit into the exiting list."""
        result_dict = {
            k: v for k, v in pat.__dict__.items()
            if k in self.cols_to_save_exit
        }

        self.exits.append(
            result_dict
        )

    def save_discard(
            self,
            donor: 'entities.Donor',
            matchl_: List['MatchList.MatchRecord']
    ) -> None:
        """Save an exit into the exiting list."""

        # Count number of center offers & any obl.
        profile_turndowns = 0
        for offer in matchl_:
            if (
                offer.patient.profile is None or
                not offer.patient.profile._check_acceptable(
                    offer.donor
                )
            ):
                profile_turndowns += 1

        if matchl_:
            result_dict = reduce(
                lambda a, b: {**a, **b},
                [o.__dict__ for o in [matchl_[0], donor] if o is not None]
            )
        else:
            result_dict = reduce(
                lambda a, b: {**a, **b},
                [o.__dict__ for o in [donor] if o is not None]
            )
        dict_to_return = {
            k: v for k, v in result_dict.items()
            if k in self.cols_to_save_discard
        }

        # Add relevant discard information.
        dict_to_return[cn.N_OFFERS] = len(matchl_)
        dict_to_return[cn.N_PROFILE_TURNDOWNS] = profile_turndowns

        self.discards.append(
            dict_to_return
        )

    def save_match_list(
        self,
        matchl_: 'MatchList.MatchList',
        write_every_n: int = int(1e6)
    ) -> None:
        """Save match list information"""
        self.match_lists.extend(
            matchl_.return_match_info()
        )
        if len(self.match_lists) >= write_every_n:
            self.match_lists_to_file()
            self.match_lists.clear()

    def return_exits(self) -> pd.DataFrame:
        """Return exits (RM) as a pd.DataFrame"""
        data_ = pd.DataFrame.from_records(
            self.exits,
            columns=self.cols_to_save_exit
        )
        cols = [
            col for col in self.cols_to_save_exit if col in data_.columns
        ]
        data_ = data_.loc[:, cols]
        data_.dropna(how='all', axis='columns', inplace=True)
        return data_

    def return_discards(self) -> pd.DataFrame:
        """Return discards as a pd.DataFrame"""
        data_ = pd.DataFrame.from_records(
            self.discards,
            columns=self.cols_to_save_discard
        )
        cols = [
            col for col in self.cols_to_save_discard if col in data_.columns
        ]
        data_ = data_.loc[:, cols]
        data_.dropna(how='all', axis='columns', inplace=True)
        return data_

    def return_patient_info(
        self,
        patients: Dict[int, 'entities.Patient'],
        selected_cols: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """Return patient info"""

        if selected_cols is None:
            selected_cols = self.cols_to_save_patients

        data_ = pd.DataFrame.from_records(
            list(
                p.__dict__ |
                {
                    cn.ET_MMP: p.get_et_mmp(),
                    cn.VPRA: p.vpra
                } for p in patients.values()
            ),
            columns=selected_cols
        )
        cols = [
            col for col in selected_cols if col in data_.columns
        ]
        data_.dropna(how='all', axis='columns', inplace=True)
        data_ = data_.loc[:, cols]
        data_.loc[:, cn.R_AGE_LISTING] = (
            data_.loc[:, cn.TIME_REGISTRATION] - data_.loc[:, cn.R_DOB]
        ).apply(lambda dt: round_to_decimals(dt / timedelta(days=365), 1))
        return data_

    def save_waitlist_snapshot(
        self,
        patients,
        time: float
    ) -> None:
        data_ = self.return_patient_info(
            patients=patients,
            selected_cols=self.cols_to_save_snapshot
        )
        data_.loc[:, 'snapshot_date'] = (
            timedelta(days=time) + self.sim_start_date
        )
        self.snapshots_waitlist.append(data_)

    def return_all_matchlists(self) -> pd.DataFrame:
        """Return all match lists"""
        data_ = pd.DataFrame.from_records(
            self.match_lists,
            columns=es.MATCH_INFO_COLS
        )
        data_ = data_.loc[:, list(es.MATCH_INFO_COLS)]
        id_cols = [c for c in data_.columns if c.startswith('id')]
        data_ = data_.astype(
            {id_col: 'Int64' for id_col in id_cols}
        )
        data_ = data_.dropna(
            subset=data_.columns.difference(['cn.ACCEPTANCE_REASON']),
            how='all'
        )
        try:
            data_ = data_.astype({cn.ACCEPTED: bool})
            data_.loc[:, cn.ACCEPTED] = (
                data_.loc[:, cn.ACCEPTANCE_REASON].isin(
                    es.ACCEPTANCE_CODES
                )
            )
        except Exception as e:
            print('No acceptance reason column in the following dataset:')
            print(data_)
            data_.loc[:, cn.ACCEPTANCE_REASON] = cn.T3
            data_.loc[:, cn.ACCEPTED] = (
                data_.loc[:, cn.ACCEPTANCE_REASON].isin(
                    es.ACCEPTANCE_CODES
                )
            )

        return data_

    def return_transplantations(
            self,
            patients: Dict[int, 'entities.Patient'],
            cens_date: datetime,
            save_posttxp: bool = True
    ) -> pd.DataFrame:
        """Return transplantation information as DataFrame."""
        data_ = pd.DataFrame.from_records(
            self.transplantations,
            columns=self.cols_to_save_txp + es.OUTPUT_COLS_EXIT_CONSTRUCTED
        )
        cols = [
            col for col in data_.columns
            if col in self.cols_to_save_txp + es.OUTPUT_COLS_EXIT_CONSTRUCTED
        ]
        data_ = data_.loc[:, cols]
        data_.dropna(how='all', axis='columns', inplace=True)
        data_ = data_.astype({cn.ID_REGISTRATION: object})

        if save_posttxp and self.posttransplant:
            # Add post-transplant information
            # (graft failure & patient failure dates)
            d_posttxp = pd.DataFrame.from_records(
                self.posttransplant
            )
            d_posttxp = d_posttxp.astype(
                {cn.ID_REREGISTRATION: object}
            )

            data_ = data_.merge(
                d_posttxp,
                on=[cn.ID_DONOR, cn.ID_RECIPIENT],
                suffixes=(None, '_posttxp'),
                how='left'
            )

            # Add outcomes for re-registrations (whether retransplanted or not)
            d_post_retx = data_.loc[
                :, [cn.ID_REGISTRATION, cn.MATCH_DATE]
            ].dropna().rename(
                columns={
                    cn.ID_REGISTRATION: cn.ID_REREGISTRATION,
                    cn.MATCH_DATE: cn.DATE_RETRANSPLANT
                }
            )
            d_post_retx = d_post_retx.astype(
                {cn.ID_REREGISTRATION: object}
            )
            data_ = data_.merge(
                d_post_retx,
                on=[cn.ID_REREGISTRATION],
                how='left'
            )

            # Add exit statuses for non-retransplanted patients
            d_sim_exits = self.return_exits().loc[
                :, [cn.ID_REGISTRATION, cn.EXIT_DATE, cn.EXIT_STATUS]
            ].rename(
                columns={cn.ID_REGISTRATION: cn.ID_REREGISTRATION}
            )
            d_sim_exits = d_sim_exits.astype({cn.ID_REREGISTRATION: object})

            data_ = data_.merge(
                d_sim_exits,
                on=[cn.ID_REREGISTRATION],
                suffixes=(None, '_rereg'),
                how='left'
            )

            # Add final urgency code for relisted patients.
            data_.loc[
                data_.loc[:, cn.ID_REREGISTRATION].notna() &
                data_.loc[:, cn.DATE_RETRANSPLANT].isna() &
                data_.loc[:, cn.EXIT_STATUS_REREG].isna(),
                cn.EXIT_STATUS_REREG
            ] = data_.loc[
                data_.loc[:, cn.ID_REREGISTRATION].notna() &
                data_.loc[:, cn.DATE_RETRANSPLANT].isna() &
                data_.loc[:, cn.EXIT_STATUS_REREG].isna(),
                cn.ID_REREGISTRATION
            ].apply(
                lambda id: patients[id].urgency_code
            )

            # Construct time to cens
            data_.loc[
                :, 'cens_date'
            ] = (
                data_.loc[
                    :, [
                        cn.PATIENT_FAILURE_DATE,
                        cn.DATE_RETRANSPLANT,
                        'exit_date_rereg'
                    ]
                ].min(axis=1, skipna=True).
                astype('datetime64[ns]').fillna(cens_date)
            )

            data_.loc[
                :, cn.TIME_TO_CENS
            ] = (
                data_.loc[:, 'cens_date'] -
                data_.loc[:, cn.MATCH_DATE]
            ).apply(
                lambda td: td / timedelta(days=1)
            )

            # Construct time-to-patient failure columns
            data_.loc[:, cn.PATIENT_FAILURE] = (
                data_.loc[:, cn.PATIENT_FAILURE_DATE].notna() |
                data_.loc[:, cn.EXIT_STATUS_REREG].isin(['D'])
            ).astype(int)
            data_.loc[
                data_.loc[:, cn.EXIT_STATUS_REREG] == 'D',
                cn.PATIENT_FAILURE_DATE
            ] = data_.loc[
                :, 'exit_date_rereg'
            ]

            data_.loc[:, cn.TIME_TO_PATIENT_FAILURE] = (
                pd.to_datetime(data_.loc[:, cn.PATIENT_FAILURE_DATE]) -
                pd.to_datetime(data_.loc[:, cn.MATCH_DATE])
            ).dt.days

            # Construct time-to-retransplant and time-to-rereg columns
            data_.loc[:, cn.TIME_TO_RETX] = (
                pd.to_datetime(data_.loc[:, cn.DATE_RETRANSPLANT]) -
                pd.to_datetime(data_.loc[:, cn.MATCH_DATE])
            ).dt.days
            data_.loc[:, cn.TIME_TO_REREG] = (
                pd.to_datetime(data_.loc[:, cn.PATIENT_RELISTING_DATE]) -
                pd.to_datetime(data_.loc[:, cn.MATCH_DATE])
            ).dt.days

            # Construct retransplant indicator column
            data_.loc[:, cn.RETRANSPLANTED] = (
                data_.loc[:, cn.DATE_RETRANSPLANT].notna()
            ).astype(int)

            # Construct time-to-event column
            data_.loc[:, cn.EVENT] = (
                data_.loc[:, [cn.RETRANSPLANTED, cn.PATIENT_FAILURE]]
            ).max(axis=1).astype(int)
            data_.loc[:, cn.TIME_TO_EVENT] = (
                data_.loc[:, [cn.TIME_TO_RETX, cn.TIME_TO_PATIENT_FAILURE]]
            ).min(axis=1, skipna=True)

            # Fill unknown with censoring time (administrative only)
            data_.loc[:, cn.TIME_TO_EVENT] = np.ceil(
                data_.loc[:, cn.TIME_TO_EVENT].fillna(
                    data_.loc[:, cn.TIME_TO_CENS]
                ) * 1000
            ) / 1000
            data_.loc[:, cn.TIME_TO_REREG] = np.ceil(
                data_.loc[:, cn.TIME_TO_REREG].fillna(
                    data_.loc[:, cn.TIME_TO_CENS]
                ) * 1000
            ) / 1000
            data_.loc[:, cn.TIME_TO_PATIENT_FAILURE] = np.ceil(
                data_.loc[:, cn.TIME_TO_PATIENT_FAILURE].fillna(
                    data_.loc[:, cn.TIME_TO_CENS]
                ) * 1000
            ) / 1000

        return data_

    def transplantations_to_file(
            self,
            patients: Dict[int, 'entities.Patient'],
            cens_date: datetime,
            save_posttxp: bool = True,
            file: Optional[str] = None
    ) -> None:
        """Save transplantations to file"""
        if file is None:
            file = self.path_txp
        data_ = self.return_transplantations(
            patients=patients,
            cens_date=cens_date,
            save_posttxp=save_posttxp
        )
        outdir = os.path.dirname(file)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        data_.to_csv(
            path_or_buf=file,
            index=False,
            date_format=DEFAULT_DATE_TIME_HMS
        )

    def exits_to_file(self, file: Optional[str] = None):
        """Save exits to file"""
        if file is None:
            file = self.path_exits
        data_ = self.return_exits()
        data_.rename(
            columns={cn.URGENCY_REASON: cn.REMOVAL_REASON},
            inplace=True
        )
        outdir = os.path.dirname(file)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        data_.to_csv(
            path_or_buf=file,
            index=False,
            date_format=DEFAULT_DATE_TIME_HMS
        )

    def discards_to_file(self, file: Optional[str] = None) -> None:
        """Save patient deaths to file"""
        if file is None:
            file = self.path_dsc
        data_ = self.return_discards()
        outdir = os.path.dirname(file)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        data_.to_csv(
            path_or_buf=file,
            index=False,
            date_format=DEFAULT_DATE_TIME_HMS
        )

    def patients_to_file(
        self,
        patients: Dict[int, 'entities.Patient'],
        file: Optional[str] = None
    ) -> None:
        """Save patient deaths to file"""
        if file is None:
            file = self.path_pat
        data_ = self.return_patient_info(
            patients=patients
        )
        outdir = os.path.dirname(file)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        data_.to_csv(
            path_or_buf=file,
            index=False,
            date_format=DEFAULT_DATE_TIME_HMS
        )

    def waitlist_snapshots_to_file(
        self,
        file: Optional[str] = None
    ) -> None:
        if file is None:
            file = self.path_snapshot
        if file and self.snapshots_waitlist:

            d_snapshots = pd.concat(
                self.snapshots_waitlist
            )
            d_snapshots.to_csv(
                path_or_buf=file,
                index=False,
                date_format=DEFAULT_DATE_TIME_HMS
            )

    def match_lists_to_file(
            self,
            file: Optional[str] = None,
            compression: Optional[str] = None,
            write_all: Optional[bool] = False
    ) -> None:
        """Save match lists to file"""
        file = file if file else self.path_ml
        if compression:
            if isinstance(file, str) and not file.endswith(compression):
                file += f'.{compression}'

        data_ = self.return_all_matchlists()
        if not write_all:
            data_ = data_.loc[
                (data_.loc[:, cn.ACCEPTANCE_REASON] != 'FP') &
                data_.loc[:, cn.ACCEPTANCE_REASON].notna(),
                :
            ]

        outdir = os.path.dirname(file)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        csvfile = pathlib.Path(file)
        data_.to_csv(
            path_or_buf=file,
            mode='a',
            compression=compression,
            index=False,
            date_format=DEFAULT_DATE_TIME_HMS,
            header=not csvfile.exists()
        )

    def gzip_csv_file(
            self,
            file: str
    ) -> None:
        """Save match lists to file"""
        with open(file, 'rb') as f_in:
            with gzip.open(f'{file}.gzip', 'wb', compresslevel=5) as f_out:
                shutil.copyfileobj(f_in, f_out, length=128 * 1024)
        os.remove(file)

    def save_outcomes_to_file(
        self,
        patients: Dict[int, 'entities.Patient'],
        cens_date: datetime,
        save_posttxp: bool = True
    ) -> None:
        """Save all outcomes to file"""
        self.discards_to_file()
        self.transplantations_to_file(
            patients=patients,
            cens_date=cens_date,
            save_posttxp=save_posttxp
        )
        self.exits_to_file()
        self.waitlist_snapshots_to_file()
        if patients:
            self.patients_to_file(
                patients
            )
        if self.save_ml:
            self.match_lists_to_file()
            self.gzip_csv_file(
                file=self.path_ml
            )
