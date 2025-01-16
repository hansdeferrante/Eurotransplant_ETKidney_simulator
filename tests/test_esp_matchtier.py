#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri 11-02-2022

@author: H.C. de Ferrante
"""

import sys
sys.path.append('./')

from datetime import timedelta
import unittest
from simulator.code.entities import Patient, Donor
import simulator.magic_values.column_names as cn
import simulator.code.utils.read_input_files as rdr
import simulator.magic_values.etkidney_simulator_settings as es
from simulator.code.HLA.HLASystem import HLASystem
import pandas as pd
import os
from collections import defaultdict
from simulator.code.matchlist.MatchListETKASandESP import MatchListESP
from simulator.code.utils.load_simulation_entities import (
    load_balances
)


class TestESPMatchlist(unittest.TestCase):
    """Test whether match qualities are correctly
        constructed for the ESP program.
    """

    def test_esp_match_lists(
            self,
            num_records_to_test=10000,
            max_fraction_wrong=.01,
            verbose=False
    ):
        """
        Unit test to check whether tiers are correctly constructed
        for ESP match lists.
        """

        pd.set_option('display.max_rows', 50)
        ss = rdr.read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )

        # Load HLA system and match lists
        ss.needed_broad_mismatches = ('hla_b', 'hla_a')
        ss.needed_split_mismatches = ('hla_dr',)
        hla_system = HLASystem(ss)
        bal_system = load_balances(ss)

        df = rdr._read_with_datetime_cols(
            input_path='data/test/recent_esp_match_list.csv',
            dtps={
                'ID_MTR': int,
                'ID_DONOR': int,
                'PATIENT_MATCH_RANK_LAYER': str,
                'DONOR_CENTER': str,
                "PATIENT_CENTER": str,
                'D_HLA_FULL': str,
                "PATIENT_HLA": str,
                "MATCH_QUALITY": str,
                "WFMR_R_AGE": int,
                "DONOR_COUNTRY": str,
                "PATIENT_COUNTRY": str,
                "PATIENT_OFFERED": int,
                "WFTS_MATCH_DATE": object,
                "DONOR_REGION": str,
                "PATIENT_REGION": str,
                "WL_DATE_FIRST_DIALYSIS": object,
                "WLKI_PROGRAMME": str,
                "D_CENTER_SUBREGION_ESP": str
            },
            casecols=True,
            datecols=["WL_DATE_FIRST_DIALYSIS", "WFTS_MATCH_DATE"]
        )

        df = df.loc[df.patient_hla.str.len() > 0, :]
        df = df.loc[df.d_hla_full.str.len() > 0, :]
        df = df.loc[df.match_quality != 'MAN', :]
        df.loc[:, 'match_tier'] = (
            df.loc[:, 'patient_match_rank_layer'].str.slice(0, 1)
        )

        rcrd = {}
        patients = [
            Patient(
                id_recipient=1,
                r_dob=(
                    rcrd['wfts_match_date'] -
                    timedelta(days=es.DAYS_PER_YEAR_FLOAT * rcrd['wfmr_r_age'])
                ),
                recipient_country=rcrd['patient_country'],
                recipient_region=rcrd['patient_region'],
                recipient_center=rcrd['patient_center'],
                bloodgroup='O',
                hla=pat_hla,
                listing_date=(
                    rcrd['wfts_match_date'] - 10 * es.DAYS_PER_YEAR
                ),
                urgency_code='T',
                sim_set=ss,
                hla_system=hla_system,
                date_first_dial=rcrd['wl_date_first_dialysis'],
                kidney_program=rcrd['wlki_programme']
            )
            for pat_hla, rcrd in zip(
                rdr.fix_hla_string(df.patient_hla),
                df.to_dict(orient='records')
            )
        ]

        donors = [
            Donor(
                id_donor=1255,
                donor_country=rcrd[cn.DONOR_COUNTRY],
                donor_region=rcrd[cn.DONOR_REGION],
                donor_center=rcrd[cn.DONOR_CENTER],
                donor_subregion_esp=rcrd[cn.D_ESP_SUBREGION],
                bloodgroup='O',
                reporting_date=rcrd['wfts_match_date'],
                donor_dcd=False,
                weight=80,
                hla=don_hla,
                hla_system=hla_system,
                age=40,
                death_cause_group='CVA',
                malignancy=False,
                tumor_history=False,
                donor_marginal_free_text=False,
                drug_abuse=False,
                n_kidneys_available=2,
                diabetes=False,
                cardiac_arrest=False,
                last_creat=0.5,
                urine_protein=0,
                smoker=rcrd.get(cn.D_SMOKING, False),
                cmv=rcrd.get(cn.D_CMV, False),
                hbsag=rcrd.get(cn.D_HBSAG, False),
                hcvab=rcrd.get(cn.D_HCVAB, False),
                hbcab=rcrd.get(cn.D_HBCAB, False),
                sepsis=rcrd.get(cn.D_SEPSIS, False),
                meningitis=rcrd.get(cn.D_MENINGITIS, False),
                hypertension=False,
                euthanasia=rcrd.get(cn.D_EUTHANASIA, False),
                rescue=False,
                sim_set=ss
            ) for don_hla, rcrd in zip(
                rdr.fix_hla_string(df.d_hla_full),
                df.to_dict(orient='records')
            )
        ]

        errors = defaultdict(int)

        for d, p, mq, _, _, _, _, _ in zip(
            donors, patients, df.match_tier,
            df.wlki_programme, df.id_mtr, df.patient_offered,
            df.wfts_match_date, df.d_center_subregion_esp
        ):
            try:
                match_list_esp = MatchListESP(
                    patients=(p,),
                    donor=d,
                    match_date=d.reporting_date.to_pydatetime(),
                    sim_start_date=d.sim_set.SIM_START_DATE,
                    hla_system=hla_system,
                    bal_system=bal_system,
                    calc_points=ss.calc_esp_score,
                    store_score_components=(ss.STORE_SCORE_COMPONENTS)
                )
            except Exception as e:
                print('Failed to make ESP match list for this record')
                print(p.__dict__)
                print(e)
                exit()

            mr = match_list_esp.match_list[0]

            if mr.__dict__[cn.MTCH_TIER] != mq:
                errors[mr.donor_alloc_center] += 1

        mean_correct = round(sum(errors.values()) / len(donors), 4)
        print(
            f'{mean_correct * 100}% records have incorrect '
            f'match tier in ESP'
        )
        print(errors)
        assert mean_correct <= max_fraction_wrong


if __name__ == '__main__':
    unittest.main()
