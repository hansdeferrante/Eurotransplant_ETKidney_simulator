#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri 11-02-2022


@author: H.C. de Ferrante
"""

import sys
sys.path.append('./')
if True:  # noqa E402
    from datetime import timedelta
    import unittest
    from simulator.code.entities import Patient, Donor
    import simulator.magic_values.column_names as cn
    import simulator.code.utils.read_input_files as rdr
    import simulator.magic_values.etkidney_simulator_settings as es
    from simulator.code.HLA.HLASystem import HLASystem
    import pandas as pd
    import os
    import simulator.magic_values.magic_values_rules as mr
    import random
    from statistics import mean


class TestConstructMatchQuality(unittest.TestCase):
    """Test whether match qualities are correctly constructed
    """

    def test_etkas_match_lists(
            self,
            num_records_to_test=10000,
            min_fraction_correct=.995,
            verbose=False
    ):
        """
            Unit test to test whether mismatches are correctly constructed
            for ETKAS match lists. There may be some errors here, because
            the patient HLA at the start of the month was taken (instead of the
            last known HLA). Match qualities manually checked (verbose=True)
            shows that match codes are actually correctly constructed.
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
        DUMMY_DATE = pd.Timestamp('2000-01-01')

        # Read in random offers
        random.seed(14)
        n = 882053
        min_n = round(n / 2)
        records_to_skip = (
            sorted(range(1, min_n)) +
            sorted(
                random.sample(
                    range(min_n, n),
                    (n - min_n) - num_records_to_test)
            )
        )

        df = pd.read_csv(
            '../raw_data/processed/all_etkas_match_lists_with_mq_vpra.csv',
            encoding='utf-8',
            skiprows=records_to_skip,
            low_memory=False
        )
        df = df.loc[
            :,
            [
                'D_HLA_FULL', 'PATIENT_HLA', 'MATCH_QUALITY',
                'ID_MTR', 'PATIENT_OFFERED', 'WFTS_MATCH_DATE'
            ]
        ]

        df = df.loc[df.PATIENT_HLA.str.len() > 0, :]
        df = df.loc[df.D_HLA_FULL.str.len() > 0, :]
        df = df.loc[df.MATCH_QUALITY != 'MAN', :]

        rcrd = {}
        patients = [
            Patient(
                id_recipient=1,
                r_dob=DUMMY_DATE - 45 * es.DAYS_PER_YEAR,
                recipient_country='Belgium',
                recipient_region='region',
                recipient_center='BLMTP',
                bloodgroup='O',
                hla=pat_hla,
                listing_date=DUMMY_DATE,
                urgency_code='T',
                sim_set=ss,
                hla_system=hla_system
            )
            for pat_hla in rdr.fix_hla_string(df.PATIENT_HLA)
        ]

        donors = [
            Donor(
                id_donor=1255,
                donor_country='Belgium',
                donor_region=None,
                donor_center='BLMTP',
                bloodgroup='O',
                reporting_date=DUMMY_DATE,
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
                sim_set=ss,
                donor_subregion_esp='Belgium'
            ) for don_hla in rdr.fix_hla_string(df.D_HLA_FULL)
        ]

        match_qualities = list()
        correct_match = list()

        for d, p, mq, id_mtr, id_pat, date in zip(
            donors, patients, df.MATCH_QUALITY, df.ID_MTR, df.PATIENT_OFFERED,
            df.WFTS_MATCH_DATE
        ):
            if p.hla is None or d.hla is None:
                continue

            mm_dict = hla_system.count_mismatches(
                d_hla=d.hla,
                p_hla=p.hla
            )

            mm_str = "".join(str(mm_dict[k]) for k in (
                mr.MMB_HLA_A, mr.MMB_HLA_B, mr.MMS_HLA_DR
            )
            ) if mm_dict else None

            if mm_str and mq != mm_str:
                correct_match.append(0)
                if verbose:
                    print(id_mtr)
                    print(id_pat)
                    print(date)

                    print(f'Original: {mq} != Constructed: {mm_str}')
                    print(f'Donor broads: {d.hla.broads}')
                    print(f'Donor splits: {d.hla.splits}')
                    print(f'Patient broads: {p.hla.broads}')
                    print(f'Patient splits: {p.hla.splits}')
                    print(p.hla)
            else:
                correct_match.append(1)

            match_qualities.append(mm_str)

        mean_correct = round(mean(correct_match), 4)
        print(
            f'{mean_correct * 100}% records have correct '
            f'match quality for ETKAS match records'
        )
        assert mean_correct > min_fraction_correct

    def test_transplant_quality(
            self,
            num_records_to_test=10000,
            min_fraction_correct=.995,
            verbose=False
    ):
        """Unit test to test whether mismatches are correctly constructed
            for ETKAS match lists. There may be some errors here, because
            the patient HLA at the start of the month was taken (instead of the
            last known HLA). Match qualities manually checked (verbose=True)
            shows that match codes are actually correctly constructed.
        """

        pd.set_option('display.max_rows', 50)
        ss = rdr.read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )

        # Load HLA system and match lists
        ss.needed_broad_mismatches = ('hla_b', 'hla_a',)
        ss.needed_split_mismatches = ('hla_dr',)
        hla_system = HLASystem(ss)
        DUMMY_DATE = pd.Timestamp('2000-01-01')

        # Read in random offers
        random.seed(14)
        n = 44839
        min_n = round(n / 2)
        records_to_skip = (
            sorted(range(1, min_n)) +
            sorted(
                random.sample(
                    range(min_n, n),
                    (n - min_n) - num_records_to_test)
            )
        )
        df = pd.read_csv(
            '../raw_data/transplantations.csv',
            encoding='utf-8',
            skiprows=records_to_skip,
            low_memory=False
        )

        df = df.loc[
            :,
            [
                'DONOR_HLA', 'CURRENT_PATIENT_HLA', 'MATCH_QUALITY',
                'ID_MTR', 'ID_REC', 'DATE_TXP'
            ]
        ]

        df = df.loc[df.CURRENT_PATIENT_HLA.str.len() > 0, :]
        df = df.loc[df.DONOR_HLA.str.len() > 0, :]
        df = df.loc[df.MATCH_QUALITY.fillna('MAN') != 'MAN', :]
        df = df.loc[
            df.loc[:, 'CURRENT_PATIENT_HLA'].str.
            contains('A.*B.*DR.*', regex=True, na=True)
        ]
        df = df.loc[
            df.loc[:, 'DONOR_HLA'].str.
            contains('A.*B.*DR.*', regex=True, na=True)
        ]

        rcrd = {}
        patients = [
            Patient(
                id_recipient=1,
                r_dob=DUMMY_DATE - 45 * es.DAYS_PER_YEAR,
                recipient_country='Belgium',
                recipient_region='region',
                recipient_center='BLMTP',
                bloodgroup='O',
                hla=pat_hla,
                listing_date=DUMMY_DATE,
                urgency_code='T',
                sim_set=ss,
                hla_system=hla_system
            )
            for pat_hla in rdr.fix_hla_string(df.CURRENT_PATIENT_HLA)
        ]

        donors = [
            Donor(
                id_donor=1255,
                donor_country='Belgium',
                donor_region=None,
                donor_center='BLMTP',
                bloodgroup='O',
                reporting_date=DUMMY_DATE,
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
                sim_set=ss,
                donor_subregion_esp='Belgium'
            )
            for don_hla in rdr.fix_hla_string(df.DONOR_HLA)
        ]

        match_qualities = list()
        correct_match = list()

        for d, p, mq, id_mtr, id_pat, date in zip(
            donors, patients, df.MATCH_QUALITY, df.ID_MTR, df.ID_REC,
            df.DATE_TXP
        ):
            if p.hla is None or d.hla is None:
                continue

            mm_dict = hla_system.count_mismatches(
                d_hla=d.hla,
                p_hla=p.hla
            )

            mm_str = "".join(str(mm_dict[k]) for k in (
                mr.MMB_HLA_A, mr.MMB_HLA_B, mr.MMS_HLA_DR
            )
            ) if mm_dict else None

            if mm_str and mq != mm_str:
                correct_match.append(0)
                if verbose:
                    print(id_mtr)
                    print(id_pat)
                    print(date)
                    print(p.hla)
                    print(d.hla)
                    print(f'Original: {mq} != Constructed: {mm_str}')
                    print(f'Donor broads: {d.hla.broads}')
                    print(f'Donor splits: {d.hla.splits}')
                    print(f'Patient broads: {p.hla.broads}')
                    print(f'Patient splits: {p.hla.splits}')
                    print(p.hla)
            else:
                correct_match.append(1)

            match_qualities.append(mm_str)

        mean_correct = round(mean(correct_match), 4)
        print(
            f'{mean_correct * 100}% records correct match quality '
            f'for transplantations'
        )
        assert mean_correct > min_fraction_correct


if __name__ == '__main__':
    unittest.main()
