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
    from simulator.code.HLA.HLASystem import HLASystem, Unacceptables
    import pandas as pd
    import os
    import simulator.magic_values.magic_values_rules as mr
    import random
    from statistics import mean
    from simulator.code.HLA.MMPSystem import MMPSystem


class TestConstructMMP(unittest.TestCase):
    """Test whether match qualities are correctly constructed
    """

    def test_etkas_match_lists(
            self,
            num_records_to_test=1000,
            min_fraction_correct=0.90,
            verbose=False
    ):
        """
            Unit test to test whether mismatch probs are correctly constructed
            for ETKAS match lists. There may be some errors here, because
            the last known patient MMP at the start of the month was taken.
            Sometimes, these are dated because patients updated their vPRA.
        """

        pd.set_option('display.max_rows', 100)
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
        mmp_system = MMPSystem(ss, hla_system=hla_system)
        DUMMY_DATE = pd.Timestamp('2000-01-01')

        # Read in random offers
        random.seed(14)
        n = 882053
        min_n = round(4 / 5 * n)
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
                'PATIENT_OFFERED', 'PATIENT_HLA',
                'WFTS_MATCH_DATE', 'GRAFT_BLOODGROUP',
                'UNACCEPTABLE_ANTIGENS', 'WFMR_XC_MMP_SPLIT'
            ]
        ]

        df = df.loc[df.PATIENT_HLA.str.len() > 0, :]
        df = df.loc[~df.WFMR_XC_MMP_SPLIT.isna(), :]
        df = df.reset_index()

        rcrd = {}

        def set_unacceptables(pat: Patient, unacc: str):
            pat.unacceptable_antigens = Unacceptables(hla_system, unacc)
            return pat
        patients = [
            set_unacceptables(
                Patient(
                    id_recipient=1,
                    r_dob=DUMMY_DATE - 45 * es.DAYS_PER_YEAR,
                    recipient_country='Belgium',
                    recipient_region='region',
                    recipient_center='BLMTP',
                    bloodgroup=pat_bg,
                    hla=pat_hla,
                    listing_date=DUMMY_DATE,
                    urgency_code='T',
                    sim_set=ss,
                    hla_system=hla_system,
                    mmp_system=mmp_system
                ),
                pat_unacc
            )
            for pat_hla, pat_bg, pat_unacc in zip(
                rdr.fix_hla_string(df.PATIENT_HLA),
                df.GRAFT_BLOODGROUP,
                df.UNACCEPTABLE_ANTIGENS
            )
        ]

        df.loc[:, 'constructed_mmp'] = list(
            pat.get_et_mmp()
            for pat in patients
        )
        df.loc[:, 'reduced_hla'] = list(
            pat.hla.reduced_hla
            for pat in patients
        )
        df.loc[:, 'vpra'] = list(
            pat.vpra
            for pat in patients
        )
        df.loc[:, 'diff'] = (
            df.loc[:, 'constructed_mmp'] -
            df.loc[:, 'WFMR_XC_MMP_SPLIT']
        )

        df.loc[:, 'diff_within_10pct'] = df.loc[
            :,
            'diff'
        ].between(-10, 10)

        frac_correct = df.diff_within_10pct.mean()
        assert frac_correct > min_fraction_correct, (
            f'Only {round(frac_correct * 100, 1)}% of '
            f'calculated MMPs are within 10 pct points '
            f'of database-MMP'
        )
        print(
            f'{round(frac_correct * 100, 1)}% of calculated MMPs '
            f'are within 10 pct points of database-MMP'
        )
        if verbose:
            print(
                df.sort_values(
                    by=['diff']
                ).loc[
                    :,
                    [
                        'vpra', 'GRAFT_BLOODGROUP',
                        'reduced_hla',
                        'WFMR_XC_MMP_SPLIT',
                        'constructed_mmp',
                        'UNACCEPTABLE_ANTIGENS',
                        'diff_within_10pct'
                    ]
                ]
            )


if __name__ == '__main__':
    unittest.main()
