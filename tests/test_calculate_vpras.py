#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri 11-02-2022


@author: H.C. de Ferrante
"""

import pandas as pd
import os
from numpy import random
import unittest
import sys
sys.path.append('./')
if True:  # noqa E402
    import simulator.code.utils.read_input_files as rdr
    import simulator.magic_values.etkidney_simulator_settings as es
    from simulator.code.HLA.HLASystem import HLASystem


class TestCalculateVPRA(unittest.TestCase):
    """Test whether VPRA is correctly constructed
    """

    def test_vpra(
            self,
            num_records_to_test=20000,
            min_fraction_correct=.995,
            verbose=False
    ):
        """
            Unit test to test whether VPRA is correctly
            calculated.
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

        # Read in all transplantations
        df = pd.read_csv(
            '../raw_data/transplantations.csv',
            encoding='utf-8',
            low_memory=False
        )

        df = df.loc[
            :,
            [
                'ID_REC', 'CURRENT_PATIENT_HLA',
                'VPRA', 'UNACCEPTABLE_ANTIGENS'
            ]
        ]
        df = df.loc[df.UNACCEPTABLE_ANTIGENS.str.len() > 0, :]
        df = df.loc[~ df.VPRA.isna(), :]

        random.seed(1)
        df = df.sample(500)

        calculated_vpras = (
            {
                'unacc': unacc,
                'calculated_vpra': hla_system.calculate_vpra_from_string(
                    pat_hla
                ),
                'database_vpra': pat_vpra
            } for unacc, pat_hla, pat_vpra in zip(
                df.UNACCEPTABLE_ANTIGENS,
                rdr.fix_hla_string(df.UNACCEPTABLE_ANTIGENS),
                df.VPRA
            )
        )
        df_vpras = pd.DataFrame.from_records(calculated_vpras)
        df_vpras.loc[:, 'vpra_difference'] = (
            df_vpras.calculated_vpra * 100 -
            df_vpras.database_vpra
        )

        med_diff = df_vpras.vpra_difference.median()
        print(
            f'Median difference in calculated vs. '
            f'database vPRA: {med_diff}%'
        )
        assert (abs(med_diff)) < 1, (
            f'Calculated vPRA appears to be miscalibrated.'
            f'Median difference is {med_diff}'
        )


if __name__ == '__main__':
    unittest.main()
