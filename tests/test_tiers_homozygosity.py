import sys
import unittest
import pandas as pd
import os
from itertools import groupby

sys.path.append('./')

from datetime import timedelta
from simulator.code.utils.load_simulation_entities import load_balances
from simulator.code.entities import Patient, Donor
import simulator.magic_values.column_names as cn
from simulator.magic_values.inputfile_settings import DTYPE_OFFERLIST
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.code.utils.read_input_files as rdr
from simulator.code.matchlist.MatchListETKASandESP import \
    MatchListCurrentETKAS
from simulator.code.HLA.HLASystem import HLASystem, Unacceptables
from simulator.code.HLA.MMPSystem import MMPSystem

from simulator.code.utils.read_input_files import \
    read_sim_settings


class TestKidneyMatchListRanks(unittest.TestCase):
    """ Test whether kidney match lists yield the exact
        same order as match lists exported from the ET
        data warehouse
    """

    def test_mrl(self):
        """Test whether kidney match list ranks are exactly reproduced.
        """

        pd.set_option('display.max_rows', 50)
        ss = read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )

        # Read in the real MatchLists
        n = 882053
        records_to_read = [0] + list(range(318351, 318512))
        records_to_skip = (
            x for x in range(n)
            if x not in records_to_read
        )

        data_types = DTYPE_OFFERLIST
        data_types['wfmr_xc_mism'] = 'object'
        data_types['wfmr_xc_bonus_0mm'] = 'object'
        data_types['patient_offered'] = 'str'
        df_ml = rdr._read_with_datetime_cols(
            input_path='data/test/match_list_test.csv',
            dtps=DTYPE_OFFERLIST,
            datecols=[cn.WFTS_MATCH_DATE, cn.PATIENT_DOB],
            casecols=True
        )
        df_ml.columns = map(str.lower, df_ml.columns)
        df_ml = df_ml.loc[
            :, [
                cn.ID_MTR, 'patient_offered', cn.WFTS_MATCH_DATE,
                cn.PATIENT_DOB, cn.MATCH_RANK_LAYER, cn.PATIENT_HLA_LITERAL,
                cn.D_HLA_FULL, 'donor_center', 'patient_center',
                'unacceptable_antigens',
                "wfmr_xc_wait", "wfmr_xc_mism", "wfmr_xc_bonus_0mm",
                "wfmr_xc_mmp", "wfmr_xc_mmp_split", "wfmr_xc_balance",
                "wfmr_xc_balance_reg", "wfmr_xc_dist", "wfmr_xc_bonus_paed",
                "wfmr_xc_bonus_hu", cn.PATIENT_COUNTRY, cn.PATIENT_REGION,
                cn.DONOR_REGION, cn.DONOR_COUNTRY, cn.GRAFT_BLOODGROUP
            ]
        ]
        hla_system = HLASystem(ss)
        mmp_system = MMPSystem(ss, hla_system)
        k = 0

        def set_unacceptables(pat: Patient, unacc: str):
            pat.unacceptable_antigens = Unacceptables(hla_system, unacc)
            return pat

        for _, df_sel in groupby(
            df_ml.to_dict(orient='records'),
            lambda x: (x[cn.ID_MTR])
        ):
            k += 1
            if k % 100 == 0:
                print(f'Finished checking {k} matchlists.')

            patient_list = []

            # Construct patient list
            rcrds = list(df_sel)
            ss.SIM_START_DATE = rcrds[0][cn.WFTS_MATCH_DATE].to_pydatetime()
            bal_system = load_balances(ss)
            CURRENT_DATE = ss.SIM_START_DATE
            patient_list = [
                set_unacceptables(
                    Patient(
                        id_recipient=i,
                        r_dob=rcrd[cn.PATIENT_DOB],
                        recipient_country=rcrd[cn.PATIENT_COUNTRY],
                        recipient_region=rcrd[cn.PATIENT_REGION],
                        recipient_center=rcrd[cn.PATIENT_CENTER],
                        bloodgroup=rcrd[cn.GRAFT_BLOODGROUP],
                        hla=rcrd[cn.PATIENT_HLA_LITERAL],
                        listing_date=CURRENT_DATE,
                        urgency_code='T',
                        sim_set=ss,
                        hla_system=hla_system,
                        mmp_system=mmp_system,
                        date_first_dial=(
                            CURRENT_DATE - timedelta(
                                days=(rcrd['wfmr_xc_wait'] / 33.33 * 365)
                            )
                        )
                    ),
                    rcrd['unacceptable_antigens']
                ) for i, rcrd in enumerate(rcrds)
            ]

            # Construct donor list
            donors = [
                Donor(
                    id_donor=1255,
                    donor_country=rcrd[cn.DONOR_COUNTRY],
                    donor_region=rcrd[cn.DONOR_REGION],
                    donor_center=rcrd[cn.DONOR_CENTER],
                    bloodgroup=rcrd[cn.GRAFT_BLOODGROUP],
                    reporting_date=CURRENT_DATE,
                    donor_dcd=False,
                    weight=80,
                    hla=rcrd[cn.D_HLA_FULL],
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
                    donor_subregion_esp=rcrd[cn.DONOR_CENTER]
                ) for rcrd in rcrds[0:1]
            ]

            ml = MatchListCurrentETKAS(
                patients=patient_list,
                donor=donors[0],
                match_date=CURRENT_DATE,
                hla_system=hla_system,
                bal_system=bal_system,
                calc_points=ss.calc_etkas_score,
                sim_start_date=ss.SIM_START_DATE,
                store_score_components=True
            )

            mtch_tiers = [
                mr.mtch_tier_str
                for mr in ml.return_match_list()
            ]
            orig_mtch_tiers = df_ml.loc[:, cn.MATCH_RANK_LAYER].values

            assert all(orig_mtch_tiers == mtch_tiers), (
                f'Match tiers are incorrectly returned'
            )


if __name__ == '__main__':
    unittest.main()
