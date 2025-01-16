import sys
import unittest
import os
from datetime import datetime
import pandas as pd

sys.path.append('./')

from simulator.code.AcceptanceModule import AcceptanceModule
from simulator.code.entities import Donor
from simulator.code.HLA.HLASystem import HLASystem
import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es
from simulator.code.utils.read_input_files import \
    read_sim_settings, read_travel_times


class TestCurrentAllocation(unittest.TestCase):
    """Test whether current allocation is done correctly.
    """

    def test_rescue_trigger(self):
        """Test center acceptances
        """

        print(
            'Testing whether acceptance model is'
            ' correct for center/obligations.'
        )

        # Read in simulation settings
        ss = read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )
        hla_system = HLASystem(sim_set=ss)
        travel_time_dict = read_travel_times()

        d_test_donors = pd.read_csv('data/test/triggered_rescue.csv')
        acc_module = AcceptanceModule(
            seed=1,
            patient_acc_policy='LR',
            center_acc_policy='LR',
            verbose=0,
            simulate_rescue=True,
            simulate_random_effects=False
        )

        dummy_date = datetime(year=2000, month=1, day=1)

        dummy_hla = 'A19 A28 B12 B13 DR7 DR53'

        n_bg_incompatible = 0
        for id, rcrd in d_test_donors.iterrows():

            don = Donor(
                id_donor=1255,
                donor_country=rcrd[cn.DONOR_COUNTRY],
                donor_region='Bel_1',
                donor_center=rcrd[cn.DONOR_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                reporting_date=dummy_date,
                weight=rcrd[cn.DONOR_WEIGHT],
                height=rcrd[cn.DONOR_HEIGHT],
                donor_dcd=rcrd[cn.D_DCD],
                age=rcrd['donor_age'],
                death_cause_group=rcrd[cn.DEATH_CAUSE_GROUP],
                malignancy=rcrd[cn.D_MALIGNANCY],
                tumor_history=False,
                donor_marginal_free_text=False,
                drug_abuse=rcrd[cn.D_DRUG_ABUSE],
                n_kidneys_available=2,
                hla=dummy_hla,
                hla_system=hla_system,
                diabetes=rcrd[cn.D_DIABETES],
                cardiac_arrest=rcrd[cn.D_CARREST],
                last_creat=rcrd[cn.D_LAST_CREAT],
                urine_protein=rcrd[cn.D_URINE_PROTEIN],
                smoker=rcrd.get(cn.D_SMOKING, False),
                cmv=rcrd.get(cn.D_CMV, False),
                hbsag=rcrd.get(cn.D_HBSAG, False),
                hcvab=rcrd.get(cn.D_HCVAB, False),
                hbcab=rcrd.get(cn.D_HBCAB, False),
                sepsis=rcrd.get(cn.D_SEPSIS, False),
                meningitis=rcrd.get(cn.D_MENINGITIS),
                hypertension=rcrd[cn.D_HYPERTENSION],
                euthanasia=rcrd.get(cn.D_EUTHANASIA, False),
                rescue=False,
                sim_set=ss,
                donor_subregion_esp=rcrd[cn.DONOR_CENTER]
            )

            pd.set_option('display.max_rows', None)
            phat = acc_module.predict_rescue_prob(
                don,
                kth_offer=10,
                verbose=0
            )

            diff_prob = (
                phat - (1 - rcrd.phat_10)
            )
            assert abs(diff_prob) <= 1e-2, \
                (
                    f'Turndown probability incorrect for test case '
                    f'{id} ({diff_prob})'
            )

            n_offers = acc_module.simulate_offers_until_nonstandard_alloc(
                donor=don,
                r_prob=phat
            )
            assert n_offers <= 10, \
                (
                    f"At most 10 offer should've been "
                    f"generated for {id} (not {n_offers})"
                )


if __name__ == '__main__':
    unittest.main()
