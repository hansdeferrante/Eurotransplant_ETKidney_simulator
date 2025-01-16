import sys
import unittest
import os
from datetime import timedelta, datetime
import pandas as pd

sys.path.append('./')

from simulator.code.PostTransplantPredictor import PostTransplantPredictor
from simulator.code.entities import Patient, Donor
from simulator.code.HLA.HLASystem import HLASystem
from simulator.code.HLA.MMPSystem import MMPSystem
import simulator.magic_values.inputfile_settings as ifs
import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es
from simulator.code.utils.read_input_files import \
    read_sim_settings, fix_hla_string
from simulator.code.utils.load_simulation_entities import load_balances
from simulator.code.matchlist.MatchListETKASandESP import \
    MatchListCurrentETKAS


class TestCurrentAllocation(unittest.TestCase):
    """Test whether current allocation is done correctly.
    """

    def verify_posttxp_surv(
            self,
            path_test_cases: str,
            verbosity: int = 2,
            slack: float = 5e-2
    ) -> None:
        """Test whether ML have acceptable donor/recipient combinations.
        """

        # Read in simulation settings
        ss = read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )

        d_txps = pd.read_csv(path_test_cases)
        d_txps['date_transplanted'] = d_txps['date_txp'].apply(
            lambda x: datetime.strptime(
                x,
                ifs.DEFAULT_DATE_TIME_HMS
            )
        ).dt.to_pydatetime()
        d_txps[cn.DATE_FIRST_DIAL] = (
            d_txps['date_transplanted'] -
            d_txps[cn.YEARS_ON_DIAL].map(
                lambda x: timedelta(days=x) * es.DAYS_PER_YEAR_FLOAT
            )
        )
        d_txps.current_patient_hla = fix_hla_string(
            d_txps.current_patient_hla
        )
        hla_system = HLASystem(ss)
        mmp_system = MMPSystem(ss, hla_system)
        bal_system = load_balances(ss, update_balances=False)
        ptp = PostTransplantPredictor(
            offset_ids_transplants=99999,
            seed=14
        )

        n_bg_incompatible = 0
        for id, rcrd in d_txps.iterrows():

            pat = Patient(
                id_recipient=1,
                date_first_dial=rcrd[cn.DATE_FIRST_DIAL],
                r_dob=(
                    rcrd['date_transplanted'] -
                    timedelta(days=rcrd[cn.R_MATCH_AGE] * 365.24)
                ),
                recipient_country=rcrd[cn.PATIENT_COUNTRY],
                recipient_region='Netherlands',
                recipient_center=rcrd[cn.PATIENT_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                listing_date=rcrd['date_transplanted'] - timedelta(days=100),
                urgency_code='T',
                sex=rcrd[cn.R_SEX],
                sim_set=ss,
                time_since_prev_txp=rcrd[cn.TIME_SINCE_PREV_TXP],
                type_retx=(
                    cn.RETRANSPLANT_DURING_SIM if rcrd[cn.IS_RETRANSPLANT]
                    else cn.NO_RETRANSPLANT
                ),
                hla=rcrd[cn.CURRENT_PATIENT_HLA],
                hla_system=hla_system,
                mmp_system=mmp_system
            )

            don = Donor(
                id_donor=1255,
                donor_country=rcrd[cn.DONOR_COUNTRY],
                donor_region='Germany',
                donor_center=rcrd[cn.DONOR_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                reporting_date=rcrd['date_transplanted'],
                weight=80,
                height=180,
                donor_dcd=rcrd[cn.D_DCD],
                age=rcrd['donor_age'],
                death_cause_group=rcrd[cn.DEATH_CAUSE_GROUP],
                malignancy=rcrd[cn.D_MALIGNANCY],
                tumor_history=rcrd[cn.D_TUMOR_HISTORY],
                donor_marginal_free_text=False,
                drug_abuse=rcrd[cn.D_DRUG_ABUSE],
                n_kidneys_available=2,
                hla=rcrd[cn.DONOR_HLA],
                hla_system=hla_system,
                diabetes=rcrd[cn.D_DIABETES],
                cardiac_arrest=(
                    1 if rcrd[cn.D_CARREST] == 'Yes' else 0
                ),
                last_creat=rcrd[cn.D_LAST_CREAT],
                smoker=False,
                hbsag=rcrd[cn.D_HBSAG],
                hcvab=rcrd[cn.D_HCVAB],
                hbcab=rcrd[cn.D_HBCAB],
                sepsis=rcrd[cn.D_SEPSIS],
                meningitis=rcrd[cn.D_MENINGITIS],
                hypertension=rcrd[cn.D_HYPERTENSION],
                euthanasia=False,
                rescue=rcrd[cn.RESCUE],
                urine_protein=False,
                cmv=False,
                sim_set=ss,
                donor_subregion_esp=rcrd[cn.DONOR_CENTER]
            )

            ml = MatchListCurrentETKAS(
                patients=[pat],
                donor=don,
                match_date=rcrd['date_transplanted'].to_pydatetime(),
                hla_system=hla_system,
                bal_system=bal_system,
                calc_points=ss.calc_etkas_score,
                sim_start_date=ss.SIM_START_DATE
            )

            # Print match info.
            if verbosity >= 1:
                print(f'***** Record {id} *******')
                print(f'Original scale param: {rcrd["b"]}')
                print(f'Original shape param: {rcrd["a"]}')
                print(f'Original surv prob: {rcrd["surv_prob"]}')

            if ml.return_match_list():
                offer = ml.return_match_list()[0]
                if rcrd[cn.RESCUE]:
                    offer.set_acceptance(reason=cn.T3)
                else:
                    offer.set_acceptance(reason=cn.T1)
                if hasattr(offer, '_initialize_acceptance_information'):
                    offer._initialize_acceptance_information()
                    offer._initialize_posttxp_information(ptp)

                surv_prob = ptp.calculate_survival(
                    offer=offer,
                    time=rcrd['pred_time'],
                    verbosity=verbosity
                )

                assert (surv_prob - rcrd['surv_prob']) <= slack, \
                    f'Survival probability incorrect for test case {id}'
            else:
                n_bg_incompatible += 1

        assert n_bg_incompatible < d_txps.shape[0] / 100 * 5, \
            f'{n_bg_incompatible}/{d_txps.shape[0]} empty ' \
            f'match lists. Due to BG incompatibility?'

    @unittest.skip('Skip')
    def test_posttxp(self) -> None:
        self.verify_posttxp_surv(
            verbosity=0,
            path_test_cases=es.POSTTXP_SURV_TESTPATHS[None],
            slack=0.01
        )

    def test_sim_posttxp_surv(
            self,
            verbosity=0
    ) -> None:
        """Assess whether code works for posttxp
        """

        # Read in simulation settings
        ss = read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )

        d_txps = pd.read_csv(es.POSTTXP_SURV_TESTPATHS[None])
        d_txps['date_transplanted'] = d_txps['date_txp'].apply(
            lambda x: datetime.strptime(
                x,
                ifs.DEFAULT_DATE_TIME_HMS
            )
        ).dt.to_pydatetime()
        d_txps[cn.DATE_FIRST_DIAL] = (
            d_txps['date_transplanted'] -
            d_txps[cn.YEARS_ON_DIAL].map(
                lambda x: timedelta(days=x) * es.DAYS_PER_YEAR_FLOAT
            )
        )
        d_txps.current_patient_hla = fix_hla_string(
            d_txps.current_patient_hla
        )
        hla_system = HLASystem(ss)
        mmp_system = MMPSystem(ss, hla_system)
        bal_system = load_balances(ss, update_balances=False)
        ptp = PostTransplantPredictor(
            offset_ids_transplants=99999,
            seed=14
        )

        n_bg_incompatible = 0
        for id, rcrd in d_txps.iterrows():

            pat = Patient(
                id_recipient=1,
                date_first_dial=rcrd[cn.DATE_FIRST_DIAL],
                r_dob=(
                    rcrd['date_transplanted'] -
                    timedelta(days=rcrd[cn.R_MATCH_AGE] * 365.24)
                ),
                recipient_country=rcrd[cn.PATIENT_COUNTRY],
                recipient_region='Netherlands',
                recipient_center=rcrd[cn.PATIENT_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                listing_date=rcrd['date_transplanted'] - timedelta(days=100),
                urgency_code='T',
                sex=rcrd[cn.R_SEX],
                sim_set=ss,
                time_since_prev_txp=rcrd[cn.TIME_SINCE_PREV_TXP],
                type_retx=(
                    cn.RETRANSPLANT_DURING_SIM if rcrd[cn.IS_RETRANSPLANT]
                    else cn.NO_RETRANSPLANT
                ),
                hla=rcrd[cn.CURRENT_PATIENT_HLA],
                hla_system=hla_system,
                mmp_system=mmp_system
            )

            don = Donor(
                id_donor=1255,
                donor_country=rcrd[cn.DONOR_COUNTRY],
                donor_region='Germany',
                donor_center=rcrd[cn.DONOR_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                reporting_date=rcrd['date_transplanted'],
                weight=80,
                height=180,
                donor_dcd=rcrd[cn.D_DCD],
                age=rcrd['donor_age'],
                death_cause_group=rcrd[cn.DEATH_CAUSE_GROUP],
                malignancy=rcrd[cn.D_MALIGNANCY],
                tumor_history=rcrd[cn.D_TUMOR_HISTORY],
                donor_marginal_free_text=False,
                drug_abuse=rcrd[cn.D_DRUG_ABUSE],
                n_kidneys_available=2,
                hla=rcrd[cn.DONOR_HLA],
                hla_system=hla_system,
                diabetes=rcrd[cn.D_DIABETES],
                cardiac_arrest=(
                    1 if rcrd[cn.D_CARREST] == 'Yes' else 0
                ),
                last_creat=rcrd[cn.D_LAST_CREAT],
                smoker=False,
                hbsag=rcrd[cn.D_HBSAG],
                hcvab=rcrd[cn.D_HCVAB],
                hbcab=rcrd[cn.D_HBCAB],
                sepsis=rcrd[cn.D_SEPSIS],
                meningitis=rcrd[cn.D_MENINGITIS],
                hypertension=rcrd[cn.D_HYPERTENSION],
                euthanasia=False,
                rescue=rcrd[cn.RESCUE],
                urine_protein=False,
                cmv=False,
                sim_set=ss,
                donor_subregion_esp=rcrd[cn.DONOR_CENTER]
            )

            ml = MatchListCurrentETKAS(
                patients=[pat],
                donor=don,
                match_date=rcrd['date_transplanted'].to_pydatetime(),
                hla_system=hla_system,
                bal_system=bal_system,
                calc_points=ss.calc_etkas_score,
                sim_start_date=ss.SIM_START_DATE
            )

            # Print match info.
            if verbosity >= 1:
                print(f'***** Record {id} *******')
                print(f'Original scale param: {rcrd["b"]}')
                print(f'Original shape param: {rcrd["a"]}')
                print(f'Original surv prob: {rcrd["surv_prob"]}')

            if ml.return_match_list():
                offer = ml.return_match_list()[0]
                if rcrd[cn.RESCUE]:
                    offer.set_acceptance(reason=cn.T3)
                else:
                    offer.set_acceptance(reason=cn.T1)
                if hasattr(offer, '_initialize_acceptance_information'):
                    offer._initialize_acceptance_information()
                    offer._initialize_posttxp_information(ptp)

                # Simulate post-transplant survival
                date_fail, date_relist, cause_fail = (
                    ptp.simulate_posttransplant(
                        offer=offer,
                        current_date=(
                            rcrd['date_transplanted'].to_pydatetime()
                        )
                    )
                )
            else:
                n_bg_incompatible += 1

        assert n_bg_incompatible < d_txps.shape[0] / 100 * 5, \
            f'{n_bg_incompatible}/{d_txps.shape[0]} empty ' \
            f'match lists. Due to BG incompatibility?'


if __name__ == '__main__':
    unittest.main()
