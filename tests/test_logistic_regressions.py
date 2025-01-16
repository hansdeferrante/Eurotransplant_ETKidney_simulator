import sys
import unittest
import os
from datetime import timedelta, datetime
import pandas as pd
from random import randint

sys.path.append('./')


from simulator.code.AcceptanceModule import AcceptanceModule
from simulator.code.utils.load_simulation_entities import load_balances
from simulator.code.entities import Patient, Donor
from simulator.code.HLA.HLASystem import HLASystem, Unacceptables
from simulator.code.HLA.MMPSystem import MMPSystem
import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es
from simulator.code.utils.read_input_files import \
    read_sim_settings, read_travel_times
from simulator.code.matchlist.MatchListETKASandESP import \
    MatchListCurrentETKAS
from simulator.code.matchlist.MatchListETKASandESP import \
    MatchListESP


class TestCurrentAllocation(unittest.TestCase):
    """Test whether current allocation is done correctly.
    """

    # @unittest.skip('Skip for now')
    def acceptance_centers_test(
        self, test_data_file: str,
        verbose: bool = False
    ):
        """Test center acceptances
        """

        if 'ESP' in test_data_file:
            ESP = True
        else:
            ESP = False

        # Read in simulation settings
        ss = read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )
        hla_system = HLASystem(sim_set=ss)
        mmp_system = MMPSystem(hla_system=hla_system, sim_set=ss)
        bal_system = load_balances(ss)
        travel_time_dict = read_travel_times()

        d_acc_records = pd.read_csv(test_data_file)
        acc_module = AcceptanceModule(
            seed=1,
            patient_acc_policy='LR',
            center_acc_policy='LR',
            verbose=0,
            simulate_random_effects=False
        )

        dummy_date = datetime(year=2000, month=1, day=1)

        dummy_hla = 'A19 A28 B12 B13 DR7 DR53'

        # d_acc_records = d_acc_records.loc[617:,]

        n_bg_incompatible = 0
        for id, rcrd in d_acc_records.iterrows():

            pat = Patient(
                id_recipient=randint(0, 10000),
                id_reg=randint(0, 10000),
                r_dob=(
                    dummy_date -
                    40 * es.DAYS_PER_YEAR
                ),
                recipient_country=rcrd[cn.PATIENT_COUNTRY],
                recipient_region=rcrd[cn.PATIENT_REGION],
                recipient_center=rcrd[cn.PATIENT_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                listing_date=dummy_date,
                urgency_code=cn.T,
                sex='Male',
                sim_set=ss,
                type_retx=(
                    cn.NO_RETRANSPLANT
                ),
                hla_system=hla_system,
                mmp_system=mmp_system,
                hla=dummy_hla,
                date_first_dial=dummy_date - es.DAYS_PER_YEAR,
                time_since_prev_txp=6 * 300
            )

            don = Donor(
                id_donor=1255,
                donor_country=rcrd[cn.DONOR_COUNTRY],
                donor_region=rcrd[cn.DONOR_REGION],
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
                smoker=rcrd.get(cn.D_SMOKING, False),
                hbsag=rcrd.get(cn.D_HBSAG, False),
                hcvab=rcrd.get(cn.D_HCVAB, False),
                hbcab=rcrd.get(cn.D_HBCAB, False),
                sepsis=rcrd.get(cn.D_SEPSIS, False),
                meningitis=rcrd.get(cn.D_MENINGITIS),
                hypertension=rcrd[cn.D_HYPERTENSION],
                euthanasia=rcrd.get(cn.D_EUTHANASIA, False),
                rescue=rcrd.get(cn.RESCUE),
                cmv=rcrd.get(cn.D_CMV, 0),
                urine_protein=rcrd.get(cn.D_URINE_PROTEIN, 0),
                sim_set=ss,
                donor_subregion_esp=[cn.DONOR_COUNTRY]
            )

            if ESP:
                ml = MatchListESP(
                    patients=[pat],
                    donor=don,
                    match_date=dummy_date,
                    hla_system=hla_system,
                    bal_system=bal_system,
                    calc_points=ss.calc_esp_score,
                    sim_start_date=ss.SIM_START_DATE,
                    travel_time_dict=travel_time_dict
                )
            else:
                ml = MatchListCurrentETKAS(
                    patients=[
                        pat
                    ],
                    donor=don,
                    match_date=dummy_date,
                    hla_system=hla_system,
                    bal_system=bal_system,
                    calc_points=ss.calc_etkas_score,
                    sim_start_date=ss.SIM_START_DATE,
                    travel_time_dict=travel_time_dict
                )

            pd.set_option('display.max_rows', None)

            # Check for patient offer whether prob is correct.
            if ml.return_match_list():
                offer = ml.return_match_list()[0]
                if hasattr(offer, '_initialize_acceptance_information'):
                    offer._initialize_acceptance_information()

                offer.__dict__[cn.PROFILE_COMPATIBLE] = (
                    rcrd.get(cn.PROFILE_COMPATIBLE, True)
                )
                offer.__dict__[cn.K_PREVIOUS_CENTER_REJECTIONS] = (
                    rcrd[cn.K_PREVIOUS_CENTER_REJECTIONS]
                )
                offer.__dict__[cn.ZERO_MISMATCH] = (
                    rcrd[cn.ZERO_MISMATCH]
                )
                offer.__dict__[cn.INTERNATIONAL_RESCUE] = (
                    rcrd[cn.INTERNATIONAL_RESCUE]
                )

                diff_prob = (
                    acc_module.calc_prob_center_willing_to_accept(
                        offer,
                        verbose=0,
                        selected_model='etkas_cd' if not ESP else 'esp_cd'
                    ) - rcrd.phat
                )
                assert abs(diff_prob) <= 1e-2, (
                    f'Turndown probability incorrect '
                    f'for test case {id} ({diff_prob})'
                )

            else:
                n_bg_incompatible += 1

        assert n_bg_incompatible < d_acc_records.shape[0] / 100 * 2, \
            f'{n_bg_incompatible}/{d_acc_records.shape[0]} empty match lists'

    def test_acceptance_centers(self):
        """Test whether acceptance in regular allocation is OK"""

        for test_settings in (
            'esp_cd', 'etkas_cd',
        ):
            print(f'Testing acceptance for: {test_settings}')
            self.acceptance_centers_test(
                es.LR_TEST_FILES.get(test_settings)
            )

    def acceptance_recipients_test(
            self,
            test_data_file: str,
            verbose: bool = 2
    ) -> None:
        """Test whether ML have acceptable donor/recipient combinations.
        """

        if 'ESP' in test_data_file:
            ESP = True
        else:
            ESP = False

        # Read in simulation settings
        ss = read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )

        d_acc_records = pd.read_csv(test_data_file)
        travel_time_dict = read_travel_times()
        acc_module = AcceptanceModule(
            seed=1,
            patient_acc_policy='LR',
            center_acc_policy='LR',
            verbose=0,
            simulate_random_effects=False
        )

        hla_system = HLASystem(sim_set=ss)
        mmp_system = MMPSystem(sim_set=ss, hla_system=hla_system)
        bal_system = load_balances(ss)

        dummy_date = datetime(year=2000, month=1, day=1)

        n_bg_incompatible = 0
        for id, rcrd in enumerate(d_acc_records.to_dict(orient='records')):

            pat = Patient(
                id_recipient=randint(0, 10000),
                id_reg=randint(0, 10000),
                r_dob=(
                    dummy_date -
                    timedelta(days=rcrd[cn.R_MATCH_AGE] *
                              es.DAYS_PER_YEAR_FLOAT + 15)
                ),
                recipient_country=rcrd[cn.PATIENT_COUNTRY],
                recipient_region=rcrd[cn.PATIENT_REGION],
                recipient_center=rcrd[cn.PATIENT_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                listing_date=dummy_date,
                urgency_code=(
                    cn.HU if rcrd.get(cn.PATIENT_IS_HU, False)
                    else cn.T
                ),
                sex='Male',
                sim_set=ss,
                type_retx=(
                    cn.NO_RETRANSPLANT if int(rcrd['retransplant']) == 0
                    else cn.RETRANSPLANT
                ),
                hla_system=hla_system,
                mmp_system=mmp_system,
                hla=rcrd[cn.PATIENT_HLA_LITERAL],
                date_first_dial=(
                    dummy_date -
                    timedelta(days=rcrd[cn.YEARS_ON_DIAL] *
                              es.DAYS_PER_YEAR_FLOAT)
                ),
                time_since_prev_txp=rcrd[cn.TIME_SINCE_PREV_TXP]
            )
            pat._vpra = rcrd[cn.VPRA]
            pat.unacceptable_antigens = Unacceptables(
                pat.hla_system,
                unacc_string=rcrd.get(cn.UNACC_ANT, '')
            )

            don = Donor(
                id_donor=1255,
                donor_country=rcrd[cn.DONOR_COUNTRY],
                donor_region=rcrd[cn.DONOR_REGION],
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
                hla=rcrd[cn.D_HLA_FULL],
                hla_system=hla_system,
                diabetes=rcrd[cn.D_DIABETES],
                cardiac_arrest=rcrd[cn.D_CARREST],
                last_creat=rcrd[cn.D_LAST_CREAT],
                smoker=rcrd[cn.D_SMOKING],
                hbsag=rcrd[cn.D_HBSAG],
                hcvab=rcrd[cn.D_HCVAB],
                hbcab=rcrd[cn.D_HBCAB],
                sepsis=rcrd[cn.D_SEPSIS],
                meningitis=rcrd[cn.D_MENINGITIS],
                hypertension=rcrd[cn.D_HYPERTENSION],
                euthanasia=rcrd[cn.D_EUTHANASIA],
                rescue=rcrd[cn.TXP_RESCUE],
                cmv=rcrd.get(cn.D_CMV, 0),
                urine_protein=rcrd.get(cn.D_URINE_PROTEIN, 0),
                sim_set=ss,
                donor_subregion_esp=[cn.DONOR_COUNTRY]
            )

            if ESP:
                ml = MatchListESP(
                    patients=[pat],
                    donor=don,
                    match_date=dummy_date,
                    hla_system=hla_system,
                    bal_system=bal_system,
                    calc_points=ss.calc_esp_score,
                    sim_start_date=ss.SIM_START_DATE,
                    travel_time_dict=travel_time_dict
                )
            else:
                ml = MatchListCurrentETKAS(
                    patients=[
                        pat
                    ],
                    donor=don,
                    match_date=dummy_date,
                    hla_system=hla_system,
                    bal_system=bal_system,
                    calc_points=ss.calc_etkas_score,
                    sim_start_date=ss.SIM_START_DATE,
                    travel_time_dict=travel_time_dict
                )

            # Check for patient offer whether prob is correct.
            if ml.return_match_list():
                offer = ml.return_match_list()[0]
                if cn.PROFILE_COMPATIBLE in rcrd:
                    offer.__dict__[cn.PROFILE_COMPATIBLE] = (
                        rcrd[cn.PROFILE_COMPATIBLE]
                    )
                else:
                    offer.__dict__[cn.PROFILE_COMPATIBLE] = True

                offer._initialize_acceptance_information()
                diff_prob = (
                    acc_module.calculate_prob_patient_accept(
                        offer, verbose=verbose
                    ) - rcrd['phat']
                )
                assert abs(diff_prob) <= 1e-2, (
                    f'Turndown probability incorrect for test '
                    f'case {id} phat: {round(rcrd["phat"], 3)}, '
                    f'delta: {round(diff_prob, 3)}).\n'
                )

            else:
                n_bg_incompatible += 1

        assert n_bg_incompatible < d_acc_records.shape[0] / 100 * 10 * 1.5, \
            f'{n_bg_incompatible}/{d_acc_records.shape[0]} empty ' \
            f'match lists.'

    def test_acceptance_patients(self):
        """Test whether acceptance in regular allocation is OK"""

        for test_settings in (
            'etkas_rd', 'esp_rd',
        ):
            print(f'Testing acceptance for: {test_settings}')
            self.acceptance_recipients_test(
                es.LR_TEST_FILES.get(test_settings),
                verbose=0
            )

    # @unittest.skip('Skipping dkt test')
    def test_dkt(self, verbose: bool = False):
        """Test whether dkt transplantation after acceptance
        """

        # Read in simulation settings
        ss = read_sim_settings(
            os.path.join(
                es.DIR_SIM_SETTINGS,
                'sim_settings_test.yaml'
            )
        )
        hla_system = HLASystem(sim_set=ss)
        mmp_system = MMPSystem(sim_set=ss, hla_system=hla_system)
        bal_system = load_balances(ss)

        d_acc_records = pd.read_csv('data/test/acceptance_dkt.csv')
        acc_module = AcceptanceModule(
            seed=1,
            patient_acc_policy='LR',
            center_acc_policy='LR',
            verbose=0,
            simulate_random_effects=False
        )

        dummy_date = datetime(year=2000, month=1, day=1)

        dummy_hla = 'A19 A28 B12 B13 DR7 DR53'

        n_bg_incompatible = 0
        travel_time_dict = read_travel_times()
        for id, rcrd in d_acc_records.iterrows():

            pat = Patient(
                id_recipient=randint(0, 10000),
                id_reg=randint(0, 10000),
                r_dob=(
                    dummy_date -
                    timedelta(days=rcrd[cn.R_MATCH_AGE] *
                              es.DAYS_PER_YEAR_FLOAT + 10)
                ),
                recipient_country=rcrd[cn.PATIENT_COUNTRY],
                recipient_region=rcrd[cn.PATIENT_REGION],
                recipient_center=rcrd[cn.PATIENT_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                listing_date=dummy_date,
                urgency_code=cn.T,
                sex='Male',
                sim_set=ss,
                type_retx=(
                    cn.NO_RETRANSPLANT
                ),
                hla_system=hla_system,
                mmp_system=mmp_system,
                hla=dummy_hla,
                date_first_dial=dummy_date -
                timedelta(days=1 * es.DAYS_PER_YEAR_FLOAT),
                time_since_prev_txp=6 * 300
            )

            don = Donor(
                id_donor=1255,
                donor_country=rcrd[cn.DONOR_COUNTRY],
                donor_region=rcrd[cn.DONOR_REGION],
                donor_center=rcrd[cn.DONOR_CENTER],
                bloodgroup=rcrd[cn.D_BLOODGROUP],
                reporting_date=dummy_date,
                weight=rcrd.get(cn.DONOR_WEIGHT, 70),
                height=rcrd.get(cn.DONOR_HEIGHT, 180),
                donor_dcd=rcrd.get(cn.D_DCD, False),
                age=rcrd['donor_age'],
                death_cause_group=rcrd.get(cn.DEATH_CAUSE_GROUP, 'Anoxia'),
                malignancy=rcrd.get(cn.D_MALIGNANCY, False),
                tumor_history=False,
                donor_marginal_free_text=False,
                drug_abuse=rcrd.get(cn.D_DRUG_ABUSE, False),
                n_kidneys_available=2,
                hla=dummy_hla,
                hla_system=hla_system,
                diabetes=rcrd.get(cn.D_DIABETES, False),
                cardiac_arrest=rcrd.get(cn.D_CARREST, False),
                last_creat=rcrd.get(cn.D_LAST_CREAT, 1),
                smoker=rcrd.get(cn.D_SMOKING, False),
                hbsag=rcrd.get(cn.D_HBSAG, False),
                hcvab=rcrd.get(cn.D_HCVAB, False),
                hbcab=rcrd.get(cn.D_HBCAB, False),
                sepsis=rcrd.get(cn.D_SEPSIS, False),
                meningitis=rcrd.get(cn.D_MENINGITIS),
                hypertension=rcrd.get(cn.D_HYPERTENSION, False),
                euthanasia=rcrd.get(cn.D_EUTHANASIA, False),
                rescue=False,
                cmv=rcrd.get(cn.D_CMV, False),
                urine_protein=rcrd.get(cn.D_URINE_PROTEIN, 0),
                sim_set=ss,
                donor_subregion_esp=[cn.DONOR_COUNTRY]
            )

            ml = MatchListCurrentETKAS(
                patients=[
                    pat
                ],
                donor=don,
                match_date=dummy_date,
                hla_system=hla_system,
                bal_system=bal_system,
                calc_points=ss.calc_etkas_score,
                sim_start_date=ss.SIM_START_DATE,
                travel_time_dict=travel_time_dict
            )

            # Check for patient offer whether prob is correct.
            if ml.return_match_list():
                offer = ml.return_match_list()[0]
                offer.donor.rescue = (
                    rcrd[cn.ACCEPTANCE_REASON] == 'T3'
                )
                offer.set_acceptance(rcrd[cn.ACCEPTANCE_REASON])

                if cn.PROFILE_COMPATIBLE in rcrd:
                    offer.__dict__[cn.PROFILE_COMPATIBLE] = (
                        rcrd[cn.PROFILE_COMPATIBLE]
                    )
                else:
                    offer.__dict__[cn.PROFILE_COMPATIBLE] = True

                offer._initialize_acceptance_information()

                diff_prob = (
                    acc_module.calc_prob_dkt(
                        offer, verbose=verbose
                    ) - rcrd['phat']
                )
                assert abs(diff_prob) <= 1e-2, (
                    f'Turndown probability incorrect for test '
                    f'case {id} ' f'(phat: {round(rcrd["phat"], 3)}, '
                    f'delta: {round(diff_prob, 3)})'
                )

            else:
                n_bg_incompatible += 1


if __name__ == '__main__':
    unittest.main()
