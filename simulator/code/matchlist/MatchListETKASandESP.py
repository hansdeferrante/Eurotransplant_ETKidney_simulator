from typing import List, Any, Tuple, Dict, Optional
import pandas as pd

import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
from simulator.code.PostTransplantPredictor import PostTransplantPredictor
from simulator.code.matchlist.MatchList import MatchList, MatchRecord
from simulator.code.HLA.HLASystem import HLASystem
from simulator.code.entities import (
    Patient, Profile
)
from simulator.code.BalanceSystem import BalanceSystem
from simulator.code.utils.utils import round_to_decimals
from simulator.code.matchlist.ScoringFunction import MatchPointFunction


class PatientMatchRecord:
    """ Class implementing a patient 'MatchRecord'. Points
        are calculated for patient-match records based on
        patient items only. Patient-match records are used
        to maintain a sorted list of active patients per
        blood group, which speeds up sorting generated match lists.
    """

    def __init__(
        self, patient: Patient,
        match_time: float,
        calc_points: MatchPointFunction,
        bal_system: BalanceSystem
    ) -> None:

        self.id_reg = patient.id_registration

        # National balance points
        self.__dict__[cn.BALANCE_NAT] = bal_system.return_national_balances(
            normalize=True, current_time=match_time, group_values=('median',)
        )[patient.__dict__[cn.RECIPIENT_COUNTRY]]

        # Determine years the patient has been on dialysis
        if (dial_time := patient.get_dial_time_sim_start()) is not None:
            if match_time > dial_time:
                self.__dict__[cn.YEARS_ON_DIAL] = (
                    match_time -
                    dial_time
                ) / es.DAYS_PER_YEAR_FLOAT
            else:
                self.__dict__[cn.YEARS_ON_DIAL] = 0
        else:
            self.__dict__[cn.YEARS_ON_DIAL] = 0

        # Mismatch frequency (counted)
        if patient.sim_set.determine_mismatchfreqs:
            self.__dict__.update(patient.get_hla_mismatchfreqs())

        # Mismatch probability
        if cn.ET_MMP in calc_points.variables:
            self.__dict__[cn.ET_MMP] = patient.get_et_mmp()

        # Add vPRA
        self.__dict__[cn.VPRA] = patient.vpra

        # Add previously accrued wait-time for re-transplant candidates
        self.__dict__[cn.PREVIOUS_WT] = patient.__dict__[cn.PREVIOUS_WT]
        self.__dict__[cn.R_PED] = patient.__dict__.get('_pediatric', False)

        self.total_match_points = calc_points.calc_patient_score(
            self.__dict__
        )

    def __str__(self):
        return f'Patient {self.id_reg} with {self.total_match_points} pts'


class MatchRecordETKAS(MatchRecord):
    """
    Class which implements a match record for ETKAS.

    Attributes
    ----------
    patient : Patient
        Patient information.
    donor : Donor
        Donor information.
    match_date : datetime
        Date that match list is generated.
    """

    def __init__(
        self,
        patient: Patient,
        *args,
        **kwargs
    ) -> None:

        # Construct match records for the patient.
        super(MatchRecordETKAS, self).__init__(
            patient=patient,
            *args, **kwargs
        )

        self.type_record = 'ETKAS'
        bal_system: BalanceSystem = kwargs['bal_system']
        hla_system: HLASystem = kwargs['hla_system']
        calc_points: MatchPointFunction = kwargs['calc_points']

        # Determine whether patient is pediatric.
        self._determine_pediatric()
        self._determine_match_distance()

        # Determine whether the candidate is on dialysis,
        # and how long the candidate has been on dialysis
        if (dial_time := self.patient.get_dial_time_sim_start()) is not None:
            if self.match_time > dial_time:
                self.__dict__[cn.YEARS_ON_DIAL] = (
                    self.match_time -
                    dial_time
                ) / es.DAYS_PER_YEAR_FLOAT
                self.__dict__[cn.ON_DIAL] = 1
            else:
                self.__dict__[cn.YEARS_ON_DIAL] = 0
                self.__dict__[cn.ON_DIAL] = 0
        else:
            self.__dict__[cn.YEARS_ON_DIAL] = 0
            self.__dict__[cn.ON_DIAL] = 0

        # If there are 0 mismatches between donor and candidate,
        # determine whether the donor is fully homozygous, and if so,
        # the homozygosity level of the recipient
        if self.__dict__[cn.ZERO_MISMATCH]:
            self.__dict__[cn.D_FULLY_HOMOZYGOUS] = (
                self.donor.hla.fully_homozygous
            )
            if self.__dict__[cn.D_FULLY_HOMOZYGOUS]:
                self.__dict__[cn.R_HOMOZYGOSITY_LEVEL] = (
                    self.patient.hla.homozygosity_level
                )
            else:
                self.__dict__[cn.R_HOMOZYGOSITY_LEVEL] = 0

        # Determine match tier
        self._check_zipped_rules(
            attr_name=cn.MTCH_TIER,
            rule_dict=self.sim_set['ETKAS_MATCH_TIERS'],
            value_required=True,
            default=es.DEFAULT_ETKAS_TIER
        )

        # Mismatch frequency (counted) or probability (analytical)
        # formula
        if self.sim_set.determine_mismatchfreqs:
            self.__dict__.update(patient.get_hla_mismatchfreqs())
        if cn.ET_MMP in calc_points.variables:
            self.__dict__[cn.ET_MMP] = patient.get_et_mmp()
        self.__dict__[cn.VPRA] = patient.vpra

        # Add balance points
        self.add_balance_points(bal_system=bal_system)

        # Determine additional needed features, if applicable.
        if not calc_points._vars_to_construct_initialized:
            calc_points.set_vars_to_construct(
                self.__dict__
            )
        self._mq = None
        if calc_points.vars_to_construct:
            for var in calc_points.vars_to_construct:
                try:
                    self.__dict__[var] = (
                        es.VARIABLE_CONSTRUCTION_FUNCTIONS[var](self)
                    )
                except Exception as e:
                    raise Exception(
                        f'Could not construct {var} for '
                        f'{self.__class__.__name__}'
                    )

        # Calculate match points. Also by component if asked for
        # in simulation settings.
        if self.store_score_components:
            self.sc = calc_points.calc_score_components(
                self.__dict__
            )
            if self.sc:
                self.total_match_points = sum(
                    self.sc.values()
                )
                self.__dict__.update(self.sc)
        else:
            self.total_match_points = calc_points.calc_score(
                self.__dict__
            )
            self.sc = None

        # Set match tuple
        self._match_tuple = None

    def add_balance_points(self, bal_system: BalanceSystem) -> None:
        """Add balance points to the match record."""
        sd = self.__dict__

        if bal_system.group_vars:
            bal_group = tuple(sd[k] for k in bal_system.group_vars)
        else:
            bal_group = ('1',)

        # National balance should always be inserted
        sd[cn.BALANCE_NAT] = bal_system.return_national_balances(
            normalize=True, current_time=self.match_time,
            group_values=bal_group
        )[sd[cn.RECIPIENT_COUNTRY]]

        # Regional balance may be inserted if recipient equals donor
        # country, and if they are in the regional balance countries
        if (
            sd[cn.D_ALLOC_COUNTRY] in es.COUNTRIES_REGIONAL_BALANCES and
            sd[cn.RECIPIENT_COUNTRY] in es.COUNTRIES_REGIONAL_BALANCES
        ):
            if (
                reg_bals := bal_system.return_regional_balances(
                    normalize=True,
                    current_time=self.match_time,
                    group_values=bal_group
                )
            ):
                sd[cn.BALANCE_REG] = (
                    reg_bals[sd[cn.RECIPIENT_COUNTRY]].get(
                        sd[cn.RECIPIENT_CENTER], 0
                    )
                )
            else:
                sd[cn.BALANCE_REG] = 0
        else:
            sd[cn.BALANCE_REG] = 0

    def _initialize_acceptance_information(self) -> None:
        """Initialize acceptance information for the match record."""
        # Initialize travel times
        if self.center_travel_times:
            self.__dict__.update(
                self.center_travel_times[
                    self.__dict__[cn.RECIPIENT_CENTER]
                ]
            )

        # Copy over patient and donor information, needed for acceptance.
        self.__dict__.update(self.donor.offer_inherit_cols)
        self.__dict__.update(self.patient.offer_inherit_cols)
        self.__dict__[cn.VPRA] = self.patient.vpra
        self.__dict__[cn.ANY_UNACC] = (
            1 if self.patient.unacceptable_antigens.unacceptables
            else 0
        )

        # Information copied over manually
        self.__dict__[cn.D_MALIGNANCY] = self.donor.malignancy
        self.__dict__[cn.D_DRUG_ABUSE] = self.donor.drug_abuse
        self.__dict__[cn.RESCUE] = self.donor.rescue
        self.__dict__[cn.INTERNATIONAL_RESCUE] = (
            (self.__dict__[cn.MATCH_DISTANCE] == cn.INT) and
            self.donor.rescue == 1
        )
        self.__dict__[cn.NONLOCAL_ESP_DONOR] = (
            self.__dict__[cn.MATCH_DISTANCE] != cn.L and
            self.__dict__[cn.D_AGE] >= self.donor.age_esp_eligible
        )
        self.__dict__[cn.VPRA_PERCENT] = self.__dict__[cn.VPRA] * 100

        # Information
        self.__dict__[cn.DONOR_AGE_75P] = int(
            self.__dict__[cn.D_AGE] >= 75
        )

        # Determine match abroad (part of acceptance)
        self._determine_match_abroad()

    def _initialize_posttxp_information(
        self, ptp: PostTransplantPredictor
    ) -> None:
        """Initialize post-transplant information for the match record."""
        # Date relative to 2014
        self.__dict__[cn.YEAR_TXP_RT2014] = (
            self.__dict__[cn.MATCH_DATE].year - 2014
        )

        self._calculate_posttxp_survival(ptp=ptp)

    def _calculate_posttxp_survival(
        self, ptp: PostTransplantPredictor
    ) -> None:
        """Calculate post-transplant survival probabilities."""
        for window in es.WINDOWS_TRANSPLANT_PROBS:
            self.__dict__[f'{es.PREFIX_SURVIVAL}_{window}'] = (
                round_to_decimals(
                    ptp.calculate_survival(
                        offer=self,
                        time=window
                    ),
                    3
                )
            )

    def _determine_pediatric(self) -> None:
        """Determine whether donor/patient are pediatric."""
        pat_dict = self.patient.__dict__

        # If patient is pediatric, or pediatric status is unknown,
        # update pediatric status
        if pat_dict['_pediatric'] is True or pat_dict['_pediatric'] is None:
            self.__dict__[cn.R_PED] = self.patient.get_pediatric(
                ped_fun=self.sim_set.check_r_pediatric,
                match_age=self.unrounded_match_age
            )
        else:
            self.__dict__[cn.R_PED] = pat_dict['_pediatric']

        if self.donor.__dict__[cn.D_PED] is None:
            self.__dict__[cn.D_PED] = self.donor.get_pediatric(
                ped_fun=self.sim_set.check_d_pediatric
            )
        else:
            self.__dict__[cn.D_PED] = self.donor.__dict__[cn.D_PED]

    def _determine_match_abroad(self) -> None:
        """Determine if the match is abroad."""
        self.__dict__[cn.MATCH_ABROAD] = (
            0 if (
                self.__dict__[cn.PATIENT_COUNTRY] ==
                self.__dict__[cn.DONOR_COUNTRY]
            )
            else 1
        )

    @property
    def match_tuple(self):
        if not self._match_tuple:
            self._match_tuple = self.return_match_tuple()
        return self._match_tuple

    @property
    def mq(self) -> Tuple[int, ...]:
        if self._mq is not None:
            return self._mq
        else:
            self._mq = tuple(
                self.__dict__[k] for k in es.MISMATCH_STR_DEFINITION
            )
        return self._mq

    @property
    def match_quality_compatible(self):
        if self._mq_compatible:
            return self._mq_compatible
        if type(self.patient.profile) is Profile:
            self._mq_compatible = self.patient.profile._check_hla_acceptable(
                mq=self.mq
            )
        else:
            self._mq_compatible = True
        return self._mq_compatible

    @property
    def offerable(self):
        if self.no_unacceptable_antigens is False:
            return False
        elif not self.donor.rescue and not self.other_profile_compatible:
            return False
        elif not self.donor.rescue and not self.match_quality_compatible:
            return False
        else:
            if (
                self.no_unacceptable_antigens and
                (
                    self.donor.rescue or
                    (
                        self.other_profile_compatible and
                        self.match_quality_compatible
                    )
                )
            ):
                return True

    def __lt__(self, other):
        """Compare match records for sorting."""
        return self.match_tuple > other.match_tuple


class MatchListCurrentETKAS(MatchList):
    """
    Class implementing a match list for the Current ETKAS system.

    Attributes
    ----------
    match_list : List[MatchRecordETKAS]
        List of match records.
    """

    def __init__(
            self,
            sort: bool = False,
            record_class=MatchRecordETKAS,
            store_score_components: bool = False,
            travel_time_dict: Optional[Dict[str, Any]] = None,
            *args,
            **kwargs
    ) -> None:
        super(
            MatchListCurrentETKAS, self
        ).__init__(
            sort=sort,
            record_class=record_class,
            attr_order_match=es.DEFAULT_ETKAS_ATTR_ORDER,
            store_score_components=store_score_components,
            travel_time_dict=travel_time_dict,
            *args,
            **kwargs
        )

        self.match_list.sort()
        self.sorted = True

        self.ext_alloc_priority = False

        for rank, match_record in enumerate(self.match_list):
            if type(match_record) is MatchRecordETKAS:
                match_record.add_patient_rank(rank)

    def return_match_list(
            self
    ) -> List[MatchRecord]:
        """Return the match list."""
        return self.match_list

    def return_match_info(
        self
    ) -> List[Dict[str, Any]]:
        """Return match information as a list of dictionaries."""
        return [
            matchr.return_match_info() for matchr in self.match_list
        ]


class MatchRecordESP(MatchRecord):
    """
    Class which implements a match record for current ESP.

    Attributes
    ----------
    patient : Patient
        Patient information.
    donor : Donor
        Donor information.
    match_date : datetime
        Date that match list is generated.
    """

    def __init__(
        self,
        patient: Patient,
        *args,
        **kwargs
    ) -> None:

        # Construct general match records for the patient.
        super(MatchRecordESP, self).__init__(
            patient=patient,
            *args, **kwargs
        )

        self.type_record = 'ESP'
        calc_points = kwargs['calc_points']

        # Determine general match distances
        self._determine_match_distance()

        # Determine specific match distance definitions for ESP
        self.determine_esp_patient_location()

        self.__dict__[cn.VPRA] = patient.vpra

        # Determine age
        self.__dict__[cn.PATIENT_ESP_AGED] = (
            self.__dict__[cn.R_MATCH_AGE] >=
            self.patient.age_esp_eligible
        )

        # Determine years the patient has been on dialysis
        if (dial_time := self.patient.get_dial_time_sim_start()) is not None:
            if self.match_time > dial_time:
                self.__dict__[cn.YEARS_ON_DIAL] = (
                    self.match_time -
                    dial_time
                ) / es.DAYS_PER_YEAR_FLOAT
                self.__dict__[cn.ON_DIAL] = 1
            else:
                self.__dict__[cn.YEARS_ON_DIAL] = 0
                self.__dict__[cn.ON_DIAL] = 0
        else:
            self.__dict__[cn.YEARS_ON_DIAL] = 0
            self.__dict__[cn.ON_DIAL] = 0

        # Add previously accrued wait-time for re-transplant candidates
        self.__dict__[cn.PREVIOUS_WT] = patient.__dict__[cn.PREVIOUS_WT]

        # Calculate match points (as well as components if necessary)
        if self.store_score_components:
            self.sc = calc_points.calc_score_components(
                self.__dict__
            )
            if self.sc:
                self.total_match_points = sum(
                    self.sc.values()
                )
                self.__dict__.update(self.sc)
        else:
            self.total_match_points = calc_points.calc_score(
                self.__dict__
            )
            self.sc = None

        self.determine_esp_tier()

        # Set match tuple
        self._match_tuple = None
        self._mq = None

    def _initialize_acceptance_information(self) -> None:
        """Initialize acceptance information for the match record."""
        # Initialize travel times to recipient centers
        if self.center_travel_times:
            self.__dict__.update(
                self.center_travel_times[
                    self.__dict__[cn.RECIPIENT_CENTER]
                ]
            )

        # Copy over patient and donor information, needed for acceptance.
        self.__dict__.update(self.donor.offer_inherit_cols)
        self.__dict__.update(self.patient.offer_inherit_cols)

        # Information copied over manually
        self.__dict__[cn.D_MALIGNANCY] = self.donor.malignancy
        self.__dict__[cn.D_DRUG_ABUSE] = self.donor.drug_abuse
        self.__dict__[cn.RESCUE] = self.donor.rescue
        self.__dict__[cn.INTERNATIONAL_RESCUE] = (
            (self.__dict__[cn.MATCH_DISTANCE] == cn.INT) and
            self.donor.rescue == 1
        )
        self.__dict__[cn.NONLOCAL_DCD_DONOR] = (
            (self.__dict__[cn.MATCH_DISTANCE] == cn.INT) &
            self.__dict__[cn.GRAFT_DCD]
        )
        self.__dict__[cn.VPRA_PERCENT] = self.patient.vpra * 100
        self.__dict__[cn.ANY_UNACC] = (
            1 if self.patient.unacceptable_antigens.unacceptables
            else 0
        )
        self.__dict__[cn.DONOR_AGE_75P] = int(
            self.__dict__[cn.D_AGE] >= 75
        )

        # Determine whether donor is from abroad
        self._determine_match_abroad()

    def _initialize_posttxp_information(
        self, ptp: PostTransplantPredictor
    ) -> None:
        """Initialize post-transplant information for the match record."""
        # Date relative to 2014
        self.__dict__[cn.YEAR_TXP_RT2014] = (
            self.__dict__[cn.MATCH_DATE].year - 2014
        )

        self._calculate_posttxp_survival(ptp=ptp)

    def _calculate_posttxp_survival(
        self, ptp: PostTransplantPredictor
    ) -> None:
        """Calculate post-transplant survival probabilities."""
        for window in es.WINDOWS_TRANSPLANT_PROBS:
            self.__dict__[f'{es.PREFIX_SURVIVAL}_{window}'] = (
                round_to_decimals(
                    ptp.calculate_survival(
                        offer=self,
                        time=window
                    ),
                    3
                )
            )

    def determine_esp_patient_location(self) -> None:
        """Determine the patient's location for ESP."""
        # Determine patient location for ESP

        self.__dict__[cn.MATCH_SUBREGIONAL] = int(
            self._same_esp_subregion(
                self.sim_set.ESP_SUBREGIONS
            )
        )
        match_geography = self.__dict__[cn.GEOGRAPHY_MATCH]
        rec_country = self.__dict__[cn.RECIPIENT_COUNTRY]

        if match_geography == cn.L and (
            rec_country in {mgr.AUSTRIA, mgr.BELGIUM}
        ):
            self.__dict__[cn.ESP_MATCH_LOCATION] = cn.L
        elif self.__dict__[cn.MATCH_SUBREGIONAL]:
            self.__dict__[cn.ESP_MATCH_LOCATION] = cn.S
        elif match_geography == mgr.R and rec_country == mgr.GERMANY:
            self.__dict__[cn.ESP_MATCH_LOCATION] = mgr.R
        elif self.__dict__[cn.MATCH_INTERNATIONAL] == 1:
            self.__dict__[cn.ESP_MATCH_LOCATION] = cn.A
        else:
            self.__dict__[cn.ESP_MATCH_LOCATION] = cn.H

        self.__dict__[cn.MATCH_LSR] = int(
            (self.__dict__[cn.ESP_MATCH_LOCATION] == cn.L) |
            (self.__dict__[cn.ESP_MATCH_LOCATION] == mgr.R) |
            (self.__dict__[cn.MATCH_SUBREGIONAL])
        )

    def determine_esp_tier(self) -> None:
        """Determine the ESP tier."""
        # Determine the ESP tier
        if self.__dict__[cn.DONOR_COUNTRY] == mgr.GERMANY:
            self.__dict__[cn.KIDNEY_PROGRAM] = (
                self.patient.__dict__[cn.KIDNEY_PROGRAM]
                if self.patient.__dict__[cn.KIDNEY_PROGRAM]
                else mgr.ETKAS
            )

        # Determine match rank tier (A-G)
        self._check_zipped_rules(
            attr_name=cn.MTCH_TIER,
            rule_dict=(
                self.sim_set['ESP_MATCH_TIERS'][
                    self.__dict__[cn.DONOR_COUNTRY]
                ]
            ),
            value_required=True
        )

        self.__dict__[cn.MTCH_TIER_REV] = mgr.ESP_TIER_REVERSAL_DICT[
            self.__dict__[cn.MTCH_TIER]
        ]

    def _same_esp_subregion(self, dict_esp_subregion) -> bool:
        """Check if donor and recipient are in the same ESP subregion."""
        if (
            self.__dict__[cn.D_ALLOC_COUNTRY] in mgr.COUNTRIES_SUBREGIONS and
            self.allocation_national
        ):
            d_ctr = self.donor.donor_subregion_esp
            r_ctr = self.__dict__[cn.RECIPIENT_CENTER]

            if (
                r_ctr in dict_esp_subregion and
                d_ctr in dict_esp_subregion
            ):

                if dict_esp_subregion[r_ctr] == dict_esp_subregion[d_ctr]:
                    return True
        return False

    def _determine_match_abroad(self) -> None:
        """Determine if the match is abroad."""
        if self.__dict__[cn.GEOGRAPHY_MATCH] == 'A':
            self.__dict__[cn.MATCH_ABROAD] = 1
        else:
            self.__dict__[cn.MATCH_ABROAD] = 0

    @property
    def match_tuple(self):
        if not self._match_tuple:
            self._match_tuple = self.return_match_tuple()
        return self._match_tuple

    @property
    def mq(self) -> Tuple[int, ...]:
        if self._mq is not None:
            return self._mq
        else:
            self._mq = tuple(
                self.__dict__[k] for k in es.MISMATCH_STR_DEFINITION
            )
        return self._mq

    @property
    def match_quality_compatible(self):
        if self._mq_compatible:
            return self._mq_compatible
        if type(self.patient.profile) is Profile:
            self._mq_compatible = self.patient.profile._check_hla_acceptable(
                mq=self.mq
            )
        else:
            self._mq_compatible = True
        return self._mq_compatible

    @property
    def offerable(self):
        if self.no_unacceptable_antigens is False:
            return False
        elif not self.donor.rescue and self.other_profile_compatible is False:
            return False
        elif (
            not self.donor.rescue and
            not self.__dict__[cn.MTCH_TIER] in mgr.FMR_TIERS_ESP
        ):
            return False
        elif (
            self.donor.rescue and
            self.patient.profile and
            (self.patient.profile.extended_esp == 0)
        ):
            return False
        else:
            return True

    def __lt__(self, other):
        """Compare match records for sorting."""
        return self.match_tuple > other.match_tuple

    def __str__(self):
        """Return a string representation of the match record."""
        points_str = (
            f'{str(self.total_match_points).rjust(4, " ")} '
            f'total match points'
        )
        try:
            fdial = self.patient.date_first_dial.strftime("%Y-%m-%d")
        except Exception as e:
            fdial = 'none'
        if cn.MTCH_TIER in self.__dict__:
            tier_str = f' in tier {self.mtch_tier_str.ljust(2)}'
        else:
            tier_str = ''
        accr = (
            self.__dict__.get(cn.ACCEPTANCE_REASON)
            if self.__dict__.get(cn.ACCEPTANCE_REASON)
            else ""
        )
        pat_age = self.__dict__[cn.R_MATCH_AGE]
        return (
            f'{self.determine_mismatch_string()} offer{tier_str} to '
            f'{str(self.__dict__[cn.ID_RECIPIENT]).rjust(6, "0")} '
            f'({self.__dict__[cn.RECIPIENT_CENTER]}) '
            f'on {self.date_match.strftime("%Y-%m-%d")} '
            f'from {self.__dict__[cn.D_ALLOC_CENTER]} '
            f'with date first dial: {fdial} '
            f'with {points_str} and age {pat_age} '
            f'({type(self).__name__}, {accr})'
        )


class MatchListESP(MatchList):
    """
    Class implementing a match list for the Current ESP system.

    Attributes
    ----------
    match_list : List[MatchRecordESP]
        List of match records.
    """

    def __init__(
        self,
        sort: bool = False,
        record_class=MatchRecordESP,
        store_score_components: bool = False,
        travel_time_dict: Optional[Dict[str, Any]] = None,
        *args,
        **kwargs
    ) -> None:
        super(
            MatchListESP, self
        ).__init__(
            sort=sort,
            record_class=record_class,
            attr_order_match=es.DEFAULT_ESP_ATTR_ORDER,
            store_score_components=store_score_components,
            travel_time_dict=travel_time_dict,
            *args,
            **kwargs
        )
        self.ext_alloc_priority = False

        self.match_list.sort()
        self.sorted = True
        if travel_time_dict:
            self.center_travel_times = travel_time_dict[
                self.donor.__dict__[cn.DONOR_CENTER]
            ]
        else:
            self.center_travel_times = None

        for rank, match_record in enumerate(self.match_list):
            if type(match_record) is MatchRecordESP:
                match_record.add_patient_rank(rank)

    def return_match_list(
            self
    ) -> List[MatchRecord]:
        """Return the match list."""
        return [m for m in self.match_list]

    def return_match_info(
        self
    ) -> List[Dict[str, Any]]:
        """Return match information as a list of dictionaries."""
        return [
            matchr.return_match_info() for matchr in self.match_list
        ]
