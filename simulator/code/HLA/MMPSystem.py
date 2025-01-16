#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

import yaml
import numpy as np
from typing import (
    List, Dict, Tuple, Optional, Callable, Set
)
from datetime import date
from math import isnan, prod
from collections import defaultdict

import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
import simulator.magic_values.column_names as cn
import simulator.code.utils.read_input_files as rdr

from simulator.code.HLA.HLASystem import HLASystem, HLAProfile
from simulator.code.entities import Patient, Donor
from simulator.code.utils.utils import (
    round_to_decimals, DotDict
)


class MMPSystem:
    """
        Class which implements several different types of
        definitions on how difficult a patient is to match.
        This includes ET's mismatch probability, a favorably
        matched haplotype frequency, and Gilk's matchability.
    ...

    Attributes   #noqa
    ----------


    Methods
    -------

    """

    def __init__(
        self,
        sim_set: DotDict,
        hla_system: HLASystem
    ):
        self.sim_set = sim_set
        self.hla_system = hla_system

        # Set loci for calculation of the mismatch probability
        self.mmp_loci_splits = set(
            sim_set.get('LOCI_MMP_SPLIT', (mgr.HLA_DR,))
        )
        self.mmp_loci_broads = set(
            sim_set.get('LOCI_MMP_BROAD', (mgr.HLA_A, mgr.HLA_B,))
        )
        self.mmp_loci = self.mmp_loci_splits | self.mmp_loci_broads

        # Load necessary allele frequencies for calculation
        # of mismatch probability
        self.load_bg_and_allele_frequencies(sim_set=sim_set)

        self._donor_pool_hlas = None
        self._structured_donor_pool = None
        self._match_potentials = None

    def load_bg_and_allele_frequencies(self, sim_set: DotDict):
        # Load allele frequencies and blood type frequencies
        # from external files.
        with open(
            sim_set.PATH_ALLELE_FREQUENCIES_SPLIT, "r", encoding='utf-8'
        ) as file:
            self.allele_frequencies_split: Dict[
                str, Dict[str, float]
            ] = yaml.load(
                file, Loader=yaml.FullLoader
            )
        with open(
            sim_set.PATH_ALLELE_FREQUENCIES_BROAD, "r", encoding='utf-8'
        ) as file:
            self.allele_frequencies_broad: Dict[
                str, Dict[str, float]
            ] = yaml.load(
                file, Loader=yaml.FullLoader
            )
        with open(
            sim_set.PATH_BG_FREQUENCIES, "r", encoding='utf-8'
        ) as file:
            self.bg_frequencies: Dict[str, float] = yaml.load(
                file, Loader=yaml.FullLoader
            )

    def _calculate_et1994_zero_or_one_mismatch(
            self,
            hlas: Dict[str, Set[str]],
            freqs: Dict[str, Dict[str, float]]
    ) -> Tuple[float, float]:
        """ This function calculates based on HLA frequency tables
            the probability that a randomly arriving donor has
            0 or 1 mismatch. Note that this function assumes that
            HLA antigens are independently distributed, i.e., this function
            ignores genetic linkage.
        """

        # Calculate per locus the probability that a random HLA matches to it.
        prob_match_per_locus = {
            locus: sum(
                freqs[locus].get(antigen, 1e-4)
                for antigen in list(antigens)
            ) for locus, antigens in hlas.items()
            if locus in self.mmp_loci
        }

        # Probability of 0 mismatches in total
        prob_both_match = {
            locus: prob_match**2
            for locus, prob_match in prob_match_per_locus.items()
        }
        MMP0 = 1
        for prob in prob_both_match.values():
            MMP0 = MMP0 * prob

        # Probability of exactly 1 mismatch per locus
        prob_1heterozygous_mismatch = {
            locus: x * (1 - x) + (1 - x) * x
            for locus, x in prob_match_per_locus.items()
        }

        # Probability of a homozygous mismatch
        all_hlas = set().union(*hlas)
        prob_1homozygous_mismatch = {
            locus: sum(
                freq ** 2 * (antigen not in all_hlas)
                for antigen, freq in freqs_per_locus.items()
            ) for locus, freqs_per_locus in freqs.items()
        }

        # Probability that there is no mismatch on other loci
        try:
            prob_no_mismatch_on_other_loci = {
                locus: MMP0 / prob_match**2 for
                locus, prob_match in prob_match_per_locus.items()
            }
        except Exception as e:
            prob_no_mismatch_on_other_loci = {
                locus: prod(
                    v for k, v in prob_match_per_locus.items()
                    if k != locus
                )**2
                for locus in prob_match_per_locus.keys()
            }

        prob_exactly_one_mismatch = {}
        for locus, prob_nomm_other_loci in (
            prob_no_mismatch_on_other_loci.items()
        ):
            prob_exactly_one_mismatch[locus] = (
                prob_nomm_other_loci *
                (
                    prob_1homozygous_mismatch[locus] +
                    prob_1heterozygous_mismatch[locus]
                )
            )

        MMP1 = sum(prob_exactly_one_mismatch.values())

        return MMP0, MMP1

    def _calculate_donor_mismatch_probability(
        self,
        prob_hla_match: float,
        bg: str,
        vpra: float,
        n_donors: int = 1000
    ) -> float:
        """Function which calculates the donor mismatch probability, based on
           a HLA match probability, candidate bloodgroup, and the candidate
           vPRA. This function assumes that unacceptable antigens, a HLA match,
           and BG are independently distributed. It does not ignore
           genetic linkage in calculating mismatch probabilities
        """
        prob_bg_match = self.bg_frequencies[bg]
        return round_to_decimals(
            100 * (
                (1 - (prob_bg_match * (1 - vpra) * (prob_hla_match)))**n_donors
            ),
            4
        )

    def _calculate_hla_mismatch_frequency(
        self,
        prob_hla_match: float,
        n_donors: int = 1000
    ) -> float:
        """Function which calculates the HLA mismatch frequency. This only
        takes into account probability of finding a good HLA match among
        the next 1000 donors"""
        return round_to_decimals(
            100 * ((1 - prob_hla_match)**n_donors),
            4
        )

    def calculate_et1994_mmp_and_matchfrequency(
        self,
        bg: str,
        vpra: float,
        n_donors: int = 1000,
        **kwargs
    ) -> Tuple[float, float]:
        """This function calculates the ET donor mismatch probability.
        It returns the ET donor MMP, and ET HLA-mismatch frequency.
        """

        MMP0, MMP1 = self._calculate_et1994_zero_or_one_mismatch(
            **kwargs
        )
        return (
            self._calculate_donor_mismatch_probability(
                prob_hla_match=MMP0 + MMP1,
                bg=bg,
                vpra=vpra,
                n_donors=n_donors
            ),
            self._calculate_hla_mismatch_frequency(
                prob_hla_match=MMP0 + MMP1,
                n_donors=n_donors
            )
        )

    def calculate_broad_split_mmps(
        self, p: Patient
    ) -> Tuple[Dict[str, float], Dict[str, float]]:
        """Calculate ET donor mismatch probabilities
            and ET donor mismatch frequencies.

        Note that both these definitions assume:
            (i) independence between HLA on the loci,
            (ii) independence between vPRA, BG, and HLA
        """
        mmps = dict()

        match_hlas_broad = {
            locus: p.hla.broads[i_loc]
            for i_loc, locus in enumerate(mgr.ALL_HLA_LOCI)
            if locus in self.mmp_loci
        }
        et_mmp_broad, et_matchfreq_broad = (
            self.calculate_et1994_mmp_and_matchfrequency(
                hlas=match_hlas_broad,
                bg=p.__dict__[cn.R_BLOODGROUP],
                vpra=p.vpra,
                freqs=self.allele_frequencies_broad
            )
        )

        match_hlas_split = {
            **match_hlas_broad,
            **{
                locus: p.hla.splits[i_loc]
                for i_loc, locus in enumerate(mgr.ALL_HLA_LOCI)
                if locus in self.mmp_loci_splits
            }
        }
        et_mmp_split, et_matchfreq_split = (
            self.calculate_et1994_mmp_and_matchfrequency(
                hlas=match_hlas_split,
                bg=p.__dict__[cn.R_BLOODGROUP],
                vpra=p.vpra,
                freqs=self.allele_frequencies_split
            )
        )

        mmps = {
            cn.ET_MMP_BROAD: et_mmp_broad,
            cn.ET_MMP_SPLIT: et_mmp_split
        }

        match_freqs = {
            cn.ET_MMP_BROAD: et_matchfreq_broad,
            cn.ET_MMP_SPLIT: et_matchfreq_split
        }

        return mmps, match_freqs

    def return_mm_dicts(
        self, p_hla: HLAProfile
    ) -> List[Optional[Dict[str, int]]]:

        return list(
            self.hla_system.count_mismatches(d_hla=d.hla, p_hla=p_hla)
            for d in self.structured_donor_pool_hlas
        )

    def calculate_hla_match_probability(
        self, p_hla: HLAProfile,
        hla_match_pot_definitions: Dict[str, Callable]
    ) -> Dict[str, float]:
        """
            Calculate HLA-match probability, i.e. the probability
            that a patient has a favorable match with the donor.

            The dictionary hla_match_pot_definitions contains
            as values Callables, which should indicate whether
            a particular match is favorably matched.

            For instance, the 1ABDR definition returns 1
            if there is at most 1 mismatch on the HLA-A, HLA-B,
            and HLA-DR loci.
        """

        mm_dicts = self.return_mm_dicts(
            p_hla=p_hla
        )

        hmps = {
            hmp_name: sum(
                hmp_fun(mm_dict) for mm_dict in mm_dicts
                if hmp_fun(mm_dict) is not None
            ) / sum(1 for mm_dict in mm_dicts if hmp_fun(mm_dict) is not None)
            if sum(
                1 for mm_dict in mm_dicts if hmp_fun(mm_dict) is not None
            ) > 0
            else np.nan
            for hmp_name, hmp_fun in hla_match_pot_definitions.items()
        }

        return hmps

    @property
    def structured_donor_pool_hlas(self):
        if self._structured_donor_pool is not None:
            return self._structured_donor_pool
        elif self.sim_set.PATH_DONOR_POOL is not None:
            df_don_pool = rdr.read_donor_pool(self.sim_set.PATH_DONOR_POOL)
            df_don_pool_hlas = df_don_pool.loc[
                :,
                [cn.ID_DONOR, cn.D_BLOODGROUP, cn.DONOR_HLA]
            ].drop_duplicates()
            self._structured_donor_pool = list(
                Donor.from_dummy_donor(
                    bloodgroup=bg,
                    hla=s,
                    hla_system=self.hla_system,
                    sim_set=self.sim_set
                )
                for s, bg in zip(
                    df_don_pool_hlas.donor_hla,
                    df_don_pool_hlas.d_bloodgroup
                )
            )
            return self._structured_donor_pool
        else:
            raise Exception(
                "Specify a path to a donor pool to "
                "calculate vPRAs in the yml file."
            )

    @property
    def match_potentials(self):
        if self._match_potentials is not None:
            return self._match_potentials
        elif self.sim_set.PATH_MATCH_POTENTIALS is not None:
            df_match_potentials = rdr.read_hla_match_potentials(
                self.sim_set.PATH_MATCH_POTENTIALS
            )
            self._match_potentials = (
                df_match_potentials.set_index(
                    cn.PATIENT_HLA
                ).to_dict(orient='index')
            )
            return self._match_potentials
        else:
            raise Exception(
                "Specify a path to match potentials to "
                "calculate vPRAs in the yml file."
            )
