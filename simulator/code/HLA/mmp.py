from __future__ import annotations

from math import prod
from typing import Dict, Optional, Set, Tuple

import yaml

import simulator.magic_values.column_names as cn
import simulator.magic_values.magic_values_rules as mgr
from simulator.code.utils.utils import round_to_decimals

from .typing import HLATyping


class MMPService:
    """ Calculator for MMP-like quantities. These are always calculated
        from the antigen frequencies
    """
    def __init__(self, sim_set):
        self.sim_set = sim_set
        self.mmp_loci_splits = set(sim_set.get("LOCI_MMP_SPLIT", (mgr.HLA_DR,)))
        self.mmp_loci_broads = set(
            sim_set.get("LOCI_MMP_BROAD", (mgr.HLA_A, mgr.HLA_B))
        )
        self.mmp_loci = self.mmp_loci_splits | self.mmp_loci_broads

        self._allele_frequencies_split: Optional[Dict[str, Dict[str, float]]] = None
        self._allele_frequencies_broad: Optional[Dict[str, Dict[str, float]]] = None
        self._bg_frequencies: Optional[Dict[str, float]] = None

    def _load_bg_and_allele_frequencies(self) -> None:
        if (
            self._allele_frequencies_split is not None
            and self._allele_frequencies_broad is not None
            and self._bg_frequencies is not None
        ):
            return

        with open(self.sim_set.PATH_ALLELE_FREQUENCIES_SPLIT, "r", encoding="utf-8") as file:
            self._allele_frequencies_split = yaml.load(file, Loader=yaml.FullLoader)
        with open(self.sim_set.PATH_ALLELE_FREQUENCIES_BROAD, "r", encoding="utf-8") as file:
            self._allele_frequencies_broad = yaml.load(file, Loader=yaml.FullLoader)
        with open(self.sim_set.PATH_BG_FREQUENCIES, "r", encoding="utf-8") as file:
            self._bg_frequencies = yaml.load(file, Loader=yaml.FullLoader)

    @property
    def allele_frequencies_split(self) -> Dict[str, Dict[str, float]]:
        self._load_bg_and_allele_frequencies()
        return self._allele_frequencies_split  # type: ignore[return-value]

    @property
    def allele_frequencies_broad(self) -> Dict[str, Dict[str, float]]:
        self._load_bg_and_allele_frequencies()
        return self._allele_frequencies_broad  # type: ignore[return-value]

    @property
    def bg_frequencies(self) -> Dict[str, float]:
        self._load_bg_and_allele_frequencies()
        return self._bg_frequencies  # type: ignore[return-value]

    def _calculate_et1994_zero_or_one_mismatch(
        self,
        hlas: Dict[str, Set[str]],
        freqs: Dict[str, Dict[str, float]],
    ) -> Tuple[float, float]:
        prob_match_per_locus = {
            locus: sum(freqs[locus].get(antigen, 1e-4) for antigen in list(antigens))
            for locus, antigens in hlas.items()
            if locus in self.mmp_loci
        }

        prob_both_match = {
            locus: prob_match ** 2
            for locus, prob_match in prob_match_per_locus.items()
        }
        mmp0 = 1.0
        for prob in prob_both_match.values():
            mmp0 *= prob

        prob_1heterozygous_mismatch = {
            locus: x * (1 - x) + (1 - x) * x
            for locus, x in prob_match_per_locus.items()
        }

        all_hlas = set().union(*hlas.values()) if hlas else set()
        prob_1homozygous_mismatch = {
            locus: sum(
                freq ** 2 * (antigen not in all_hlas)
                for antigen, freq in freqs_per_locus.items()
            )
            for locus, freqs_per_locus in freqs.items()
        }

        try:
            prob_no_mismatch_on_other_loci = {
                locus: mmp0 / (prob_match ** 2)
                for locus, prob_match in prob_match_per_locus.items()
            }
        except Exception:
            prob_no_mismatch_on_other_loci = {
                locus: prod(
                    value
                    for key, value in prob_match_per_locus.items()
                    if key != locus
                ) ** 2
                for locus in prob_match_per_locus.keys()
            }

        prob_exactly_one_mismatch = {}
        for locus, prob_nomm_other_loci in prob_no_mismatch_on_other_loci.items():
            prob_exactly_one_mismatch[locus] = (
                prob_nomm_other_loci
                * (
                    prob_1homozygous_mismatch[locus]
                    + prob_1heterozygous_mismatch[locus]
                )
            )

        mmp1 = sum(prob_exactly_one_mismatch.values())
        return mmp0, mmp1

    def _calculate_donor_mismatch_probability(
        self,
        prob_hla_match: float,
        bg: str,
        vpra: float,
        n_donors: int = 1000,
    ) -> float:
        prob_bg_match = self.bg_frequencies[bg]
        return round_to_decimals(
            100 * ((1 - (prob_bg_match * (1 - vpra) * prob_hla_match)) ** n_donors),
            4,
        )

    @staticmethod
    def _calculate_hla_mismatch_frequency(
        prob_hla_match: float,
        n_donors: int = 1000,
    ) -> float:
        return round_to_decimals(100 * ((1 - prob_hla_match) ** n_donors), 4)

    def calculate_et1994_mmp_and_matchfrequency(
        self,
        *,
        hlas: Dict[str, Set[str]],
        bg: str,
        vpra: float,
        freqs: Dict[str, Dict[str, float]],
        n_donors: int = 1000,
    ) -> Tuple[float, float]:
        mmp0, mmp1 = self._calculate_et1994_zero_or_one_mismatch(hlas=hlas, freqs=freqs)
        prob_hla_match = mmp0 + mmp1
        return (
            self._calculate_donor_mismatch_probability(
                prob_hla_match=prob_hla_match,
                bg=bg,
                vpra=vpra,
                n_donors=n_donors,
            ),
            self._calculate_hla_mismatch_frequency(
                prob_hla_match=prob_hla_match,
                n_donors=n_donors,
            ),
        )

    @staticmethod
    def _candidate_match_hlas(
        *,
        candidate_profile: HLATyping,
        mmp_loci: Set[str],
        mmp_loci_splits: Set[str],
    ) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
        match_hlas_broad = {
            locus: candidate_profile.broads[i_loc]
            for i_loc, locus in enumerate(mgr.ALL_HLA_LOCI)
            if locus in mmp_loci
        }
        match_hlas_split = {
            **match_hlas_broad,
            **{
                locus: candidate_profile.splits[i_loc]
                for i_loc, locus in enumerate(mgr.ALL_HLA_LOCI)
                if locus in mmp_loci_splits
            },
        }
        return match_hlas_broad, match_hlas_split

    def calculate_et_mmp(
        self,
        *,
        candidate_profile: HLATyping,
        candidate_bloodgroup: str,
        vpra: float,
    ) -> Optional[float]:
        try:
            _, match_hlas_split = self._candidate_match_hlas(
                candidate_profile=candidate_profile,
                mmp_loci=self.mmp_loci,
                mmp_loci_splits=self.mmp_loci_splits,
            )
            et_mmp, _ = self.calculate_et1994_mmp_and_matchfrequency(
                hlas=match_hlas_split,
                bg=candidate_bloodgroup,
                vpra=vpra,
                freqs=self.allele_frequencies_split,
            )
            return et_mmp
        except Exception:
            return None

    def calculate_broad_split_mmps(
        self,
        *,
        candidate_profile: HLATyping,
        candidate_bloodgroup: str,
        vpra: float,
    ) -> Tuple[Dict[str, float], Dict[str, float]]:
        match_hlas_broad, match_hlas_split = self._candidate_match_hlas(
            candidate_profile=candidate_profile,
            mmp_loci=self.mmp_loci,
            mmp_loci_splits=self.mmp_loci_splits,
        )
        et_mmp_broad, et_matchfreq_broad = self.calculate_et1994_mmp_and_matchfrequency(
            hlas=match_hlas_broad,
            bg=candidate_bloodgroup,
            vpra=vpra,
            freqs=self.allele_frequencies_broad,
        )

        et_mmp_split, et_matchfreq_split = self.calculate_et1994_mmp_and_matchfrequency(
            hlas=match_hlas_split,
            bg=candidate_bloodgroup,
            vpra=vpra,
            freqs=self.allele_frequencies_split,
        )

        mmps = {
            cn.ET_MMP_BROAD: et_mmp_broad,
            cn.ET_MMP_SPLIT: et_mmp_split,
        }
        match_freqs = {
            cn.ET_MMP_BROAD: et_matchfreq_broad,
            cn.ET_MMP_SPLIT: et_matchfreq_split,
        }
        return mmps, match_freqs
