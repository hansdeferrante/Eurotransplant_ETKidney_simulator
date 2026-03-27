from __future__ import annotations

from collections import defaultdict
from math import nan
from typing import Any, Callable, Dict, FrozenSet, Iterable, List, Optional, Set, Tuple, Union

import numpy as np

import simulator.code.utils.read_input_files as rdr
import simulator.magic_values.column_names as cn
import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr

from .matching import HLAMatcher, MatchRules, PreparedTyping
from .mmp import MMPService
from .ontology import HLAOntology
from .typing import HLATyping
from .unacceptables import Unacceptables


class HLAStatsAPI:
    def __init__(
        self,
        sim_set,
        *,
        log_missing_codes: bool = False,
    ):
        self.sim_set = sim_set
        self.log_missing_codes = log_missing_codes

        self.ontology = HLAOntology(sim_set)
        self.rules = MatchRules.from_sim_set(sim_set)
        self.matcher = HLAMatcher(self.ontology, self.rules)
        self.mmp = MMPService(sim_set=sim_set)

        self._donor_pool_hlas: Optional[List[FrozenSet[str]]] = None
        self._structured_donor_pool: Optional[List[HLATyping]] = None
        self._match_potentials: Optional[Dict[str, Dict[str, float]]] = None

        self.nan_value = nan
        self.api = self

        self.needed_broad_mismatches = set(self.rules.needed_broad_mismatches)
        self.needed_split_mismatches = set(self.rules.needed_split_mismatches)
        self._required_loci = tuple(
            loc for loc in mgr.ALL_HLA_LOCI if loc in self.rules.required_loci
        )
        self.required_loci_indices = tuple(
            i for i, loc in enumerate(mgr.ALL_HLA_LOCI)
            if loc in self.rules.required_loci
        )
        self._all_loci_tuple = tuple(mgr.ALL_HLA_LOCI)
        self._all_loci_set = set(self._all_loci_tuple)
        self._split_indices = tuple(
            i for i, loc in enumerate(mgr.ALL_HLA_LOCI)
            if loc in self.needed_split_mismatches
        )

        self.hz_loci_broad = self.needed_broad_mismatches - self.needed_split_mismatches
        self.hz_loci_split = self.needed_split_mismatches
        self.loci_zero_mismatch = tuple(
            sim_set.get(
                "LOCI_ZERO_MISMATCH",
                (mgr.MMB_HLA_A, mgr.MMB_HLA_B, mgr.MMS_HLA_DR),
            )
        )

        self.missing_codes = defaultdict(int)
        self.mmp_loci_splits = set(self.mmp.mmp_loci_splits)
        self.mmp_loci_broads = set(self.mmp.mmp_loci_broads)
        self.mmp_loci = set(self.mmp.mmp_loci)

    def __getattr__(self, name: str):
        # Keep existing code simple by exposing ontology attributes.
        if hasattr(self.ontology, name):
            return getattr(self.ontology, name)
        raise AttributeError(name)

    @property
    def required_loci(self):
        return self._required_loci

    @staticmethod
    def list_match_quality_order() -> Tuple[str, ...]:
        return tuple(es.HLA_MQS_STR)

    def _load_bg_and_allele_frequencies(self) -> None:
        self.mmp._load_bg_and_allele_frequencies()

    @property
    def allele_frequencies_split(self) -> Dict[str, Dict[str, float]]:
        return self.mmp.allele_frequencies_split

    @property
    def allele_frequencies_broad(self) -> Dict[str, Dict[str, float]]:
        return self.mmp.allele_frequencies_broad

    @property
    def bg_frequencies(self) -> Dict[str, float]:
        return self.mmp.bg_frequencies

    def _normalize_hla_input_string(self, hla_string: str) -> str:
        return self.ontology._normalize_hla_input_string(hla_string)

    def lookup_alleles(self, input_generator, codes_to_broad, codes_to_split):
        for allele in input_generator:
            if (allele in codes_to_broad) and (allele in codes_to_split):
                yield (codes_to_broad[allele][1], codes_to_split[allele][1])
            elif allele in codes_to_broad:
                yield (codes_to_broad[allele][1], None)

    def return_all_antigens(self, input_string: str, expand: bool = False) -> List[str]:
        input_string_upp = input_string.upper()
        if self.log_missing_codes:
            for code in input_string_upp.split():
                if code not in self.ontology.all_match_codes:
                    self.missing_codes[code] += 1
        return self.ontology.return_all_antigens(input_string=input_string, expand=expand)

    def _count_broad_mismatches(
        self,
        d_hla: HLATyping,
        p_hla: HLATyping,
        loci: Optional[Set[str]] = None,
        locus_prefix: Optional[str] = "mmb_",
    ) -> Dict[str, Optional[int]]:
        loci_set = self._all_loci_set if loci is None else set(loci)
        prefix = "" if locus_prefix is None else locus_prefix
        return {
            f"{prefix}{sel_locus}": len(d_hla.broads[i] - p_hla.broads[i])
            for i, sel_locus in enumerate(self._all_loci_tuple)
            if sel_locus in loci_set
        }

    def _count_split_mismatches(
        self,
        d_hla: HLATyping,
        p_hla: HLATyping,
        loci: Optional[Set[str]] = None,
        locus_prefix: Optional[str] = "mms_",
    ) -> Dict[str, Optional[int]]:
        loci_set = self._all_loci_set if loci is None else set(loci)
        prefix = "" if locus_prefix is None else locus_prefix

        mm = {
            f"{prefix}{sel_locus}": len(d_hla.splits[i] - p_hla.splits[i])
            for i, sel_locus in enumerate(self._all_loci_tuple)
            if sel_locus in loci_set
        }

        for i_loc, locus in enumerate(self._all_loci_tuple):
            if locus in loci_set and (
                d_hla.has_unknown_splits[i_loc] or p_hla.has_unknown_splits[i_loc]
            ):
                common_broads = d_hla.broads[i_loc] & p_hla.broads[i_loc]
                p_match_hlas = set(p_hla.broads[i_loc])
                d_match_hlas = set(d_hla.broads[i_loc])

                for cb in common_broads:
                    if cb in self.ontology.splittable_broads[locus]:
                        p_splits = (
                            self.ontology.broads_to_splits[locus][cb] & p_hla.splits[i_loc]
                        )
                        if p_splits:
                            d_splits = (
                                self.ontology.broads_to_splits[locus][cb] & d_hla.splits[i_loc]
                            )
                            if d_splits:
                                p_match_hlas = (p_match_hlas | p_splits) - {cb}
                                d_match_hlas = (d_match_hlas | d_splits) - {cb}

                mm[f"{prefix}{locus}"] = len(d_match_hlas - p_match_hlas) if p_match_hlas else None

        return mm

    def count_mismatches_tuple(
        self,
        d_hla: HLATyping,
        p_hla: HLATyping,
        safely: bool = True,
    ) -> Optional[Tuple[Optional[int], ...]]:
        return self.matcher.count_mismatches(d_hla=d_hla, p_hla=p_hla, safely=safely)

    def count_mismatches_prepared(
        self,
        d_prepared: PreparedTyping,
        p_prepared: PreparedTyping,
        ignore_ambiguities: bool = True,
        safely: bool = True,
        nanv=nan,
    ) -> Optional[Dict[str, Union[int, float, None]]]:
        mm_tuple = self.matcher.count_mismatches_prepared(
            d=d_prepared,
            p=p_prepared,
            safely=safely,
        )
        if mm_tuple is None:
            return None

        mm = {
            key: mm_tuple[i]
            for i, key in enumerate(self.matcher.mismatch_output_keys)
        }

        if not ignore_ambiguities:
            d_hla = d_prepared.source_typing
            p_hla = p_prepared.source_typing
            for i in self._split_indices:
                locus = self._all_loci_tuple[i]
                key = f"mms_{locus}"
                if d_hla.has_unknown_splits[i] or p_hla.has_unknown_splits[i]:
                    mm[key] = nanv

        mm_total = 0
        for value in mm.values():
            if value is None:
                continue
            if isinstance(value, float) and np.isnan(value):
                continue
            mm_total += int(value)
        mm[cn.MM_TOTAL] = -mm_total
        return mm

    def count_mismatches(
        self,
        d_hla: HLATyping,
        p_hla: HLATyping,
        ignore_ambiguities: bool = True,
        safely: bool = True,
        nanv=nan,
    ) -> Optional[Dict[str, Union[int, float, None]]]:
        mm_tuple = self.count_mismatches_tuple(d_hla=d_hla, p_hla=p_hla, safely=safely)
        if mm_tuple is None:
            return None

        mm = {
            key: mm_tuple[i]
            for i, key in enumerate(self.matcher.mismatch_output_keys)
        }

        if not ignore_ambiguities:
            for i in self._split_indices:
                locus = self._all_loci_tuple[i]
                key = f"mms_{locus}"
                if d_hla.has_unknown_splits[i] or p_hla.has_unknown_splits[i]:
                    mm[key] = nanv

        mm_total = 0
        for value in mm.values():
            if value is None:
                continue
            if isinstance(value, float) and np.isnan(value):
                continue
            mm_total += int(value)
        mm[cn.MM_TOTAL] = -mm_total
        return mm

    def count_mismatches_dict(
        self,
        *,
        donor_hla: HLATyping,
        patient_hla: HLATyping,
        safely: bool = True,
    ) -> Optional[Dict[str, Optional[int]]]:
        return self.count_mismatches(
            d_hla=donor_hla,
            p_hla=patient_hla,
            safely=safely,
        )

    def count_broad_mismatches_dict(
        self,
        *,
        donor_hla: HLATyping,
        patient_hla: HLATyping,
        loci: Optional[Iterable[str]] = None,
        locus_prefix: str = "mmb_",
    ) -> Dict[str, int]:
        return self._count_broad_mismatches(
            d_hla=donor_hla,
            p_hla=patient_hla,
            loci=set(loci) if loci is not None else None,
            locus_prefix=locus_prefix,
        )

    def count_split_mismatches_dict(
        self,
        *,
        donor_hla: HLATyping,
        patient_hla: HLATyping,
        ignore_ambiguities: bool = True,
        loci: Optional[Iterable[str]] = None,
        locus_prefix: str = "mms_",
    ) -> Dict[str, Optional[float | int]]:
        mm = self._count_split_mismatches(
            d_hla=donor_hla,
            p_hla=patient_hla,
            loci=set(loci) if loci is not None else None,
            locus_prefix=locus_prefix,
        )
        if ignore_ambiguities:
            return mm

        mm_copy: Dict[str, Optional[float | int]] = dict(mm)
        loci_set = set(loci) if loci is not None else self._all_loci_set
        for i, loc in enumerate(self._all_loci_tuple):
            if loc in loci_set and (
                donor_hla.has_unknown_splits[i] or patient_hla.has_unknown_splits[i]
            ):
                mm_copy[f"{locus_prefix}{loc}"] = float("nan")
        return mm_copy

    @property
    def donor_pool_hlas(self) -> List[FrozenSet[str]]:
        if self._donor_pool_hlas is not None:
            return self._donor_pool_hlas
        if self.sim_set.PATH_DONOR_POOL is None:
            raise Exception(
                "Specify PATH_DONOR_POOL in simulation settings to calculate vPRAs."
            )

        df_don_pool = rdr.read_donor_pool(self.sim_set.PATH_DONOR_POOL)
        df_don_pool_hlas = df_don_pool.loc[:, [cn.ID_DONOR, cn.DONOR_HLA]].drop_duplicates()
        self._donor_pool_hlas = [
            frozenset(self.return_all_antigens(s, expand=True))
            for s in df_don_pool_hlas.donor_hla
        ]
        return self._donor_pool_hlas

    @property
    def structured_donor_pool_hlas(self) -> List[HLATyping]:
        if self._structured_donor_pool is not None:
            return self._structured_donor_pool
        if self.sim_set.PATH_DONOR_POOL is None:
            raise Exception(
                "Specify PATH_DONOR_POOL in simulation settings to calculate HLA probabilities."
            )

        df_don_pool = rdr.read_donor_pool(self.sim_set.PATH_DONOR_POOL)
        df_don_pool_hlas = df_don_pool.loc[:, [cn.ID_DONOR, cn.DONOR_HLA]].drop_duplicates()
        self._structured_donor_pool = [
            HLATyping(str(hla), ontology=self.ontology, verbose=False)
            for hla in df_don_pool_hlas.donor_hla
        ]
        return self._structured_donor_pool

    def calculate_vpra_from_string(self, unacc_str: Optional[str]) -> float:
        if not isinstance(unacc_str, str):
            return 0.0
        unacceptables = Unacceptables(ontology=self.ontology, unacc_string=unacc_str)
        return (
            sum(
                len(unacceptables.unacceptables & donor_antigens) > 0
                for donor_antigens in self.donor_pool_hlas
            ) / len(self.donor_pool_hlas)
        )

    def calculate_et1994_mmp_and_matchfrequency(
        self,
        *,
        hlas: Dict[str, Set[str]],
        bg: str,
        vpra: float,
        freqs: Dict[str, Dict[str, float]],
        n_donors: int = 1000,
    ) -> Tuple[float, float]:
        return self.mmp.calculate_et1994_mmp_and_matchfrequency(
            hlas=hlas,
            bg=bg,
            vpra=vpra,
            freqs=freqs,
            n_donors=n_donors,
        )

    def calculate_broad_split_mmps(
        self,
        *,
        candidate_profile: Optional[HLATyping] = None,
        candidate_bloodgroup: Optional[str] = None,
        vpra: Optional[float] = None,
        p: Optional[Any] = None,
    ) -> Tuple[Dict[str, float], Dict[str, float]]:
        if p is not None:
            candidate_profile = p.hla
            candidate_bloodgroup = getattr(p, "bloodgroup", None)
            if candidate_bloodgroup is None:
                candidate_bloodgroup = p.__dict__[cn.R_BLOODGROUP]
            vpra = p.vpra

        if candidate_profile is None or candidate_bloodgroup is None or vpra is None:
            raise TypeError(
                "Provide either (candidate_profile, candidate_bloodgroup, vpra) or p=Patient."
            )

        return self.mmp.calculate_broad_split_mmps(
            candidate_profile=candidate_profile,
            candidate_bloodgroup=candidate_bloodgroup,
            vpra=vpra,
        )

    def return_mm_dicts(
        self,
        patient_hla: HLATyping,
    ) -> List[Optional[Dict[str, Optional[int]]]]:
        return [
            self.count_mismatches(d_hla=d_hla, p_hla=patient_hla)
            for d_hla in self.structured_donor_pool_hlas
        ]

    def calculate_hla_match_probability(
        self,
        patient_hla: HLATyping,
        hla_match_pot_definitions: Dict[str, Callable],
    ) -> Dict[str, float]:
        mm_dicts = self.return_mm_dicts(patient_hla=patient_hla)
        hmps = {
            hmp_name: (
                sum(
                    hmp_fun(mm_dict)
                    for mm_dict in mm_dicts
                    if mm_dict is not None and hmp_fun(mm_dict) is not None
                )
                / sum(
                    1
                    for mm_dict in mm_dicts
                    if mm_dict is not None and hmp_fun(mm_dict) is not None
                )
                if sum(
                    1
                    for mm_dict in mm_dicts
                    if mm_dict is not None and hmp_fun(mm_dict) is not None
                ) > 0
                else np.nan
            )
            for hmp_name, hmp_fun in hla_match_pot_definitions.items()
        }
        return hmps

    @property
    def match_potentials(self):
        if self._match_potentials is not None:
            return self._match_potentials
        if self.sim_set.PATH_MATCH_POTENTIALS is None:
            raise Exception("Specify PATH_MATCH_POTENTIALS in simulation settings.")

        df_match_potentials = rdr.read_hla_match_potentials(self.sim_set.PATH_MATCH_POTENTIALS)
        self._match_potentials = df_match_potentials.set_index(cn.PATIENT_HLA).to_dict(orient="index")
        return self._match_potentials
