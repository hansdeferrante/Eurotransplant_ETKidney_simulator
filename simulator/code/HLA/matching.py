from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Tuple

import simulator.magic_values.magic_values_rules as mgr

from .ontology import HLAOntology
from .typing import HLATyping


def _diff_len(dset: FrozenSet[str], pset: FrozenSet[str]) -> int:
    return len(dset - pset)


@dataclass(frozen=True, slots=True)
class MatchRules:
    """
    The specified rules for HLA matching in ET.
    """
    needed_broad_mismatches: FrozenSet[str]
    needed_split_mismatches: FrozenSet[str]

    @classmethod
    def from_sim_set(cls, sim_set) -> MatchRules:
        needed_broad = (
            frozenset(sim_set.NEEDED_BROAD_MISMATCHES)
            if getattr(sim_set, "NEEDED_BROAD_MISMATCHES", None) is not None
            else frozenset({mgr.HLA_A, mgr.HLA_B})
        )
        needed_split = (
            frozenset(sim_set.NEEDED_SPLIT_MISMATCHES)
            if getattr(sim_set, "NEEDED_SPLIT_MISMATCHES", None) is not None
            else frozenset({mgr.HLA_DR})
        )
        return cls(
            needed_broad_mismatches=needed_broad,
            needed_split_mismatches=needed_split,
        )

    @property
    def required_loci(self) -> FrozenSet[str]:
        return self.needed_broad_mismatches | self.needed_split_mismatches


@dataclass(slots=True)
class PreparedEpletView:
    """Lazy eplet payload attached to PreparedTyping per eplet signature."""

    eplets_by_definition_by_locus: Dict[str, Dict[str, FrozenSet[str]]]
    universe_by_locus: Dict[str, FrozenSet[str]]
    universe_by_group: Dict[str, Dict[str, FrozenSet[str]]]
    hlamm_comparison_universe_by_locus: Dict[str, FrozenSet[str]]
    eplets_per_allele_by_mode: Dict[str, Dict[str, Dict[str, FrozenSet[str]]]]


@dataclass(slots=True)
class PreparedTyping:
    """
    Cached view of HLATyping for fast matching. The specific
    view depends on the HLA matching rules. 
    All fields are derived once per (typing, ruleset).
    """

    match_broads: Tuple[FrozenSet[str], ...]
    match_splits: Tuple[FrozenSet[str], ...]

    required_loci_known: bool
    required_loci: Tuple[str, ...]
    needed_split_loci: Tuple[str, ...]

    broads_by_locus: Tuple[FrozenSet[str], ...]
    splits_by_locus: Tuple[FrozenSet[str], ...]
    alleles_by_locus: Tuple[FrozenSet[str], ...]
    has_unknown_splits: Tuple[bool, ...]
    has_unknown_match_splits: bool
    has_any_unknown_splits: bool
    homozygosity_per_locus: Dict[str, Optional[int]]
    homozygosity_level: int
    fully_homozygous: bool
    reduced_hla_match: str
    source_typing: HLATyping
    eplet_views: Dict[str, PreparedEpletView] = field(default_factory=dict)


class HLAMatcher:
    """
    Owns matching policy and algorithms.
    Depends on ontology + MatchRules. Uses PreparedTyping cached on HLATyping.
    """

    def __init__(self, ontology: HLAOntology, rules: MatchRules):
        self.ontology = ontology
        self.rules = rules

        # Precompute locus indices once
        self._broad_ix = [
            i
            for i, loc in enumerate(mgr.ALL_HLA_LOCI)
            if loc in rules.needed_broad_mismatches
        ]
        self._split_ix = [
            i
            for i, loc in enumerate(mgr.ALL_HLA_LOCI)
            if loc in rules.needed_split_mismatches
        ]

        # Precompute output keys once (avoid f-string in hot loop)
        self._broad_keys = [f"mmb_{mgr.ALL_HLA_LOCI[i]}" for i in self._broad_ix]
        self._split_keys = [f"mms_{mgr.ALL_HLA_LOCI[i]}" for i in self._split_ix]
        self._mismatch_output_keys = tuple(self._broad_keys + self._split_keys)
        self._mismatch_key_to_index = {
            key: i for i, key in enumerate(self._mismatch_output_keys)
        }
        self._split_out_ix = tuple(
            range(len(self._broad_ix), len(self._broad_ix) + len(self._split_ix))
        )
        self._has_split_mismatches = bool(self._split_ix)

        self._k_needed_mms = len(self._broad_ix) + len(self._split_ix)
        self._required_loci_tuple = tuple(
            loc for loc in mgr.ALL_HLA_LOCI if loc in rules.required_loci
        )
        self._needed_split_loci_tuple = tuple(
            loc for loc in mgr.ALL_HLA_LOCI if loc in rules.needed_split_mismatches
        )
        self._hz_loci_broad = set(rules.needed_broad_mismatches) - set(
            rules.needed_split_mismatches
        )
        self._hz_loci_split = set(rules.needed_split_mismatches)

    @property
    def k_needed_mms(self) -> int:
        return self._k_needed_mms

    @property
    def mismatch_output_keys(self) -> Tuple[str, ...]:
        return self._mismatch_output_keys

    @property
    def mismatch_key_to_index(self) -> Dict[str, int]:
        return self._mismatch_key_to_index

    def prepare_typing_view(self, typing: HLATyping) -> PreparedTyping:
        """ Function to cache a prepared view of the HLA typing. This
            extracts the needed splits and broads.
        """
        cached = typing._prepared_typing
        if cached is not None:
            return cached

        loci_known = {
            mgr.ALL_HLA_LOCI[i]: bool(typing.broads[i])
            for i in range(len(mgr.ALL_HLA_LOCI))
        }
        required_loci_known = all(
            loci_known.get(locus, False)
            for locus in self.rules.required_loci
        )

        match_broads = tuple(typing.broads[i] for i in self._broad_ix)
        match_splits = tuple(typing.splits[i] for i in self._split_ix)
        has_any_unknown_splits = any(typing.has_unknown_splits)

        homozygosity_per_locus = (
            {
                f"hz_{locus}": value
                for locus, value in typing.broads_homozygous.items()
                if locus in self._hz_loci_broad
            }
            | {
                f"hz_{locus}": value
                for locus, value in typing.splits_homozygous.items()
                if locus in self._hz_loci_split
            }
        )
        homozygosity_level = int(
            sum(v for v in homozygosity_per_locus.values() if v is not None)
        )
        fully_homozygous = all(
            v == 1 for v in homozygosity_per_locus.values()
        )
        reduced_hla_match = " ".join(sorted(set().union(*match_broads, *match_splits)))

        prepared = PreparedTyping(
            match_broads=match_broads,
            match_splits=match_splits,
            required_loci_known=required_loci_known,
            required_loci=self._required_loci_tuple,
            needed_split_loci=self._needed_split_loci_tuple,
            broads_by_locus=typing.broads,
            splits_by_locus=typing.splits,
            alleles_by_locus=typing.alleles,
            has_unknown_splits=typing.has_unknown_splits,
            has_unknown_match_splits=any(
                typing.has_unknown_splits[i] for i in self._split_ix
            ),
            has_any_unknown_splits=has_any_unknown_splits,
            homozygosity_per_locus=homozygosity_per_locus,
            homozygosity_level=homozygosity_level,
            fully_homozygous=fully_homozygous,
            reduced_hla_match=reduced_hla_match,
            source_typing=typing,
        )
        typing._prepared_typing = prepared
        return prepared

    def _apply_unknown_split_adjustments(
        self,
        mm_values: List[Optional[int]],
        d: PreparedTyping,
        p: PreparedTyping,
    ) -> None:
        for out_i, i_loc in zip(self._split_out_ix, self._split_ix):
            locus = mgr.ALL_HLA_LOCI[i_loc]
            if not (d.has_unknown_splits[i_loc] or p.has_unknown_splits[i_loc]):
                continue

            common_broads = d.broads_by_locus[i_loc] & p.broads_by_locus[i_loc]
            p_match_hlas = set(p.broads_by_locus[i_loc])
            d_match_hlas = set(d.broads_by_locus[i_loc])

            for common_broad in common_broads:
                if common_broad in self.ontology.splittable_broads[locus]:
                    p_splits_for_broad = (
                        self.ontology.broads_to_splits[locus][common_broad]
                        & p.splits_by_locus[i_loc]
                    )
                    if p_splits_for_broad:
                        d_splits_for_broad = (
                            self.ontology.broads_to_splits[locus][common_broad]
                            & d.splits_by_locus[i_loc]
                        )
                        if d_splits_for_broad:
                            p_match_hlas = (p_match_hlas | p_splits_for_broad) - {
                                common_broad
                            }
                            d_match_hlas = (d_match_hlas | d_splits_for_broad) - {
                                common_broad
                            }

            mm_values[out_i] = _diff_len(d_match_hlas, p_match_hlas) if p_match_hlas else None

    def count_mismatches_prepared(
        self,
        d: PreparedTyping,
        p: PreparedTyping,
        safely: bool = True,
    ) -> Optional[Tuple[Optional[int], ...]]:
        if safely and (not d.required_loci_known or not p.required_loci_known):
            return None

        mm_values: List[Optional[int]] = [0] * self._k_needed_mms

        for out_i, i_loc in enumerate(self._broad_ix):
            mm_values[out_i] = _diff_len(d.broads_by_locus[i_loc], p.broads_by_locus[i_loc])

        for out_i, i_loc in zip(self._split_out_ix, self._split_ix):
            mm_values[out_i] = _diff_len(d.splits_by_locus[i_loc], p.splits_by_locus[i_loc])

        if self._has_split_mismatches and (
            d.has_unknown_match_splits or p.has_unknown_match_splits
        ):
            self._apply_unknown_split_adjustments(mm_values, d, p)
        return tuple(mm_values)

    def count_mismatches(
        self,
        d_hla: HLATyping,
        p_hla: HLATyping,
        safely: bool = True,
    ) -> Optional[Tuple[Optional[int], ...]]:
        d = self.prepare_typing_view(d_hla)
        p = self.prepare_typing_view(p_hla)
        return self.count_mismatches_prepared(d=d, p=p, safely=safely)

    def has_unacceptables(
        self,
        donor_typing: HLATyping,
        normalized_unacceptables: FrozenSet[str],
        exclude_ambiguous_xx: bool = False,
        donor_antigens: Optional[FrozenSet[str]] = None,
    ) -> bool:
        if not normalized_unacceptables:
            return False
        if donor_antigens is None:
            prepared = donor_typing.prepare_antigens()
            donor_antigens = (
                prepared.expanded_antigens_exclude_xx
                if exclude_ambiguous_xx
                else prepared.expanded_antigens_all
            )

        return not normalized_unacceptables.isdisjoint(donor_antigens)
