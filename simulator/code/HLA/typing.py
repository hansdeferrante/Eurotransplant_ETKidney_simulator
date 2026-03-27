from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, FrozenSet, List, Optional, Set, Tuple
from warnings import warn

import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
from simulator.code.utils.utils import freeze_set_values

from .ontology import HLAOntology

AntigenTuple = Tuple[Optional[str], Optional[str], Optional[str]]

if TYPE_CHECKING:
    from .api import HLAStatsAPI


@dataclass(frozen=True, slots=True)
class PreparedAntigens:
    expanded_antigens_all: FrozenSet[str]
    expanded_antigens_exclude_xx: FrozenSet[str]
    loci_present_all: FrozenSet[str]
    loci_present_exclude_xx: FrozenSet[str]


class HLATyping:
    """
    HLA typing representation. This depends on the HLA ontology used
    (= ETRL match table)
    """

    __slots__ = (
        "hla_string",
        "ontology",
        "broads",
        "splits",
        "alleles",
        "_unknown_splits",
        "_all_antigens",
        "_broads_homozygous",
        "_splits_homozygous",
        "_reduced_hla_all",
        "_verbose",
        "_prepared_typing",
        "_prepared_antigens",
    )

    def __init__(
        self,
        hla_string: str,
        ontology: Optional[HLAOntology] = None,
        *,
        hla_stats_api: Optional["HLAStatsAPI"] = None,
        verbose: bool = True,
    ):
        if ontology is None and hla_stats_api is not None:
            ontology = hla_stats_api.ontology
        if ontology is None:
            raise TypeError("Provide ontology for HLATyping")
        if not isinstance(hla_string, str):
            raise TypeError(f"hla_string must be str, got {type(hla_string)}")

        self.hla_string = hla_string
        self.ontology = ontology
        self._verbose = verbose

        for fb in es.HLA_FORBIDDEN_CHARACTERS:
            if fb in self.hla_string:
                raise ValueError(
                    f"Forbidden character {fb} in HLA string: {self.hla_string}"
                )

        # Initialize the alleles-splits-broads from the typing.
        self._recompute_structured_hla()

        # caches
        self._unknown_splits = None
        self._all_antigens = None
        self._broads_homozygous = None
        self._splits_homozygous = None
        self._reduced_hla_all = None
        self._prepared_typing = None
        self._prepared_antigens = None

    def _recompute_structured_hla(self) -> None:
        broads, splits, alleles = self.initialize_hlas()
        self.broads = tuple(broads.get(loc, frozenset()) for loc in mgr.ALL_HLA_LOCI)
        self.splits = tuple(splits.get(loc, frozenset()) for loc in mgr.ALL_HLA_LOCI)
        self.alleles = tuple(alleles.get(loc, frozenset()) for loc in mgr.ALL_HLA_LOCI)

    def _invalidate_cached_properties(self, *, clear_prepared_cache: bool = False) -> None:
        self._unknown_splits = None
        self._all_antigens = None
        self._broads_homozygous = None
        self._splits_homozygous = None
        self._reduced_hla_all = None
        if clear_prepared_cache:
            self._prepared_typing = None
            self._prepared_antigens = None

    def set_structured_hla(
        self,
        hla_string: Optional[str],
    ) -> None:
        if hla_string is None:
            raise TypeError("hla_string must be str")
        if not isinstance(hla_string, str):
            raise TypeError(f"hla_string must be str, got {type(hla_string)}")
        self.hla_string = hla_string

        for fb in es.HLA_FORBIDDEN_CHARACTERS:
            if fb in self.hla_string:
                raise ValueError(
                    f"Forbidden character {fb} in HLA string: {self.hla_string}"
                )

        self._recompute_structured_hla()
        self._invalidate_cached_properties(clear_prepared_cache=True)

    def initialize_hlas(
        self,
    ) -> Tuple[
        Dict[str, FrozenSet[str]],
        Dict[str, FrozenSet[str]],
        Dict[str, FrozenSet[str]],
    ]:
        ont = self.ontology
        alleles = defaultdict(set)
        broads = defaultdict(set)
        splits = defaultdict(set)

        for code in self.hla_string.upper().split():
            if code not in ont.all_match_codes:
                ont.missing_splits[code] += 1
                continue

            if code in ont.codes_to_broad:
                locus, broad = ont.codes_to_broad[code]
                broads[locus].add(broad)
                if broad in ont.UNSPLITTABLE_BROADS:
                    splits[locus].add(broad)

            if code in ont.codes_to_split:
                locus, split = ont.codes_to_split[code]
                if split and "?" not in split and split not in ont.UNSPLITTABLE_SPLITS:
                    splits[locus].add(split)

            if ont.return_alleles and code in ont.codes_to_allele:
                locus, _ = ont.codes_to_allele[code]
                alleles[locus].add(code)

        splits = self.filter_split_typings(splits)
        return (
            freeze_set_values(broads),
            freeze_set_values(splits),
            freeze_set_values(alleles),
        )

    def filter_split_typings(
        self,
        hla_typings: Dict[str, Set[str]],
    ) -> Dict[str, Set[str]]:
        filtered_hla_typings = {}

        for locus, antigens in hla_typings.items():
            filtered_antigens = set(antigens)

            for antigen in antigens:
                splits = self.ontology.broads_to_splits[locus][antigen] - {antigen}

                if any(split in antigens for split in splits):
                    filtered_antigens.discard(antigen)

            filtered_hla_typings[locus] = filtered_antigens

        return filtered_hla_typings

    def remove_ambiguous_splits(
        self,
        antigen_sets: List[AntigenTuple],
    ) -> List[AntigenTuple]:
        # Remove broad-only entries if more specific split/allele entries exist.
        broad_dict = defaultdict(list)
        for entry in antigen_sets:
            if (
                entry[1] not in self.ontology.recognizable_broads
                or entry[2] is not None
            ):
                broad_dict[entry[0]].append(entry)

        for entry in antigen_sets:
            if (
                entry[1] in self.ontology.recognizable_broads
                and len(broad_dict[entry[0]]) == 0
                and entry not in broad_dict[entry[0]]
            ):
                broad_dict[entry[0]].append(entry)

        return sum(broad_dict.values(), [])

    def clean_hla_string(self, s: str) -> str:
        if (sf := self.ontology.clean_alle_to_orig.get(s)):
            return sf
        if s:
            s = s.replace("CW", "Cw")
            s = s.replace("DQA", "DQA-")
            s = s.replace("DPA", "DPA-")
            s = s.replace("DP0", "DP-0")
            s = s.replace("DP1", "DP-1")
            s = s.replace("DP4", "DP-4")
        return s

    def return_structured_hla(
        self,
        clean: bool = False,
    ) -> Dict[str, List[AntigenTuple]]:
        ont = self.ontology
        broads = self.broads
        splits = self.splits
        alleles = self.alleles

        antigen_sets: Dict[str, List[AntigenTuple]] = {
            k: [] for k in mgr.ALL_HLA_LOCI
        }
        for iloc, locus in enumerate(mgr.ALL_HLA_LOCI):
            broads_in_loc = broads[iloc]
            for broad in broads_in_loc:
                # If a split is known for this broad, add split and possibly allele.
                for split in splits[iloc]:
                    if ont.codes_to_broad[split][1] == broad:
                        matching_alleles = ont.splits_to_alleles[locus][split] & alleles[iloc]
                        if matching_alleles:
                            for allele in matching_alleles:
                                corr_split = ont.codes_to_split[allele][1]
                                if corr_split == split:
                                    antigen_sets[locus].append((broad, split, allele))
                                elif corr_split in ont.UNSPLITTABLE_SPLITS:
                                    if ont.codes_to_broad[corr_split][1] == broad:
                                        antigen_sets[locus].append((broad, split, allele))
                        else:
                            antigen_sets[locus].append((broad, split, None))

                # If no split is known for this broad, include broad/allele fallback.
                if not (splits[iloc] & ont.broads_to_splits[locus][broad]):
                    matching_alleles = (
                        alleles[iloc] & ont.broads_to_alleles[locus][broad]
                    )
                    if (not matching_alleles) and (broad not in ont.unsplittable_broads[locus]):
                        antigen_sets[locus].append((broad, None, None))
                    else:
                        for allele in matching_alleles:
                            antigen_sets[locus].append((broad, None, allele))

        # Check whether all antigens were recognized.
        recognized_antigens = {
            locus: {
                antigen
                for antigens in list_of_tuples
                for antigen in antigens
                if antigen is not None
            }
            for locus, list_of_tuples in antigen_sets.items()
        }
        for iloc, (locus, antigens) in enumerate(recognized_antigens.items()):
            missing_antigens = broads[iloc].difference(antigens).union(
                splits[iloc].difference(antigens)
            )
            if alleles:
                missing_antigens = missing_antigens.union(
                    alleles[iloc].difference(antigens)
                )

            if missing_antigens:
                for antigen in missing_antigens:
                    if "XX" in antigen:
                        if antigen not in ont.unrecognized_antigens:
                            ont.unrecognized_antigens.append(antigen)
                            warn(
                                f"Antigen {antigen} is not recognized. "
                                "Likely ambiguous (XX)"
                            )
                    else:
                        raise Exception(
                            "The following antigens were not recognized: "
                            f"{missing_antigens}. Structured antigens: "
                            f"{antigen_sets}. For string: {self.hla_string}"
                        )

        # Keep max 2 entries per locus.
        for locus, ant_set in antigen_sets.items():
            if len(ant_set) > 2:
                if self._verbose:
                    warn(
                        f"Found for locus {locus} the following antigen sets:\n"
                        f" {ant_set}\nInput string: {self.hla_string}. "
                        "Trying to remove ambiguous splits."
                    )
                antigen_sets[locus] = self.remove_ambiguous_splits(ant_set)
                if len(antigen_sets[locus]) > 2:
                    raise Exception(
                        "Ambiguous splits remain: "
                        f"{antigen_sets[locus]}\n"
                        f"Original string: {self.hla_string}"
                    )
            elif (
                (len(ant_set) == 2)
                and (len(set(a[0] for a in ant_set)) == 1)
                and (len(set(a[2] for a in ant_set)) != 2)
            ):
                antigen_sets[locus] = self.remove_ambiguous_splits(ant_set)

        sorted_antigens: Dict[str, List[AntigenTuple]] = {
            locus: sorted(antigens)
            for locus, antigens in antigen_sets.items()
        }

        if clean:
            sorted_antigens = {
                locus: [
                    tuple(self.clean_hla_string(k) for k in tp)  # type: ignore[misc]
                    for tp in tps
                ]
                for locus, tps in sorted_antigens.items()
            }
        return sorted_antigens

    def return_structured_hla_wide(self, clean: bool = False) -> Dict[str, Optional[str]]:
        sorted_antigens = self.return_structured_hla(clean=clean)
        return {
            (
                f"{locus}_{i + 1}_"
                f"{'broad' if k == 0 else 'split' if k == 1 else 'allele'}"
            ): ant
            for locus, antigen_sets in sorted_antigens.items()
            for i, antigen_set in enumerate(antigen_sets)
            for k, ant in enumerate(antigen_set)
        }

    def return_clean_hla_string(self) -> str:
        clean_hla = self.return_structured_hla(clean=True)
        return " ".join(
            dict.fromkeys(
                item
                for value in clean_hla.values()
                for triplet in value
                for item in triplet
                if item is not None
            )
        )

    @property
    def has_unknown_splits(self) -> Tuple[bool, ...]:
        if self._unknown_splits is None:
            self._unknown_splits = tuple(
                len(s) < len(b)
                for s, b in zip(self.splits, self.broads)
            )
        return self._unknown_splits

    @property
    def broads_homozygous(self) -> Dict[str, Optional[int]]:
        if self._broads_homozygous is None:
            self._broads_homozygous = {
                mgr.ALL_HLA_LOCI[i]: (int(len(antigens) == 1) if antigens else None)
                for i, antigens in enumerate(self.broads)
            }
        return self._broads_homozygous

    @property
    def splits_homozygous(self) -> Dict[str, Optional[int]]:
        if self._splits_homozygous is None:
            self._splits_homozygous = {
                mgr.ALL_HLA_LOCI[i]: (int(len(antigens) == 1) if antigens else None)
                for i, antigens in enumerate(self.splits)
            }
        return self._splits_homozygous

    @property
    def all_antigens(self) -> Set[str]:
        if self._all_antigens is None:
            self._all_antigens = set().union(*self.broads)
            self._all_antigens |= set().union(*self.splits)
        return self._all_antigens

    def _loci_present_for_antigens(
        self,
        donor_antigens: FrozenSet[str],
    ) -> FrozenSet[str]:
        return frozenset(
            locus
            for locus in (self.ontology._antigen_locus(code) for code in donor_antigens)
            if locus is not None
        )

    def prepare_antigens(self) -> PreparedAntigens:
        cached = self._prepared_antigens
        if cached is not None:
            return cached

        normalized_hla = self.ontology._normalize_hla_input_string(self.hla_string)
        expanded_antigens_all = frozenset(
            self.ontology.return_all_antigens(
                normalized_hla,
                expand=True,
                exclude_ambiguous_xx=False,
            )
        )
        expanded_antigens_exclude_xx = frozenset(
            self.ontology.return_all_antigens(
                normalized_hla,
                expand=True,
                exclude_ambiguous_xx=True,
            )
        )
        prepared = PreparedAntigens(
            expanded_antigens_all=expanded_antigens_all,
            expanded_antigens_exclude_xx=expanded_antigens_exclude_xx,
            loci_present_all=self._loci_present_for_antigens(expanded_antigens_all),
            loci_present_exclude_xx=self._loci_present_for_antigens(
                expanded_antigens_exclude_xx
            ),
        )
        self._prepared_antigens = prepared
        return prepared

    @property
    def reduced_hla(self) -> str:
        # Policy-free reduced representation: all known broad/split antigens.
        if self._reduced_hla_all is None:
            self._reduced_hla_all = " ".join(sorted(self.all_antigens))
        return self._reduced_hla_all

    @property
    def reduced_hla_all(self) -> str:
        return self.reduced_hla

    def __str__(self) -> str:
        ordered_hla = self.return_structured_hla()
        parts = []
        for locus in mgr.ALL_HLA_LOCI:
            hla_entries = ordered_hla.get(locus, [])
            if hla_entries:
                formatted = []
                for entry in hla_entries:
                    formatted_entry = "-".join(filter(None, entry))
                    formatted.append(f"({formatted_entry})")
                parts.append(f"{locus}: {', '.join(formatted)}")
        return ", ".join(parts)

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        for slot in self.__slots__:
            value = getattr(self, slot)
            if slot == "ontology":
                setattr(result, slot, value)
            else:
                setattr(result, slot, deepcopy(value, memo))
        return result


class AmbiguousHLATyping(HLATyping):
    """
    Class implementing an HLA typing, allowing for ambiguity on the alleles.
    """

    __slots__ = HLATyping.__slots__ + ("manual_hla_string",)

    def __init__(
        self,
        hla_string: str,
        ontology: Optional[HLAOntology] = None,
        *,
        hla_stats_api: Optional["HLAStatsAPI"] = None,
        manual_hla_string: Optional[str] = None,
        verbose: bool = True,
    ):
        self.manual_hla_string = manual_hla_string
        super().__init__(
            hla_string=hla_string,
            ontology=ontology,
            hla_stats_api=hla_stats_api,
            verbose=verbose,
        )
        self._enforce_no_split_or_broad_ambiguity()

    def _manual_allowed_split_and_broad(
        self,
    ) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
        allowed_splits = defaultdict(set)
        allowed_broads = defaultdict(set)
        if not self.manual_hla_string:
            return allowed_splits, allowed_broads

        for code in self.manual_hla_string.split():
            if code in self.ontology.codes_to_split:
                locus, split = self.ontology.codes_to_split[code]
                if split and "?" not in split:
                    allowed_splits[locus].add(split)
            if code in self.ontology.codes_to_broad:
                locus, broad = self.ontology.codes_to_broad[code]
                if broad and "?" not in broad:
                    allowed_broads[locus].add(broad)

        return allowed_splits, allowed_broads

    def _filter_alleles_using_manual_phenotype(
        self,
        loci_with_too_many_splits: Set[str],
    ) -> str:
        if not loci_with_too_many_splits or not self.manual_hla_string:
            return self.hla_string

        allowed_splits, allowed_broads = self._manual_allowed_split_and_broad()
        filtered_tokens = []
        for token in self.hla_string.split():
            if token not in self.ontology.codes_to_allele:
                filtered_tokens.append(token)
                continue

            locus, split = self.ontology.codes_to_split[token]
            if locus not in loci_with_too_many_splits:
                filtered_tokens.append(token)
                continue

            broad = self.ontology.codes_to_broad[token][1]
            if allowed_splits.get(locus):
                if split in allowed_splits[locus]:
                    filtered_tokens.append(token)
            elif allowed_broads.get(locus):
                if broad in allowed_broads[locus]:
                    filtered_tokens.append(token)
            else:
                filtered_tokens.append(token)

        return " ".join(filtered_tokens)

    def _loci_with_n_splits_above(self, threshold: int) -> Set[str]:
        loci = set()
        for iloc, locus in enumerate(mgr.ALL_HLA_LOCI):
            if len(self.splits[iloc]) > threshold:
                loci.add(locus)
        return loci

    def _loci_with_n_broads_above(self, threshold: int) -> Set[str]:
        loci = set()
        for iloc, locus in enumerate(mgr.ALL_HLA_LOCI):
            if len(self.broads[iloc]) > threshold:
                loci.add(locus)
        return loci

    def _enforce_no_split_or_broad_ambiguity(self) -> None:
        loci_with_too_many_splits = self._loci_with_n_splits_above(2)
        if loci_with_too_many_splits and self.manual_hla_string:
            filtered_hla = self._filter_alleles_using_manual_phenotype(
                loci_with_too_many_splits=loci_with_too_many_splits
            )
            if filtered_hla != self.hla_string:
                self.set_structured_hla(
                    hla_string=filtered_hla,
                )

        split_ambiguous_loci = self._loci_with_n_splits_above(2)
        broad_ambiguous_loci = self._loci_with_n_broads_above(2)
        if split_ambiguous_loci or broad_ambiguous_loci:
            raise ValueError(
                "AmbiguousHLATyping allows ambiguity only on allele level. "
                f"Split-ambiguous loci: {sorted(split_ambiguous_loci)}. "
                f"Broad-ambiguous loci: {sorted(broad_ambiguous_loci)}. "
                f"HLA: {self.hla_string}"
            )
