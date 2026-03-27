from __future__ import annotations

from collections import defaultdict
from copy import copy
from functools import reduce
from operator import or_
from typing import Dict, FrozenSet, List, Optional, Tuple

import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
import simulator.magic_values.column_names as cn
from simulator.code.utils.read_input_files import read_hla_match_table
from simulator.code.utils.utils import DotDict


class HLAOntology:
    """
    ETRL match table. This reads in the match table of antigens
    that are recognized for allocation in Eurotransplant.
    """

    def __init__(self, sim_set: DotDict):
        self.sim_set = sim_set

        hla_match_table = read_hla_match_table(sim_set.PATH_MATCH_TABLE)
        self.orig_table = hla_match_table.copy(deep=True)
        self.orig_table[cn.ALLELE] = (
            self.orig_table[cn.ALLELE]
            .fillna("")
            .astype(str)
            .str.upper()
            .str.strip()
        )
        self.orig_table[cn.SPLIT] = (
            self.orig_table[cn.SPLIT]
            .fillna("")
            .astype(str)
            .str.upper()
            .str.strip()
        )
        self.orig_table[cn.BROAD] = (
            self.orig_table[cn.BROAD]
            .fillna("")
            .astype(str)
            .str.upper()
            .str.strip()
        )
        self.orig_match_codes = {
            code
            for code in (
                set(self.orig_table[cn.ALLELE])
                | set(self.orig_table[cn.SPLIT])
                | set(self.orig_table[cn.BROAD])
            )
            if code
        }

        # Unsplittables (Eurotransplant matches DR17+DR18 as DR3)
        if "UNSPLITTABLES" not in sim_set.keys():
            sim_set.UNSPLITTABLES = {}
        self.UNSPLITTABLES = {
            u: set(v)
            for u, v in sim_set.get("UNSPLITTABLES", {}).items()
        }
        self.UNSPLITTABLE_BROADS = set(self.UNSPLITTABLES.keys())
        self.UNSPLITTABLE_SPLITS = (
            reduce(or_, self.UNSPLITTABLES.values())
            if self.UNSPLITTABLES
            else set()
        )

        self.return_alleles = bool(getattr(sim_set, "RETURN_ALLELES", False))

        # Clean glyphs, to avoid mistypings from being recognized.
        hla_match_table.loc[:, "orig_allele"] = copy(hla_match_table.allele)
        forbidden_characters_regex = getattr(
            es,
            "HLA_FORBIDDEN_CHARACTERS_REGEX",
            es.HLA_FORBIDDEN_CHARACTERS,
        )
        for forbidden_character in forbidden_characters_regex:
            hla_match_table.replace(
                to_replace={
                    "allele": forbidden_character,
                    "split": forbidden_character,
                    "broad": forbidden_character,
                },
                value="",
                inplace=True,
                regex=True,
            )

        self.clean_alle_to_orig = {
            allele: orig
            for allele, orig in zip(
                hla_match_table.allele.str.upper(),
                hla_match_table.orig_allele,
            )
        }
        self.normalized_to_display: Dict[str, str] = {}
        for normalized, original in self.clean_alle_to_orig.items():
            norm_code = str(normalized).upper().strip()
            display_code = str(original).upper().strip()
            if norm_code and display_code and norm_code not in self.normalized_to_display:
                self.normalized_to_display[norm_code] = display_code
        for col in (cn.SPLIT, cn.BROAD):
            for raw_code in self.orig_table[col]:
                normalized = self._normalize_hla_input_string(raw_code)
                if normalized and normalized not in self.normalized_to_display:
                    self.normalized_to_display[normalized] = raw_code

        hla_match_table.allele = hla_match_table.allele.str.upper()
        hla_match_table.split = hla_match_table.split.str.upper()
        hla_match_table.broad = hla_match_table.broad.str.upper()

        self.recognizable_alleles = set(hla_match_table.allele.unique())
        self.recognizable_splits = set(hla_match_table.split.unique())
        self.recognizable_broads = set(hla_match_table.broad.unique())

        self.all_match_codes = (
            self.recognizable_alleles
            .union(self.recognizable_splits)
            .union(self.recognizable_broads)
            .union(mgr.PUBLICS)
        )

        self.alleles_by_type = defaultdict(dict)
        self.splits_by_type = defaultdict(dict)
        self.broads_by_type = defaultdict(dict)
        self.broads_to_splits = defaultdict(lambda: defaultdict(set))
        self.broads_to_alleles = defaultdict(lambda: defaultdict(set))
        self.splits_to_alleles = defaultdict(lambda: defaultdict(set))

        code_to_matchinfo = defaultdict(lambda: defaultdict(dict))

        # Which HLAs to load.
        hlas_to_load = (
            getattr(sim_set, "HLAS_TO_LOAD", None)
            or getattr(sim_set, "hlas_to_load", None)
            or getattr(es, "HLAS_TO_LOAD", None)
            or getattr(es, "hlas_to_load", None)
        )
        if hlas_to_load is None:
            raise AttributeError(
                "No HLA loci configuration found. Configure HLAS_TO_LOAD in sim_set."
            )

        for hla_locus, df_sel in hla_match_table.groupby("type"):
            if hla_locus not in hlas_to_load:
                continue

            self.alleles_by_type[hla_locus] = set(df_sel.allele.unique())
            self.splits_by_type[hla_locus] = set(df_sel.split.unique())
            self.broads_by_type[hla_locus] = set(df_sel.broad.unique())

            for broad, d_alleles in df_sel.groupby("broad"):
                self.broads_to_splits[hla_locus][broad] = set(
                    d_alleles.split.unique()
                )
                self.broads_to_alleles[hla_locus][broad] = set(
                    d_alleles.allele.unique()
                )
                for split, d_alleles_split in d_alleles.groupby("split"):
                    self.splits_to_alleles[hla_locus][split] = set(
                        d_alleles_split.allele.unique()
                    )

                for allele, split, broad in d_alleles.loc[
                    :, [cn.ALLELE, cn.SPLIT, cn.BROAD]
                ].to_records(index=False):
                    code_to_matchinfo[cn.ALLELE].update({allele: (hla_locus, split)})
                    code_to_matchinfo[cn.SPLIT].update(
                        {
                            allele: (hla_locus, split),
                            split: (hla_locus, split),
                        }
                    )
                    code_to_matchinfo[cn.BROAD].update(
                        {
                            allele: (hla_locus, broad),
                            split: (hla_locus, broad),
                            broad: (hla_locus, broad),
                        }
                    )

        # Unsplittables wiring
        for broad, splits in self.UNSPLITTABLES.items():
            locus, _broad = code_to_matchinfo[cn.BROAD][broad]
            for split in splits:
                if split in self.splits_to_alleles[locus]:
                    self.splits_to_alleles[locus][broad] = (
                        self.splits_to_alleles[locus][broad].union(
                            self.splits_to_alleles[locus][split]
                        )
                    )

        for broad, splits in self.UNSPLITTABLES.items():
            locus, _broad = code_to_matchinfo[cn.BROAD][broad]
            self.broads_to_splits[locus][broad] = self.broads_to_splits[locus][
                broad
            ].union((broad,))

        self.codes_to_allele = code_to_matchinfo[cn.ALLELE]
        self.codes_to_broad = code_to_matchinfo[cn.BROAD]
        self.codes_to_split = code_to_matchinfo[cn.SPLIT]

        self.unsplittable_broads = {
            locus: (
                {
                    hla
                    for hla in in_dict.keys()
                    if len(in_dict[hla]) == 1
                }.union(self.UNSPLITTABLE_BROADS)
            )
            for locus, in_dict in self.broads_to_splits.items()
        }
        self.splittable_broads = {
            locus: (
                {
                    hla
                    for hla in in_dict.keys()
                    if len(in_dict[hla]) > 1
                }.difference(self.UNSPLITTABLE_BROADS)
            )
            for locus, in_dict in self.broads_to_splits.items()
        }

        self.missing_splits = defaultdict(int)
        self.unrecognized_antigens: List[str] = []

    def _normalize_hla_input_string(self, hla_string: str) -> str:
        normalized = hla_string.upper()
        for forbidden_char in es.HLA_FORBIDDEN_CHARACTERS:
            normalized = normalized.replace(forbidden_char, "")
        return " ".join(normalized.split())

    def validate_raw_typing_tokens(self, input_string: str) -> Dict[str, object]:
        mapped: List[Tuple[str, str]] = []
        unknown: List[str] = []
        recognized_normalized_tokens: List[str] = []

        seen_mapped = set()
        seen_unknown = set()
        seen_recognized = set()

        for raw_token in input_string.split():
            raw = raw_token.upper().strip()
            if not raw:
                continue

            normalized = self._normalize_hla_input_string(raw)
            if not normalized:
                continue

            if raw in self.orig_match_codes:
                if normalized not in seen_recognized:
                    seen_recognized.add(normalized)
                    recognized_normalized_tokens.append(normalized)
                continue

            if normalized in self.all_match_codes:
                mapped_display = self.normalized_to_display.get(normalized, normalized)
                map_pair = (raw, mapped_display)
                if map_pair not in seen_mapped:
                    seen_mapped.add(map_pair)
                    mapped.append(map_pair)
                if normalized not in seen_recognized:
                    seen_recognized.add(normalized)
                    recognized_normalized_tokens.append(normalized)
                continue

            if raw not in seen_unknown:
                seen_unknown.add(raw)
                unknown.append(raw)

        return {
            "mapped": mapped,
            "unknown": unknown,
            "recognized_normalized_tokens": recognized_normalized_tokens,
        }

    def _antigen_locus(self, antigen: str) -> Optional[str]:
        if antigen in self.codes_to_broad:
            return self.codes_to_broad[antigen][0]
        if antigen in self.codes_to_split:
            return self.codes_to_split[antigen][0]
        if antigen in self.codes_to_allele:
            return self.codes_to_allele[antigen][0]
        return None

    def _antigen_parent(self, antigen: str) -> Optional[str]:
        if antigen in self.recognizable_alleles:
            split = self.codes_to_split.get(antigen, (None, None))[1]
            if split and split != antigen and "?" not in split:
                return split
            broad = self.codes_to_broad.get(antigen, (None, None))[1]
            if broad and broad != antigen and "?" not in broad:
                return broad
            return None

        is_split_only = (
            antigen in self.recognizable_splits
            and antigen not in self.recognizable_broads
        )
        if is_split_only:
            broad = self.codes_to_broad.get(antigen, (None, None))[1]
            if broad and broad != antigen and "?" not in broad:
                return broad
        return None

    def _antigen_children(self, antigen: str) -> FrozenSet[str]:
        locus = self._antigen_locus(antigen)
        if locus is None:
            return frozenset()

        if antigen in self.recognizable_splits:
            return frozenset(
                child
                for child in self.splits_to_alleles[locus].get(antigen, set())
                if "?" not in child
            )

        if antigen in self.recognizable_broads:
            split_children = self.broads_to_splits[locus].get(antigen, set()) - {antigen}
            allele_children = self.broads_to_alleles[locus].get(antigen, set())
            return frozenset(
                child
                for child in split_children.union(allele_children)
                if "?" not in child
            )

        return frozenset()

    def return_all_antigens(
        self,
        input_string: str,
        expand: bool = False,
        exclude_ambiguous_xx: bool = False,
    ) -> List[str]:
        input_string_upp = input_string.upper()
        for code in input_string_upp.split():
            if code not in self.all_match_codes:
                self.missing_splits[code] += 1
                if code not in self.unrecognized_antigens:
                    self.unrecognized_antigens.append(code)

        all_codes = [
            code for code in input_string_upp.split(" ")
            if code in self.all_match_codes
        ]
        if not expand:
            return all_codes

        codes_to_expand = (
            [code for code in all_codes if not code.endswith("XX")]
            if exclude_ambiguous_xx
            else all_codes
        )

        broad_and_split: List[str] = []
        for allele in codes_to_expand:
            if (allele in self.codes_to_broad) and (allele in self.codes_to_split):
                broad_and_split.extend(
                    [self.codes_to_broad[allele][1], self.codes_to_split[allele][1]]
                )
            elif allele in self.codes_to_broad:
                broad_and_split.append(self.codes_to_broad[allele][1])
            elif (allele not in self.codes_to_broad) and (
                allele not in self.codes_to_split
            ):
                continue

        all_codes.extend([x for x in broad_and_split if x is not None])
        return all_codes
