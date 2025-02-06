#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:33:44 2022

@author: H.C. de Ferrante
"""

import numpy as np
from warnings import warn
from typing import (
    List, Dict, Tuple, Optional, Set,
    FrozenSet, Generator, Union, Any
)
from copy import copy
from collections import defaultdict
from copy import deepcopy

from functools import reduce
from operator import or_

import simulator.magic_values.etkidney_simulator_settings as es
import simulator.magic_values.magic_values_rules as mgr
import simulator.magic_values.column_names as cn
import simulator.code.utils.read_input_files as rdr

from simulator.code.utils.read_input_files import read_hla_match_table
from simulator.code.utils.utils import (
    DotDict, freeze_set_values
)


class HLAProfile:
    """Class which implements an HLA profile.

    Attributes
    ----------
    hla_string : str
        Original HLA string.
    alleles : Tuple[Set[str]]
        Allele codes per locus.
    splits : Tuple[Set[str]]
        Split HLA codes per locus.
    broads : Tuple[Set[str]]
        Broad HLA codes per locus.
    hla_system : 'HLASystem'
        Reference to the HLA system used for code lookups.

    Methods
    -------
    set_structured_hla(hla_string: str)
        Processes the input HLA string.
    return_structured_hla(input_string: str)
        Returns structured HLA data for broads, splits, and alleles.
    """
    __slots__ = (
        'hla_string', 'hla_system', 'broads', 'splits', 'alleles',
        '_unknown_splits', '_reduced_hla', '_broads_homozygous',
        '_splits_homozygous', '_homozygosity_per_locus', '_fully_homozygous',
        'hla_known', 'required_loci_known', 'required_loci',
        'needed_split_mismatches', 'match_broads', 'match_splits',
        '_all_antigens', '_verbose'
    )

    def __init__(
        self, hla_string: str,
        hla_system: 'HLASystem.HLASystem', verbose: bool = True
    ):
        self.hla_string = hla_string
        self.hla_system = hla_system

        # Initialize the attributes with empty sets for each locus
        self.set_structured_hla(self.hla_string, hla_system=hla_system)

        # Initialize for caching.
        self._all_antigens = None
        self._unknown_splits = None
        self._reduced_hla = None
        self._broads_homozygous = None
        self._splits_homozygous = None
        self._homozygosity_per_locus = None
        self._fully_homozygous = None
        self._verbose = verbose

    def set_structured_hla(
        self, hla_string: Optional[str],
        hla_system: 'HLASystem.HLASystem'
    ):
        if not isinstance(hla_string, str):
            raise Exception(
                f'input string should be a string, not {hla_string}'
            )
        for fb in es.HLA_FORBIDDEN_CHARACTERS:
            if fb in self.hla_string:
                raise Exception(
                    f'Forbidden character {fb} in donor '
                    f'HLA-input string!\n({self.hla_string})'
                )

        # Format the HLA string
        broads, splits, alleles = self.initialize_hlas()
        self.broads = tuple(broads.get(loc, set()) for loc in mgr.ALL_HLA_LOCI)
        self.splits = tuple(splits.get(loc, set()) for loc in mgr.ALL_HLA_LOCI)
        self.alleles = tuple(
            alleles.get(loc, set()) for loc in mgr.ALL_HLA_LOCI
        )

        # Initialize whether any HLA is known per locus
        self.hla_known = {
            locus: any(t) for t, locus
            in zip(self.broads, mgr.ALL_HLA_LOCI)
        }
        self.required_loci_known = all(
            self.hla_known[locus] for locus
            in hla_system.required_loci
        )

        # Format also the match determinantss
        self.required_loci = self.hla_system.required_loci
        self.needed_split_mismatches = self.hla_system.needed_split_mismatches
        self.match_broads = tuple(
            broads.get(locus, set()) for locus in self.required_loci
        )
        self.match_splits = tuple(
            splits.get(locus, set()) for locus in self.needed_split_mismatches
        )

    def initialize_hlas(
        self,
    ) -> Tuple[
        Dict[str, FrozenSet[str]],
        Dict[str, FrozenSet[str]],
        Optional[Dict[str, FrozenSet[str]]]
    ]:
        """
        Initialize a candidates HLA typings, based on an input HLA string
        """
        hla_system = self.hla_system

        alleles = defaultdict(set)
        broads = defaultdict(set)
        splits = defaultdict(set)
        for code in self.hla_string.upper().split():
            if code not in hla_system.all_match_codes:
                hla_system.missing_splits[code] += 1
            else:
                if code in hla_system.codes_to_broad:
                    locus, broad = hla_system.codes_to_broad[code]
                    broads[locus].add(broad)
                    if broad in self.hla_system.UNSPLITTABLE_BROADS:
                        splits[locus].add(broad)
                if code in hla_system.codes_to_split:
                    locus, split = hla_system.codes_to_split[code]
                    if '?' in split:
                        pass
                    elif split not in self.hla_system.UNSPLITTABLE_SPLITS:
                        splits[locus].add(split)
                if (
                    hla_system.return_alleles and
                    code in hla_system.codes_to_allele
                ):
                    locus, split = hla_system.codes_to_allele[code]
                    alleles[locus].add(code)

        splits = self.filter_split_typings(splits)

        broads = freeze_set_values(broads)
        splits = freeze_set_values(splits)
        alleles = freeze_set_values(alleles)

        return broads, splits, alleles

    def filter_split_typings(
        self,
        hla_typings: Dict[str, Set[str]]
    ) -> Dict[str, Set[str]]:
        """
        Filters HLA typings to the lowest known level,
        removing broad antigens if their splits are present.

        :param hla_typings: Dictionary of HLA typings per locus
        :param broads_to_splits: Dictionary mapping loci to broad antigens
                                 to their split levels
        :return: Filtered HLA typings
        """
        filtered_hla_typings = {}

        for locus, antigens in hla_typings.items():
            filtered_antigens = set(antigens)

            for antigen in antigens:
                # Remove the antigen itself from the mapping
                splits = (
                    self.hla_system.broads_to_splits[locus][antigen] -
                    {antigen}
                )

                # If any split is present, remove the broad antigen
                if any(split in antigens for split in splits):
                    filtered_antigens.discard(antigen)

            filtered_hla_typings[locus] = filtered_antigens

        return filtered_hla_typings

    def remove_ambiguous_splits(
        self, antigen_sets: List[Tuple[Optional[str]]]
    ) -> List[Tuple[Optional[str]]]:
        """This code tries to remove ambiguous splits from an HLA typing."""

        # If the 2nd entry (the split level) is recognizable as a broad, it is
        # likely a case where the broad has splits, but some alleles are
        # not recognized as splits. For instance, B*40:05 is not
        # recognized as a split, but B40 is splittable in B60 and B61.
        # If so, remove (B40, B40, None).

        # Dictionary to hold the entries by their broad category
        broad_dict = defaultdict(list)
        for entry in antigen_sets:
            if (
                entry[1] not in self.hla_system.recognizable_broads or
                entry[2] is not None
            ):
                broad_dict[entry[0]].append(entry)

        for entry in antigen_sets:
            if (
                entry[1] in self.hla_system.recognizable_broads and
                (
                    (len(broad_dict[entry[0]]) == 0)
                ) and entry not in broad_dict[entry[0]]
            ):
                broad_dict[entry[0]].append(entry)

        return sum(broad_dict.values(), [])

    def clean_hla_string(self, s: str) -> str:
        if (sf := self.hla_system.clean_alle_to_orig.get(s)):
            return sf
        if s:
            s = s.replace('CW', 'Cw')
            s = s.replace('DQA', 'DQA-')
            s = s.replace('DPA', 'DPA-')
            s = s.replace('DP0', 'DP-0')
            s = s.replace('DP1', 'DP-1')
            s = s.replace('DP4', 'DP-4')
        return s

    def return_structured_hla(
        self,
        clean=False
    ) -> Dict[str, List[Tuple[Optional[str]]]]:
        # Take a string as an input, recognize the nested structure
        # of antigens in splits in broads, and return this structured
        # HLA as a dictionary. Values can contain 2 elements, i.e.
        # the maternal and paternal HLA phenotype.
        hla_system = self.hla_system
        broads = self.broads
        splits = self.splits
        alleles = self.alleles

        antigen_sets = {k: [] for k in mgr.ALL_HLA_LOCI}
        for iloc, locus in enumerate(mgr.ALL_HLA_LOCI):
            broads_in_loc = broads[iloc]
            locus = mgr.ALL_HLA_LOCI[iloc]
            for broad in broads_in_loc:
                # In case a split is known for the broad,
                # add split and check whether allele is present
                for split in splits[iloc]:
                    if hla_system.codes_to_broad[split][1] == broad:
                        if (
                            matching_alleles := hla_system.splits_to_alleles[
                                locus
                            ][split] & alleles[iloc]
                        ):
                            for allele in matching_alleles:
                                if (corr_split := hla_system.codes_to_split[
                                    allele
                                ][1]) == split:
                                    antigen_sets[locus].append(
                                        (broad, split, allele)
                                    )
                                elif corr_split in self.hla_system.UNSPLITTABLE_SPLITS:
                                    if hla_system.codes_to_broad[
                                        corr_split
                                    ][1] == broad:
                                        antigen_sets[locus].append(
                                            (broad, split, allele)
                                        )
                        else:
                            antigen_sets[locus].append((broad, split, None))

                # In case no split is known for the broad,
                # still check whether allele is present
                if not (
                    splits[iloc] & hla_system.broads_to_splits[locus][broad]
                ):
                    if not (
                        (
                            matching_alleles := (
                                alleles[iloc] &
                                hla_system.broads_to_alleles[locus][broad]
                            )
                        )
                    ) and (
                        broad not in hla_system.unsplittable_broads[locus]
                    ):
                        antigen_sets[locus].append((broad, None, None))
                    else:
                        for allele in matching_alleles:
                            antigen_sets[locus].append((broad, None, allele))

        # Check whether all antigens were recognized.
        recognized_antigens = {
            locus: set(
                antigen
                for antigens in list_of_tuples
                for antigen in antigens
                if antigen is not None
            ) for locus, list_of_tuples in antigen_sets.items()
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
                    if 'XX' in antigen:
                        if antigen not in hla_system.unrecognized_antigens:
                            hla_system.unrecognized_antigens.append(antigen)
                            warn(
                                f'Antigen {antigen} is not recognized.'
                                f' Likely ambiguous (XX)'
                            )
                    else:
                        raise Exception(
                            f'The following antigens were not recognized: '
                            f'{missing_antigens}. Structured antigens: '
                            f'{antigen_sets}. For string: {self.hla_string}')

        # Check whether at most 2 typings are known per locus
        for locus, ant_set in antigen_sets.items():
            if len(ant_set) > 2:
                if self._verbose:
                    warn(
                        f'Found for locus {locus} the following antigen sets:'
                        f'\n {ant_set}\nInput string: {self.hla_string}. '
                        f'Trying to remove ambiguous splits.'
                    )
                    print(self.remove_ambiguous_splits(ant_set))
                antigen_sets[locus] = self.remove_ambiguous_splits(ant_set)
                if len(antigen_sets[locus]) > 2:
                    raise Exception(
                        f'Ambiguous splits remain: '
                        f'{antigen_sets[locus]}\n'
                        f'Original string: {self.hla_string}'
                    )
            # In case alleles are not known, remove ambiguous
            # split typings
            elif (
                (len(ant_set) == 2) &
                (len(set(a[0] for a in ant_set)) == 1) &
                (len(set(a[2] for a in ant_set)) != 2)
            ):
                antigen_sets[locus] = self.remove_ambiguous_splits(ant_set)

        # Sort & convert into a dictionary structure
        sorted_antigens = {
            locus: sorted(antigens) for locus, antigens in antigen_sets.items()
        }

        if clean:
            sorted_antigens = {
                locus: [
                    tuple(
                        self.clean_hla_string(k) for k in tp
                    ) for tp in tps
                ]
                for locus, tps in sorted_antigens.items()
            }
        return sorted_antigens

    def __deepcopy__(self, memo):
        # Create a new instance of the class
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        # Copy all attributes, except for hla_system and bal_system
        for k, v in self.__dict__.items():
            if k in {'hla_system'}:
                setattr(result, k, v)  # Shallow copy
            else:
                setattr(result, k, deepcopy(v, memo))  # Deep copy

        return result

    @property
    def homozygosity_per_locus(self, hz_prefix: str = 'hz_') -> Dict[str, int]:
        # Return a dictionary of homozygosity.
        # This dictionary only considers needed split
        # and needed broad mismatches.
        if self._homozygosity_per_locus is not None:
            return self._homozygosity_per_locus

        self._homozygosity_per_locus = (
            {
                f'{hz_prefix}{k}': v
                for k, v in self.broads_homozygous.items()
                if k in self.hla_system.hz_loci_broad
            } | {
                f'{hz_prefix}{k}': v
                for k, v in self.splits_homozygous.items()
                if k in self.hla_system.hz_loci_split
            }
        )
        return self._homozygosity_per_locus

    @property
    def homozygosity_level(self) -> int:
        return int(sum(self.homozygosity_per_locus.values()))

    @property
    def fully_homozygous(
        self, hz_prefix: str = 'hz_'
    ) -> bool:
        # Determine whether an HLA is fully homozygous
        if self._fully_homozygous is not None:
            return self._fully_homozygous
        else:
            self._fully_homozygous = all(
                lc == 1 for lc in self.homozygosity_per_locus.values()
            )
        return self._fully_homozygous

    @property
    def has_unknown_splits(self) -> Tuple[bool]:
        if self._unknown_splits is None:
            self._unknown_splits = tuple(
                len(s) < len(b) for s, b in zip(self.splits, self.broads)
            )

        return self._unknown_splits

    @property
    def reduced_hla(self):
        if self._reduced_hla is None:
            sel_hla = self.match_broads + self.match_splits
            self._reduced_hla = ' '.join(
                sorted(set().union(*sel_hla))
            )

        return self._reduced_hla

    @property
    def broads_homozygous(self):
        if self._broads_homozygous is None:
            self._broads_homozygous = {
                mgr.ALL_HLA_LOCI[i]: int(len(antigens) == 1)
                if antigens else None
                for i, antigens in enumerate(self.broads)
            }
        return self._broads_homozygous

    @property
    def splits_homozygous(self):
        if self._splits_homozygous is None:
            self._splits_homozygous = {
                mgr.ALL_HLA_LOCI[i]: int(len(antigens) == 1)
                if antigens else None
                for i, antigens in enumerate(self.splits)
            }
        return self._splits_homozygous

    @property
    def all_antigens(self):
        if self._all_antigens is None:
            self._all_antigens = set().union(*self.broads)
            self._all_antigens |= set().union(*self.splits)

        return self._all_antigens

    def __str__(self) -> str:
        # Get the structured HLA data
        ordered_hla = self.return_structured_hla()

        # Prepare the output string
        result = []

        for locus in mgr.ALL_HLA_LOCI:
            locus_key = locus.lower()
            hla_entries = ordered_hla.get(locus_key, [])

            if hla_entries:
                formatted_entries = []
                for entry in hla_entries:
                    # Format each tuple in the list as a string
                    # "(A11-A11-A02XX)"
                    formatted_entry = '-'.join(filter(None, entry))
                    formatted_entries.append(f'({formatted_entry})')

                result.append(f'{locus_key}: {", ".join(formatted_entries)}')

        return ', '.join(result)


class HLASystem:
    """Class which implements the HLA system
    ...

    Attributes   #noqa
    ----------
    alleles: Set[str]
        allele codes
    splits: Set[str]
        split HLA codes
    broads: Set[str]
        broad HLA codes
    self.alleles_by_type: Dict[str, Set[str]]
        alleles organized by locus
    self.splits_by_type: Dict[str, Set[str]]
        splits organized by locus
    self.broads_by_type: Dict[str, Set[str]]
        broads organized by locus
    self.broads_to_splits: Dict[str, Dict[str, Set[str]]
        per locus mapping of broads to set of splits.
    self.codes_to_broad: Dict[str, Dict[str, str]]
        per locus mapping of HLA-code to corresponding broad
    self.codes_to_split: Dict[str, Dict[str, str]]
        per locus mapping of HLA-code to corresponding split
    self.unsplittable_broads: Dict[str, Set[str]]
        unsplittable broads organized per locus
    self.loci_zero_mismatch: Set[str]
        loci for full house match
    self.allele_frequencies_split: Dict[str, float]
        allele frequences at split level
    self.allele_frequencies_broad: Dict[str, float]
        allele frequencies at broad level

    Methods
    -------
    return_all_antigens(input_string: str) -> Generator[str, None, None]
        Returns a set of all antigens recognized against ET match tables

    """

    def __init__(
        self,
        sim_set: DotDict
    ):

        self.sim_set = sim_set
        hla_match_table = read_hla_match_table(
            sim_set.PATH_MATCH_TABLE
        )

        # Unsplittables
        if 'UNSPLITTABLES' not in sim_set.keys():
            input('No unsplittables. Are you sure?')
        self.UNSPLITTABLES = {
            u: set(v)
            for u, v in
            sim_set.get('UNSPLITTABLES', {}).items()
        }
        self.UNSPLITTABLE_BROADS = set(self.UNSPLITTABLES.keys())
        if self.UNSPLITTABLES:
            self.UNSPLITTABLE_SPLITS = reduce(or_, self.UNSPLITTABLES.values())
        else:
            self.UNSPLITTABLE_SPLITS = {}
        hla_match_table.loc[:, 'orig_allele'] = copy(hla_match_table.allele)

        if self.sim_set.RETURN_ALLELES:
            self.return_alleles = True
        else:
            self.return_alleles = False

        for forbidden_character in es.HLA_FORBIDDEN_CHARACTERS:
            hla_match_table.replace(
                to_replace={
                    'allele': forbidden_character,
                    'split': forbidden_character,
                    'broad': forbidden_character
                },
                value='',
                inplace=True,
                regex=True
            )

        self.clean_alle_to_orig = {allele: orig for allele, orig in zip(
            hla_match_table.allele.str.upper(),
            hla_match_table.orig_allele
        )
        }

        hla_match_table.allele = hla_match_table.allele.str.upper()
        hla_match_table.split = hla_match_table.split.str.upper()
        hla_match_table.broad = hla_match_table.broad.str.upper()

        self.recognizable_alleles = set(hla_match_table.allele.unique())
        self.recognizable_splits = set(hla_match_table.split.unique())
        self.recognizable_broads = set(hla_match_table.broad.unique())

        self.all_match_codes = (
            self.recognizable_alleles.union(self.recognizable_splits).union(
                self.recognizable_broads
            ).union(mgr.PUBLICS)
        )
        self.alleles_by_type = defaultdict(dict)
        self.splits_by_type = defaultdict(dict)
        self.broads_by_type = defaultdict(dict)
        self.broads_to_splits = defaultdict(lambda: defaultdict(set))
        self.broads_to_alleles = defaultdict(lambda: defaultdict(set))
        self.splits_to_alleles = defaultdict(lambda: defaultdict(set))

        self.loci_zero_mismatch = sim_set.LOCI_ZERO_MISMATCH

        code_to_matchinfo = defaultdict(
            lambda: defaultdict(dict)
        )

        if sim_set.HLAS_TO_LOAD:
            hlas_to_load = sim_set.HLAS_TO_LOAD
        else:
            hlas_to_load = es.hlas_to_load

        # Load from match table the alleles
        for hla_locus, df_sel in hla_match_table.groupby('type'):
            if hla_locus not in hlas_to_load:
                continue
            self.alleles_by_type[hla_locus] = set(df_sel.allele.unique())
            self.splits_by_type[hla_locus] = set(df_sel.split.unique())
            self.broads_by_type[hla_locus] = set(df_sel.broad.unique())
            for broad, d_alleles in df_sel.groupby('broad'):
                self.broads_to_splits[hla_locus][
                    broad
                ] = set(d_alleles.split.unique())
                self.broads_to_alleles[hla_locus][
                    broad
                ] = set(d_alleles.allele.unique())
                for split, d_alleles_split in d_alleles.groupby('split'):
                    self.splits_to_alleles[hla_locus][
                        split
                    ] = set(d_alleles_split.allele.unique())

                for allele, split, broad in d_alleles.loc[
                    :,
                    [cn.ALLELE, cn.SPLIT, cn.BROAD]
                ].to_records(index=False):
                    code_to_matchinfo[cn.ALLELE].update(
                        {
                            allele: (hla_locus, split)
                        }
                    )
                    code_to_matchinfo[cn.SPLIT].update(
                        {
                            allele: (hla_locus, split),
                            split: (hla_locus, split)
                        }
                    )
                    code_to_matchinfo[cn.BROAD].update(
                        {
                            allele: (hla_locus, broad),
                            split: (hla_locus, broad),
                            broad: (hla_locus, broad)
                        }
                    )

        # Add unsplittable broads to splits_to_alleles dictionary.
        for broad, splits in self.UNSPLITTABLES.items():
            locus, broad = code_to_matchinfo[cn.BROAD][broad]
            for split in splits:
                if split in self.splits_to_alleles[locus]:
                    self.splits_to_alleles[locus][
                        broad
                    ] = self.splits_to_alleles[locus][broad].union(
                        self.splits_to_alleles[locus][split]
                    )
        # Replace unsplittables from broads_to_splits dictionary
        for broad, splits in self.UNSPLITTABLES.items():
            locus, broad = code_to_matchinfo[cn.BROAD][broad]
            self.broads_to_splits[locus][
                broad
            ] = self.broads_to_splits[locus][broad].union((broad,))

        self.codes_to_allele = code_to_matchinfo[cn.ALLELE]
        self.codes_to_broad = code_to_matchinfo[cn.BROAD]
        self.codes_to_split = code_to_matchinfo[cn.SPLIT]

        self.unsplittable_broads = {
            locus: (
                {hla for hla in in_dict.keys() if len(in_dict[hla]) == 1}
                .union(self.UNSPLITTABLE_BROADS)
            ) for locus, in_dict in self.broads_to_splits.items()
        }
        self.splittable_broads = {
            locus: (
                {hla for hla in in_dict.keys() if len(in_dict[hla]) > 1}
                .difference(self.UNSPLITTABLE_BROADS)
            ) for locus, in_dict in self.broads_to_splits.items()
        }

        self.missing_splits = defaultdict(int)
        self.unrecognized_antigens = list()

        self.needed_broad_mismatches = (
            set(sim_set.NEEDED_BROAD_MISMATCHES)
            if sim_set.NEEDED_BROAD_MISMATCHES is not None
            else {mgr.HLA_A, mgr.HLA_B}
        )
        self.needed_split_mismatches = (
            set(sim_set.NEEDED_SPLIT_MISMATCHES)
            if sim_set.NEEDED_SPLIT_MISMATCHES is not None
            else {mgr.HLA_DR}
        )
        self.required_loci = (
            self.needed_broad_mismatches |
            self.needed_split_mismatches
        )
        if (
            common_loci := (
                self.needed_split_mismatches &
                self.needed_broad_mismatches
            )
        ):
            warn(
                f'The following loci is matched at both broad '
                f'and split level: {common_loci}. The split level '
                f'will be used for matching and homozygosity'
            )
        self.hz_loci_broad = (
            self.needed_broad_mismatches - self.needed_split_mismatches
        )
        self.hz_loci_split = self.needed_split_mismatches

        self.k_needed_mms = (
            len(self.needed_split_mismatches) +
            len(self.needed_broad_mismatches)
        )
        self._donor_pool_hlas = None

    def lookup_alleles(self, input_generator, codes_to_broad, codes_to_split):
        for allele in input_generator:
            if (allele in codes_to_broad) & (allele in codes_to_split):
                yield (codes_to_broad[allele][1], codes_to_split[allele][1])
            elif (allele in codes_to_broad):
                yield (codes_to_broad[allele][1]), None
            elif (allele not in codes_to_broad) and (allele not in codes_to_split):
                continue
            else:
                print('Not known!')
                breakpoint()

    def return_all_antigens(
        self, input_string: str, expand=False
    ) -> List[str]:
        """
            Returns all antigens as a generator.
            if Expand is set to True, any allele
            will be expanded to the complete typing
            on broads and splits.
        """
        input_string_upp = input_string.upper()
        for code in input_string_upp.split():
            if code not in self.all_match_codes:
                print(
                    f'{code} is not a valid code '
                    f'(from: {input_string})'
                )

        all_codes = [
            code for code in input_string_upp.split(' ')
            if code in self.all_match_codes
        ]
        if not expand:
            return all_codes

        # Extract the second element from tuples and skip None values
        broad_and_split = [
            item for tup in self.lookup_alleles(
                all_codes,
                codes_to_broad=self.codes_to_broad,
                codes_to_split=self.codes_to_split
            )
            for item in tup
            if item is not None
        ]

        # Extend all_codes in bulk
        all_codes.extend(broad_and_split)
        return all_codes


    def _count_broad_mismatches(
        self,
        d_hla: HLAProfile,
        p_hla: HLAProfile,
        loci: Optional[Set[str]] = None,
        locus_prefix: Optional[str] = 'mmb_'
    ) -> Dict[str, Optional[int]]:
        if loci is None:
            loci = set(mgr.ALL_HLA_LOCI)

        return {
            f'{locus_prefix}{sel_locus}': (
                len(d_hla.broads[i] - p_hla.broads[i])
            )
            for (i, sel_locus) in enumerate(mgr.ALL_HLA_LOCI)
            if sel_locus in loci
        }

    def _count_split_mismatches(
        self, d_hla: HLAProfile, p_hla: HLAProfile,
        loci: Optional[Set[str]] = None,
        locus_prefix: Optional[str] = 'mms_'
    ) -> Dict[str, Optional[int]]:
        if loci is None:
            loci = set(mgr.ALL_HLA_LOCI)

        mm = {
            f'{locus_prefix}{sel_locus}': len(
                d_hla.match_splits[i] - p_hla.match_splits[i]
            )
            for (i, sel_locus) in enumerate(self.needed_split_mismatches)
            if sel_locus in loci
        }

        for i_loc, locus in enumerate(mgr.ALL_HLA_LOCI):
            if (
                locus in loci and (
                    d_hla.has_unknown_splits[i_loc] or
                    p_hla.has_unknown_splits[i_loc]
                )
            ):
                common_broads = d_hla.broads[i_loc] & p_hla.broads[i_loc]
                r_match_hlas = p_hla.broads[i_loc] | set()
                d_match_hlas = d_hla.broads[i_loc] | set()
                for cb in common_broads:
                    if cb in self.splittable_broads[locus]:
                        rs = (
                            self.broads_to_splits[locus][cb] &
                            p_hla.splits[i_loc]
                        )
                        if rs:
                            ds = (
                                self.broads_to_splits[locus][cb] &
                                d_hla.splits[i_loc]
                            )
                            if ds:
                                r_match_hlas = (r_match_hlas | rs) - {cb}
                                d_match_hlas = (d_match_hlas | ds) - {cb}

                mm[f'{locus_prefix}{locus}'] = len(
                    d_match_hlas - r_match_hlas
                ) if r_match_hlas else None

        return mm

    def calculate_vpra_from_string(self, unacc_str: str) -> float:
        'Calculate a vPRA from an input string'
        if isinstance(unacc_str, str):
            unacceptables = Unacceptables(
                hla_system=self, unacc_string=unacc_str
            )
            return (
                sum(
                    len(unacceptables.unacceptables & a) > 0
                    for a in self.donor_pool_hlas
                ) / len(self.donor_pool_hlas)
            )
        else:
            return 0

    def count_mismatches(
        self,
        d_hla: HLAProfile,
        p_hla: HLAProfile,
        safely: bool = True
    ) -> Dict[str, int]:
        # Determine mismatches

        mm = self._count_broad_mismatches(
            d_hla=d_hla, p_hla=p_hla, loci=self.needed_broad_mismatches
        )
        mm |= self._count_split_mismatches(
            d_hla=d_hla, p_hla=p_hla, loci=self.needed_split_mismatches
        )

        if safely and len(mm) != self.k_needed_mms:
            return None
        mm_total = 0
        try:
            for v in mm.values():
                if v:
                    mm_total -= v
            mm[cn.MM_TOTAL] = mm_total
        except Exception:
            pass

        return mm

    @property
    def donor_pool_hlas(self):
        if self._donor_pool_hlas is not None:
            return self._donor_pool_hlas
        elif self.sim_set.PATH_DONOR_POOL is not None:
            df_don_pool = rdr.read_donor_pool(self.sim_set.PATH_DONOR_POOL)
            df_don_pool_hlas = df_don_pool.loc[
                :, [cn.ID_DONOR, cn.DONOR_HLA]
            ].drop_duplicates()
            df_don_pool_hlas.loc[:, cn.DONOR_HLA] = (
                rdr.fix_hla_string(df_don_pool_hlas.loc[:, cn.DONOR_HLA])
            )
            self._donor_pool_hlas = list(
                frozenset(self.return_all_antigens(s, expand=True))
                for s in df_don_pool_hlas.donor_hla
            )
            return self._donor_pool_hlas
        else:
            raise Exception(
                "Specify a path to a donor pool to "
                "calculate vPRAs in the yml file."
            )


class Unacceptables:
    """Class which implements unacceptables.

    Attributes
    ----------
    unacceptables_string

    Methods
    -------
    unacceptables
    """
    __slots__ = (
        '_unacceptables', 'unacc_string', 'hla_system'
    )

    def __init__(
        self, hla_system: HLASystem, unacc_string: Optional[str] = None
    ):
        self.unacc_string = unacc_string
        self._unacceptables = None
        self.hla_system = hla_system

    def add_unacceptables(self, new_unacceptables: FrozenSet[str]) -> None:
        """Add unacceptables to a pre-existing set of unacceptables."""
        self._unacceptables = self.unacceptables | new_unacceptables
        if self.unacc_string:
            self.unacc_string = (
                self.unacc_string + ' ' + ' '.join(new_unacceptables)
            )
        else:
            self.unacc_string = ' '.join(new_unacceptables)

    def add_missing_broad(self, broad, hla_string) -> str:

        # Split the original string into parts
        prefix = ''.join([char for char in broad if not char.isdigit()])
        parts = hla_string.split()


        # Identify the position of the first part that matches the prefix
        insert_index = None
        for i, part in enumerate(parts):
            if part.startswith(prefix):
                insert_index = i
                break

        # If 'broad' is in the list, remove it from its current position
        if broad in parts:
            parts.remove(broad)

        # Insert 'broad' before the first part that matches the prefix
        if insert_index is not None:
            parts.insert(insert_index, broad)
        else:
            # If no part with the prefix exists, append 'broad' at the end
            parts.append(broad)

        # Join the parts back into a string
        result = ' '.join(parts)

        return result

    def add_missing_broads(self) -> None:
        """Function which adds missing broads to a set of unacceptable antigens,
           if complete splits are present"""
        present_broads = set(
            self.hla_system.codes_to_broad[antigen]
            for antigen in self._unacceptables
            if antigen in self.hla_system.codes_to_broad
        )
        for locus, broad in present_broads:
            if (len(corresponding_splits := self.hla_system.broads_to_splits[locus][broad]) > 1):
                if len(corresponding_splits & self._unacceptables) == len(corresponding_splits):
                    self._unacceptables |= {broad}
                    self.unacc_string = self.add_missing_broad(broad, self.unacc_string)
        self._unacceptables = frozenset(self._unacceptables)


    @property
    def unacceptables(self) -> FrozenSet[str]:
        """property which returns a set of unacceptables,
           constructed from the unacceptable string"""
        if type(self._unacceptables) is frozenset:
            return self._unacceptables
        if type(self.unacc_string) is str:
            self._unacceptables = set(
                self.hla_system.return_all_antigens(
                    self.unacc_string
                )
            )
            self.add_missing_broads()

        else:
            self._unacceptables = frozenset()
        return self._unacceptables

    def __str__(self):
        if type(self.unacc_string) is str:
            return self.unacc_string
        else:
            return ''

    def __deepcopy__(self, memo):
        # Create a new instance of the class
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        # Copy all attributes, except for hla_system and bal_system
        for k, v in self.__dict__.items():
            if k in {'hla_system'}:
                setattr(result, k, v)  # Shallow copy
            else:
                setattr(result, k, deepcopy(v, memo))  # Deep copy

        return result
