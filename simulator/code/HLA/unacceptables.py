from __future__ import annotations

from copy import deepcopy
from typing import FrozenSet, Optional

from .ontology import HLAOntology


class Unacceptables:
    """Implements unacceptable antigens. It expands missing broads by default."""

    __slots__ = ("_unacceptables", "unacc_string", "ontology")

    def __init__(
        self,
        ontology: HLAOntology,
        unacc_string: Optional[str] = None,
    ):
        self.unacc_string = unacc_string
        self._unacceptables = None
        self.ontology = ontology

    def add_unacceptables(self, new_unacceptables: FrozenSet[str]) -> None:
        self._unacceptables = self.unacceptables | new_unacceptables
        if self.unacc_string:
            self.unacc_string = self.unacc_string + " " + " ".join(new_unacceptables)
        else:
            self.unacc_string = " ".join(new_unacceptables)

    @staticmethod
    def add_missing_broad(broad: str, hla_string: str) -> str:
        prefix = "".join([char for char in broad if not char.isdigit()])
        parts = hla_string.split()

        insert_index = None
        for i, part in enumerate(parts):
            if part.startswith(prefix):
                insert_index = i
                break

        if broad in parts:
            parts.remove(broad)

        if insert_index is not None:
            parts.insert(insert_index, broad)
        else:
            parts.append(broad)

        return " ".join(parts)

    def add_missing_broads(self) -> None:
        present_broads = {
            self.ontology.codes_to_broad[antigen]
            for antigen in self._unacceptables
            if antigen in self.ontology.codes_to_broad
        }

        for locus, broad in present_broads:
            corresponding_splits = self.ontology.broads_to_splits[locus][broad]
            if len(corresponding_splits) > 1:
                if len(corresponding_splits & self._unacceptables) == len(corresponding_splits):
                    self._unacceptables |= {broad}
                    self.unacc_string = self.add_missing_broad(broad, self.unacc_string)

        self._unacceptables = frozenset(self._unacceptables)

    @property
    def unacceptables(self) -> FrozenSet[str]:
        if type(self._unacceptables) is frozenset:
            return self._unacceptables
        if type(self.unacc_string) is str:
            normalized_unacc_string = self._normalized_unacc_string() or ""
            self._unacceptables = set(
                self.ontology.return_all_antigens(normalized_unacc_string)
            )
            self.add_missing_broads()
        else:
            self._unacceptables = frozenset()
        return self._unacceptables

    def __str__(self):
        if type(self.unacc_string) is str:
            return self.unacc_string
        return ""

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

    def _normalized_unacc_string(self) -> Optional[str]:
        if type(self.unacc_string) is not str:
            return None
        return self.ontology._normalize_hla_input_string(self.unacc_string)
