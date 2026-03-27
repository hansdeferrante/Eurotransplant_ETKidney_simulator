from __future__ import annotations

from typing import Dict, Iterable, Tuple, Union

import simulator.magic_values.etkidney_simulator_settings as es


MismatchKey = Union[str, Tuple[int, int, int]]
MismatchSource = Union[str, Iterable[int], Dict[MismatchKey, int]]


class HLAMismatchCriteria:
    __slots__ = ("index_bits",)

    def __init__(self, source: MismatchSource):
        self.index_bits = self._parse_to_index_bits(source)

    @classmethod
    def from_index_tuple(cls, index_bits: Tuple[int, ...]) -> HLAMismatchCriteria:
        obj = cls.__new__(cls)
        obj.index_bits = cls._validate_index_bits(index_bits)
        return obj

    def accepts_index(self, mq_index: int) -> bool:
        if not isinstance(mq_index, int):
            raise TypeError("mq_index must be int.")
        if mq_index < 0 or mq_index >= len(self.index_bits):
            raise ValueError(
                f"mq_index out of range [0, {len(self.index_bits) - 1}]: {mq_index}"
            )
        return bool(self.index_bits[mq_index])

    def accepts_tuple(self, mq: Tuple[int, int, int]) -> bool:
        try:
            mm_a = int(mq[0])
            mm_b = int(mq[1])
            mm_dr = int(mq[2])
            self._validate_mismatch_components(mm_a, mm_b, mm_dr)
            return bool(self.index_bits[self._index_from_tuple(mm_a, mm_b, mm_dr)])
        except Exception as exc:
            raise ValueError(f"Invalid mismatch tuple: {mq}") from exc

    def to_index_tuple(self) -> Tuple[int, ...]:
        return self.index_bits

    def to_tuple_map(self) -> Dict[Tuple[int, int, int], int]:
        return {
            mq: self.index_bits[i]
            for i, mq in enumerate(es.HLA_MQS)
        }

    @staticmethod
    def _validate_binary_value(value: int) -> int:
        if value not in (0, 1, True, False):
            raise ValueError("Match quality values must be 0/1.")
        return int(value)

    @staticmethod
    def _validate_mismatch_components(mm_a: int, mm_b: int, mm_dr: int) -> None:
        if not (0 <= mm_a <= 2 and 0 <= mm_b <= 2 and 0 <= mm_dr <= 2):
            raise ValueError("Match quality mismatch components must be in [0, 2].")

    @staticmethod
    def _index_from_tuple(mm_a: int, mm_b: int, mm_dr: int) -> int:
        return (mm_dr * 9) + (mm_b * 3) + mm_a

    @classmethod
    def _validate_index_bits(cls, bits: Tuple[int, ...]) -> Tuple[int, ...]:
        if len(bits) != len(es.HLA_MQS_STR):
            raise ValueError(
                f"Match quality must be {len(es.HLA_MQS_STR)} digits, got {len(bits)}."
            )
        return tuple(cls._validate_binary_value(v) for v in bits)

    @classmethod
    def _parse_to_index_bits(cls, source: MismatchSource) -> Tuple[int, ...]:
        if isinstance(source, str):
            compact = "".join(source.replace(",", " ").split())
            if not compact:
                raise ValueError("Match quality string is empty.")
            if any(ch not in {"0", "1"} for ch in compact):
                raise ValueError("Match quality string must contain only 0/1 digits.")
            if len(compact) != len(es.HLA_MQS_STR):
                raise ValueError(
                    f"Match quality must be {len(es.HLA_MQS_STR)} digits, got {len(compact)}."
                )
            return tuple(int(ch) for ch in compact)

        if isinstance(source, dict):
            parsed = [0] * len(es.HLA_MQS_STR)
            for key, value in source.items():
                parsed_value = cls._validate_binary_value(value)

                if isinstance(key, str):
                    mq_str = "".join(key.replace(",", " ").split())
                    if len(mq_str) != 3 or not mq_str.isdigit():
                        raise ValueError(
                            "Match quality dict keys must be 3-digit strings or (A,B,DR) tuples."
                        )
                    mm_a = int(mq_str[0])
                    mm_b = int(mq_str[1])
                    mm_dr = int(mq_str[2])
                else:
                    mm_a = int(key[0])
                    mm_b = int(key[1])
                    mm_dr = int(key[2])

                cls._validate_mismatch_components(mm_a, mm_b, mm_dr)
                parsed[cls._index_from_tuple(mm_a, mm_b, mm_dr)] = parsed_value
            return tuple(parsed)

        digits_list = []
        for value in source:
            digits_list.append(str(cls._validate_binary_value(value)))
        digits = "".join(digits_list)
        if len(digits) != len(es.HLA_MQS_STR):
            raise ValueError(
                f"Match quality must be {len(es.HLA_MQS_STR)} digits, got {len(digits)}."
            )
        return tuple(int(ch) for ch in digits)
