"""Parse MIC measurements, keeping their censoring sign.

An MIC is a value in mg/L with a sign:

    "="   exact reading          0.5
    "<="  at or below the value  <=0.25   (lowest dilution tested)
    "<"   below the value        <0.25
    ">="  at or above the value  >=32
    ">"   above the value        >16      (growth at highest dilution)

The sign matters: a censored value cannot always be placed on one side of
a breakpoint, and interpretation must say so rather than guess.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

SIGNS = ("=", "<=", "<", ">=", ">")

# Units equivalent to mg/L. Anything else is refused.
_MG_PER_L = {"mg/l", "ug/ml", "µg/ml", "μg/ml", "mcg/ml"}

_SIGN_ALIASES = {
    "≤": "<=",
    "=<": "<=",
    "≥": ">=",
    "=>": ">=",
    "==": "=",
}

_PATTERN = re.compile(
    r"""^\s*
    (?P<sign><=|>=|=<|=>|==|≤|≥|<|>|=)?   # optional sign
    \s*
    (?P<value>\d+(?:\.\d+)?|\.\d+)        # number
    \s*
    (?P<unit>[A-Za-zµμ/]+)?               # optional unit
    \s*$""",
    re.VERBOSE,
)


class MICParseError(ValueError):
    """Raised when an MIC string cannot be parsed unambiguously."""


@dataclass(frozen=True)
class MIC:
    """One MIC measurement in mg/L."""

    value: float
    sign: str = "="

    def __post_init__(self) -> None:
        if self.sign not in SIGNS:
            raise MICParseError(f"unknown sign {self.sign!r}; expected one of {SIGNS}")
        if not math.isfinite(self.value) or self.value <= 0:
            raise MICParseError(f"MIC value must be a positive finite number, got {self.value!r}")

    @property
    def is_censored(self) -> bool:
        return self.sign != "="

    def __str__(self) -> str:
        sign = "" if self.sign == "=" else self.sign
        return f"{sign}{self.value:g}"


def parse_mic(text: str | float | int) -> MIC:
    """Parse an MIC such as ``"<=0.25"``, ``"≥ 32 mg/L"`` or ``0.5``.

    Units, if given, must be equivalent to mg/L (mg/L, ug/mL, µg/mL).
    Raises MICParseError on anything else, including empty input.
    """
    if isinstance(text, bool):
        raise MICParseError("a boolean is not an MIC")
    if isinstance(text, (int, float)):
        return MIC(float(text))
    if not isinstance(text, str):
        raise MICParseError(f"cannot parse MIC from {type(text).__name__}")

    match = _PATTERN.match(text)
    if match is None:
        raise MICParseError(f"cannot parse MIC from {text!r}")

    sign = match.group("sign") or "="
    sign = _SIGN_ALIASES.get(sign, sign)

    unit = match.group("unit")
    if unit is not None and unit.lower() not in _MG_PER_L:
        raise MICParseError(f"unsupported unit {unit!r} in {text!r}; expected mg/L")

    return MIC(float(match.group("value")), sign)
