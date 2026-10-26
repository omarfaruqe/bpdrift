"""Pin a breakpoint table file by version, retrieval date and SHA-256.

bpdrift never ships breakpoint data. The user downloads the official file;
a pin records exactly which file every interpretation came from, so a
result can be reproduced or audited later.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

_CHUNK = 1 << 20


class PinError(ValueError):
    """Raised when a file does not match its pin, or a pin is malformed."""


@dataclass(frozen=True)
class Pin:
    guideline: str       # "EUCAST"
    version: str         # "16.1"
    retrieved: str       # ISO date the file was downloaded, "2026-10-26"
    sha256: str
    size_bytes: int
    filename: str

    def __post_init__(self) -> None:
        if not self.guideline.strip() or not self.version.strip():
            raise PinError("guideline and version must not be empty")
        try:
            dt.date.fromisoformat(self.retrieved)
        except ValueError as exc:
            raise PinError(f"retrieved must be an ISO date, got {self.retrieved!r}") from exc
        if len(self.sha256) != 64 or any(c not in "0123456789abcdef" for c in self.sha256):
            raise PinError("sha256 must be 64 lowercase hex characters")

    @property
    def label(self) -> str:
        return f"{self.guideline} v{self.version} ({self.retrieved}, sha256:{self.sha256[:12]})"


def sha256_of(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def pin_file(
    path: str | Path,
    version: str,
    retrieved: str | dt.date,
    guideline: str = "EUCAST",
) -> Pin:
    """Create a pin for a breakpoint file the user has downloaded."""
    p = Path(path)
    if not p.is_file():
        raise PinError(f"no such file: {p}")
    if isinstance(retrieved, dt.date):
        retrieved = retrieved.isoformat()
    return Pin(
        guideline=guideline,
        version=str(version),
        retrieved=retrieved,
        sha256=sha256_of(p),
        size_bytes=p.stat().st_size,
        filename=p.name,
    )


def verify(path: str | Path, pin: Pin) -> None:
    """Raise PinError unless the file at path is exactly the pinned file."""
    actual = sha256_of(path)
    if actual != pin.sha256:
        raise PinError(
            f"{Path(path).name} does not match pin {pin.label}: sha256 is {actual[:12]}…"
        )


def save_pin(pin: Pin, path: str | Path) -> None:
    Path(path).write_text(json.dumps(asdict(pin), indent=2) + "\n", encoding="utf-8")


def load_pin(path: str | Path) -> Pin:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return Pin(**data)
    except (json.JSONDecodeError, TypeError) as exc:
        raise PinError(f"malformed pin file {path}: {exc}") from exc
