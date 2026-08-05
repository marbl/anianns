"""Candidate selection helpers for multi-window satellite annotation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence


MAX_WINDOW_SIZES = 4


@dataclass(frozen=True)
class WindowCandidate:
    """One coarse satellite candidate and the window that detected it."""

    start: int
    end: int
    count: int
    window: int
    source: str = "diagonal"

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)

    @property
    def support_bp(self) -> int:
        return max(0, self.count) * self.window

    @property
    def density(self) -> float:
        return min(1.0, self.support_bp / self.length) if self.length else 0.0

    @property
    def selection_score(self) -> float:
        # A two-window resolution penalty makes the finer call win when spans
        # are similar, but becomes negligible when a coarse call recovers a
        # materially larger supported array.
        return min(self.length, self.support_bp) - (2 * self.window)


def normalize_window_sizes(values: int | Sequence[int]) -> tuple[int, ...]:
    """Validate, deduplicate, and sort CLI window sizes."""
    if isinstance(values, int):
        values = (values,)
    windows = tuple(sorted(set(int(value) for value in values)))
    if not windows or any(window <= 0 for window in windows):
        raise ValueError("window sizes must be positive integers")
    if len(windows) > MAX_WINDOW_SIZES:
        raise ValueError(
            f"at most {MAX_WINDOW_SIZES} window sizes may be requested in one run"
        )
    return windows


def derive_window_sizes(base_window: int) -> tuple[int, int, int]:
    """Return the half/base/double resolution pyramid for one CLI value."""
    base_window = int(base_window)
    if base_window < 2:
        raise ValueError("window size must be an integer of at least 2 bp")
    return (base_window // 2, base_window, base_window * 2)


def candidate_passes_support(candidate: WindowCandidate) -> bool:
    """Apply the historical count/length support filter at its source window."""
    if candidate.length <= 0:
        return False
    if candidate.support_bp <= candidate.length * 0.6:
        return False
    if candidate.count < 10 and candidate.support_bp <= candidate.length * 0.75:
        return False
    return True


def _overlap(first: WindowCandidate, second: WindowCandidate) -> int:
    return max(0, min(first.end, second.end) - max(first.start, second.start))


def candidates_compete(first: WindowCandidate, second: WindowCandidate) -> bool:
    """Return whether two calls are alternate resolutions of the same locus."""
    # Preserve the historical behavior within one resolution. Arbitration is
    # only intended to reconcile alternate calls introduced by multi-window
    # mode; same-window band pieces are handled by boundary refinement.
    if first.window == second.window and "periodic_diagonal" not in (
        first.source,
        second.source,
    ):
        return False
    overlap = _overlap(first, second)
    if overlap <= 0 or first.length <= 0 or second.length <= 0:
        return False
    shorter_fraction = overlap / min(first.length, second.length)
    # Full or near-full containment is an alternate call even if a finer
    # resolution split out only a small internal fragment of a long array.
    return shorter_fraction >= 0.8


def select_multi_window_candidates(
    candidates: Iterable[WindowCandidate],
) -> list[WindowCandidate]:
    """Filter candidates and retain the best nonredundant window per locus."""
    unique = {
        (
            candidate.start,
            candidate.end,
            candidate.count,
            candidate.window,
            candidate.source,
        ): candidate
        for candidate in candidates
        if candidate_passes_support(candidate)
    }
    ranked = sorted(
        unique.values(),
        key=lambda candidate: (
            -candidate.selection_score,
            -candidate.density,
            candidate.window,
            candidate.start,
            candidate.end,
        ),
    )
    selected: list[WindowCandidate] = []
    for candidate in ranked:
        if any(candidates_compete(candidate, retained) for retained in selected):
            continue
        selected.append(candidate)
    return sorted(selected, key=lambda candidate: (candidate.start, candidate.end))
