import pytest

from anianns.multi_window import (
    WindowCandidate,
    candidate_passes_support,
    derive_window_sizes,
    normalize_window_sizes,
    select_multi_window_candidates,
)


def test_normalize_window_sizes_sorts_deduplicates_and_limits_count():
    assert normalize_window_sizes([4000, 1000, 2000, 1000]) == (1000, 2000, 4000)
    with pytest.raises(ValueError, match="positive"):
        normalize_window_sizes([0, 1000])
    with pytest.raises(ValueError, match="at most 4"):
        normalize_window_sizes([500, 1000, 2000, 4000, 8000])


def test_derive_window_sizes_builds_half_base_double_pyramid():
    assert derive_window_sizes(5000) == (2500, 5000, 10_000)
    assert derive_window_sizes(5001) == (2500, 5001, 10_002)
    with pytest.raises(ValueError, match="at least 2"):
        derive_window_sizes(1)


def test_support_filter_uses_the_candidates_source_window():
    assert not candidate_passes_support(WindowCandidate(0, 10_000, 2, 1000))
    assert candidate_passes_support(WindowCandidate(0, 10_000, 3, 4000))


def test_fine_window_wins_when_supported_spans_are_comparable():
    selected = select_multi_window_candidates(
        [
            WindowCandidate(0, 10_000, 10, 1000),
            WindowCandidate(0, 10_000, 3, 4000),
        ]
    )

    assert selected == [WindowCandidate(0, 10_000, 10, 1000)]


def test_coarse_window_wins_when_it_recovers_materially_more_array():
    selected = select_multi_window_candidates(
        [
            WindowCandidate(10_000, 20_000, 10, 1000),
            WindowCandidate(8_000, 38_000, 8, 4000),
        ]
    )

    assert selected == [WindowCandidate(8_000, 38_000, 8, 4000)]


def test_distinct_calls_from_different_resolutions_are_both_retained():
    candidates = [
        WindowCandidate(0, 10_000, 10, 1000),
        WindowCandidate(50_000, 90_000, 10, 4000),
    ]

    assert select_multi_window_candidates(candidates) == candidates


def test_nested_fine_fragment_does_not_duplicate_a_long_coarse_call():
    selected = select_multi_window_candidates(
        [
            WindowCandidate(10_000, 20_000, 10, 1000),
            WindowCandidate(0, 100_000, 25, 4000),
        ]
    )

    assert selected == [WindowCandidate(0, 100_000, 25, 4000)]


def test_same_window_candidates_retain_legacy_non_reconciliation_behavior():
    candidates = [
        WindowCandidate(0, 20_000, 20, 1000),
        WindowCandidate(10_000, 30_000, 20, 1000),
    ]

    assert select_multi_window_candidates(candidates) == candidates


def test_periodic_rescue_reconciles_with_same_window_fragment():
    periodic = WindowCandidate(
        0,
        100_000,
        80,
        1000,
        source="periodic_diagonal",
    )
    selected = select_multi_window_candidates(
        [periodic, WindowCandidate(80_000, 100_000, 20, 1000)]
    )

    assert selected == [periodic]
