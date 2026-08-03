"""
Synthesis: random drainage-network topology (Shreve 1974) in ``fluvtree.generate``.

Pins the Shreve invariants (a magnitude-``m`` network has ``m`` sources and
``2m - 1`` segments), that it unfolds into a valid convergent RiverNetwork, that a
``seed`` makes it reproducible, and the optional length/area assignment.
"""

import numpy as np
import pytest

from fluvtree.generate import ShreveRandomNetwork, random_network
from fluvtree.generate.topology import _downstream_IDs


@pytest.mark.parametrize("m", [1, 2, 5, 12, 30])
def test_magnitude_gives_shreve_counts(m):
    g = ShreveRandomNetwork(magnitude=m, seed=0)
    n = len(g.upstream_segment_IDs)
    heads = sum(1 for u in g.upstream_segment_IDs if not u)
    assert n == 2 * m - 1          # m sources + (m-1) confluences
    assert heads == m
    if m >= 2:                     # m=1 is the degenerate single-segment (no links)
        assert g.N_external == m


def test_builds_valid_convergent_network():
    rn = random_network(8, seed=1)
    assert len(rn.segment_ids) == 2 * 8 - 1
    assert len(rn.head_segments()) == 8
    assert len(rn.mouth_segments()) == 1
    # convergent: every non-mouth segment has exactly one downstream neighbour
    mouths = set(rn.mouth_segments())
    for s in rn.segment_ids:
        if s not in mouths:
            assert rn.downstream_segment(s) is not None


def test_magnitude_one_is_a_single_segment():
    rn = random_network(1, seed=0)
    assert len(rn.segment_ids) == 1
    assert rn.head_segments() == rn.mouth_segments()   # the one reach is both


def test_reproducible_with_seed():
    a = ShreveRandomNetwork(magnitude=15, seed=42).downstream_segment_IDs
    b = ShreveRandomNetwork(magnitude=15, seed=42).downstream_segment_IDs
    c = ShreveRandomNetwork(magnitude=15, seed=43).downstream_segment_IDs
    assert a == b
    assert a != c                                       # different seed, different draw


def test_prebuilt_topology_is_unfolded():
    # a supplied link string reproduces the same topology it was drawn from
    drawn = ShreveRandomNetwork(magnitude=7, seed=3)
    rebuilt = ShreveRandomNetwork(magnitude=7, topology=list(drawn.links))
    assert rebuilt.upstream_segment_IDs == drawn.upstream_segment_IDs
    assert rebuilt.downstream_segment_IDs == drawn.downstream_segment_IDs


def test_segment_lengths_rescaled_to_max_length():
    g = ShreveRandomNetwork(magnitude=6, seed=0,
                            segment_length=1000.0, max_length=5000.0)
    paths = [sum(g.segment_lengths[j]
                 for j in _downstream_IDs(g.downstream_segment_IDs, i))
             for i in range(len(g.segment_lengths))]
    assert np.isclose(max(paths), 5000.0)


def test_segment_areas_from_ratio_and_supply():
    g = ShreveRandomNetwork(magnitude=5, seed=0, segment_length=1000.0,
                            segment_length_area_ratio=2.0, supply_area=1.0e6)
    assert g.segment_areas is not None
    assert len(g.segment_areas) == len(g.upstream_segment_IDs)
    # source areas only on channel heads
    for i, up in enumerate(g.upstream_segment_IDs):
        if up:
            assert g.source_areas[i] == 0.0
        else:
            assert g.source_areas[i] == 1.0e6
