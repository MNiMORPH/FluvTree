"""
Tests for the FluvTree network substrate (fluvtree.network.RiverNetwork).

Two layers:
  1. Pure-substrate topology/state assertions (no GRLP needed).
  2. Cross-validation of the traversal against GRLP ground truth -- the substrate
     is extracted from GRLP's physics-free machinery, so on a shared topology the
     upstream/downstream walks must agree with grlp.Network. Skipped if grlp is
     unavailable.
"""

import numpy as np
import pytest

from fluvtree import RiverNetwork


# --------------------------------------------------------------------------- #
# Shared small topology catalog (up/down lists, as GRLP specifies a network).
#
#   Y            : two heads -> one confluence -> outlet         (3 segments)
#   multi_level  : four heads through two confluence levels      (7 segments)
# --------------------------------------------------------------------------- #

TOPOLOGIES = {
    "Y": dict(
        up=[[], [], [0, 1]],
        down=[[2], [2], []],
    ),
    "multi_level": dict(
        up=[[], [], [], [], [0, 1], [4, 2], [5, 3]],
        down=[[4], [4], [5], [6], [5], [6], []],
    ),
}


def _rn(name):
    """RiverNetwork for a named topology."""
    t = TOPOLOGIES[name]
    return RiverNetwork.from_segment_lists(t["up"], t["down"])


# --------------------------------------------------------------------------- #
# 1. Pure substrate
# --------------------------------------------------------------------------- #

def test_Y_heads_and_mouth():
    rn = _rn("Y")
    assert sorted(rn.head_segments()) == [0, 1]
    assert rn.mouth_segments() == [2]


def test_Y_immediate_neighbours():
    rn = _rn("Y")
    assert sorted(rn.upstream_segments(2)) == [0, 1]
    assert rn.upstream_segments(0) == []
    assert rn.downstream_segment(0) == 2
    assert rn.downstream_segment(1) == 2
    assert rn.downstream_segment(2) is None


def test_Y_walks():
    rn = _rn("Y")
    assert rn.walk_downstream(0) == [0, 2]        # ordered to outlet
    assert rn.walk_downstream(2) == [2]
    assert sorted(rn.walk_upstream(2)) == [0, 1, 2]


def test_multi_level_walks():
    rn = _rn("multi_level")
    assert sorted(rn.head_segments()) == [0, 1, 2, 3]
    assert rn.mouth_segments() == [6]
    # segment 0 flows 0 -> 4 -> 5 -> 6 to the outlet
    assert rn.walk_downstream(0) == [0, 4, 5, 6]
    # everything is upstream of the mouth
    assert sorted(rn.walk_upstream(6)) == [0, 1, 2, 3, 4, 5, 6]


def test_reach_profile_assembly():
    """reach_profile hstacks the shared node endpoints around the edge interior."""
    rn = _rn("Y")
    rn.set_segment_field(0, "z", np.array([10.0, 9.0]))
    u, v = rn.edge_of(0)
    rn.set_node_field(u, "z", 10.5)   # head endpoint
    rn.set_node_field(v, "z", 8.5)    # confluence endpoint
    assert np.allclose(rn.reach_profile(0, "z"), [10.5, 10.0, 9.0, 8.5])


def test_non_convergent_topology_rejected():
    """A segment with two downstream reaches is not a convergent tree."""
    G = _rn("Y").graph
    # give segment 0's downstream node a second out-edge
    G.add_edge(("jcn", 2), ("outlet2",), seg=99)
    rn = RiverNetwork(G)
    with pytest.raises(ValueError, match="non-convergent"):
        rn.downstream_segment(0)


# --------------------------------------------------------------------------- #
# 2. Cross-validation against GRLP ground truth
# --------------------------------------------------------------------------- #

def _build_grlp_network(up, down):
    """Minimal grlp.Network with a given topology (no evolution needed)."""
    grlp = pytest.importorskip("grlp")
    n = len(up)
    n_heads = sum(1 for i in range(n) if len(up[i]) == 0)
    x = [2000.0 * np.arange(1, 5, dtype=float) for _ in range(n)]
    net = grlp.Network()
    net.initialize(
        x_bl=2000.0 * 100,
        z_bl=0.0,
        S0=[0.015] * n_heads,
        Q_s_0=None,
        upstream_segment_IDs=up,
        downstream_segment_IDs=down,
        x=x,
        z=[np.zeros(len(xi)) for xi in x],
        Q=[np.ones(len(xi)) for xi in x],
        B=[100.0 * np.ones(len(xi)) for xi in x],
    )
    net.build_graph()
    return net


@pytest.mark.parametrize("name", list(TOPOLOGIES))
def test_traversal_matches_grlp(name):
    up, down = TOPOLOGIES[name]["up"], TOPOLOGIES[name]["down"]
    net = _build_grlp_network(up, down)
    rn = RiverNetwork.from_segment_lists(up, down)
    for s in range(len(up)):
        assert set(rn.walk_upstream(s)) == set(net.find_upstream_IDs(s))
        assert rn.walk_downstream(s) == list(net.find_downstream_IDs(s))
