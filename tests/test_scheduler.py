"""
Tests for the FluvTree scheduler: run an ordered ruleset over the canonical graph.

Two layers:
  1. Composition/ordering with lightweight synthetic processes (no GRLP): the
     ruleset runs in order, each process sees the prior ones' writes within a step,
     order changes the outcome, and canonical time advances on the graph.
  2. The scheduler driving the real GRLP process still reproduces standalone GRLP.
"""

import numpy as np
import pytest

from fluvtree import RiverNetwork, Scheduler


# A trivial 2-segment chain is enough for the composition tests.
_CHAIN = dict(up=[[], [0]], down=[[1], []])


def _chain():
    return RiverNetwork.from_segment_lists(_CHAIN["up"], _CHAIN["down"])


# --------------------------------------------------------------------------- #
# Synthetic processes (test doubles) operating on a graph-level scalar/log.
# --------------------------------------------------------------------------- #

class Recorder(object):
    """Appends its label to G.graph['log'] each step."""
    def __init__(self, network, label):
        self.network = network
        self.label = label

    def step(self, dt):
        self.network.graph.graph.setdefault("log", []).append(self.label)


class AddV(object):
    """v += delta each step."""
    def __init__(self, network, delta):
        self.network = network
        self.delta = delta

    def step(self, dt):
        g = self.network.graph.graph
        g["v"] = g.get("v", 0.0) + self.delta


class MulV(object):
    """v *= factor each step."""
    def __init__(self, network, factor):
        self.network = network
        self.factor = factor

    def step(self, dt):
        g = self.network.graph.graph
        g["v"] = g.get("v", 0.0) * self.factor


# --------------------------------------------------------------------------- #
# 1. Composition / ordering
# --------------------------------------------------------------------------- #

def test_runs_ruleset_in_order_each_step():
    rn = _chain()
    sched = Scheduler(rn, [Recorder(rn, "A"), Recorder(rn, "B")])
    sched.run(nt=3, dt=1.0)
    assert rn.graph.graph["log"] == ["A", "B", "A", "B", "A", "B"]


def test_within_step_visibility_and_order_matters():
    # v starts at 1; Add(1) then Mul(2) => (1+1)*2 = 4 in one step.
    rn = _chain()
    rn.graph.graph["v"] = 1.0
    Scheduler(rn, [AddV(rn, 1.0), MulV(rn, 2.0)]).step(1.0)
    assert rn.graph.graph["v"] == 4.0

    # Reverse order => 1*2 + 1 = 3. Order changes the outcome.
    rn2 = _chain()
    rn2.graph.graph["v"] = 1.0
    Scheduler(rn2, [MulV(rn2, 2.0), AddV(rn2, 1.0)]).step(1.0)
    assert rn2.graph.graph["v"] == 3.0


def test_advances_canonical_time_on_graph():
    rn = _chain()
    sched = Scheduler(rn)
    assert rn.graph.graph["t"] == 0.0
    t = sched.run(nt=4, dt=2.5)
    assert t == 10.0
    assert sched.t == 10.0
    assert rn.graph.graph["t"] == 10.0


def test_add_returns_process_and_appends():
    rn = _chain()
    sched = Scheduler(rn)
    p = Recorder(rn, "A")
    assert sched.add(p) is p
    assert sched.processes == [p]


def test_rejects_process_bound_to_foreign_network():
    rn = _chain()
    other = _chain()
    sched = Scheduler(rn)
    with pytest.raises(ValueError, match="different RiverNetwork"):
        sched.add(Recorder(other, "X"))


# --------------------------------------------------------------------------- #
# 2. Scheduler driving the real GRLP process
# --------------------------------------------------------------------------- #

def test_scheduler_drives_grlp_process_matches_standalone():
    grlp = pytest.importorskip("grlp")
    from fluvtree.processes import GRLP, build_network

    _D = 2000.0
    up, down = [[], [], [0, 1]], [[2], [2], []]
    x = [_D * np.arange(1, 5.0), _D * np.arange(1, 5.0), _D * np.arange(5, 9.0)]
    Q = [5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)]
    z0 = [np.linspace(60.0, 20.0, len(xi)) for xi in x]
    B = [100.0 * np.ones(len(xi)) for xi in x]
    S0, x_bl, z_bl = 0.015, _D * 9, 0.0

    def bdf2(net):
        net.set_time_integration(2)
        net.set_iteration_tolerance(1.0e-4)

    # ground truth
    gt = grlp.Network()
    gt.initialize(x_bl=x_bl, z_bl=z_bl, S0=[S0] * 2, Q_s_0=None,
                  upstream_segment_IDs=up, downstream_segment_IDs=down,
                  x=[xi.copy() for xi in x], z=[zi.copy() for zi in z0],
                  Q=[qi.copy() for qi in Q], B=[bi.copy() for bi in B])
    bdf2(gt)
    gt.get_z_lengths()

    # via the scheduler
    rn = build_network(up, down, x, z0, Q, B, S0, x_bl, z_bl)
    sched = Scheduler(rn, [GRLP(rn, configure=bdf2)])

    dt = 3.15e10
    for _ in range(6):
        gt.evolve_threshold_width_river_network(nt=1, dt=dt)
    sched.run(nt=6, dt=dt)

    for i, s in enumerate(rn.segment_ids):
        assert np.allclose(gt.segments[i].z, rn.get_segment_field(s, "z"),
                           rtol=0, atol=1e-9)
    assert rn.graph.graph["t"] == pytest.approx(6 * dt)
