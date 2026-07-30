"""
n = 1 stream-power incision (fluvtree.solvers.advection) via StreamPower.

Checks the physics against closed form: steady-state slope-discharge
``S = U/(K Q^m)``, an upstream-migrating knickpoint, the Scheduler one-liner, and
unconditional implicit stability.
"""

import numpy as np
import pytest

from fluvtree import Scheduler
from fluvtree.processes import StreamPower, build_network


_D = 1000.0
K, M = 1e-9, 0.5
U = 1e-3 / 3.15e7          # 1 mm/yr in m/s


def _reach(nx=20, Q=10.0):
    x = _D * np.arange(1, nx + 1)
    return build_network(
        [[]], [[]], [x], [np.zeros(nx)], [Q * np.ones(nx)], [np.ones(nx)],
        S0=0.01, x_bl=_D * (nx + 1), z_bl=0.0)


def test_single_reach_steady_state_slope_discharge():
    rn = _reach(Q=10.0)
    StreamPower(rn, K=K, m=M, U=U).step(dt=1e12, nt=300)
    x = rn.get_segment_field(0, "x")
    z = rn.get_segment_field(0, "z")
    S_analytic = U / (K * 10.0 ** M)
    z_analytic = 0.0 + S_analytic * (_D * 21 - x)
    assert np.allclose(z, z_analytic, rtol=0, atol=1e-6)


def test_network_steady_state_slopes_follow_discharge():
    # Y: two heads (Q=5) -> confluence -> outlet (Q=10). Steady slope S=U/(K Q^m),
    # so the higher-discharge trunk is gentler.
    up, down = [[], [], [0, 1]], [[2], [2], []]
    x = [_D * np.arange(1, 5.0), _D * np.arange(1, 5.0), _D * np.arange(5, 9.0)]
    Q = [5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)]
    rn = build_network(up, down, x, [np.zeros(4)] * 3, Q,
                       [np.ones(4)] * 3, S0=0.01, x_bl=_D * 9, z_bl=0.0)
    StreamPower(rn, K=K, m=M, U=U).step(dt=1e12, nt=800)
    # higher-discharge trunk (Q=10) is gentler than the heads (Q=5)
    S_head = -np.diff(rn.get_segment_field(0, "z")) / np.diff(rn.get_segment_field(0, "x"))
    S_trunk = -np.diff(rn.get_segment_field(2, "z")) / np.diff(rn.get_segment_field(2, "x"))
    assert S_trunk.mean() < S_head.mean()
    for s, q in ((0, 5.0), (1, 5.0), (2, 10.0)):
        z = rn.get_segment_field(s, "z")
        xs = rn.get_segment_field(s, "x")
        S = -np.diff(z) / np.diff(xs)
        assert np.allclose(S, U / (K * q ** M), rtol=1e-4, atol=0)


def test_knickpoint_migrates_upstream():
    rn = _reach(Q=10.0)
    proc = StreamPower(rn, K=K, m=M, U=U)
    proc.step(dt=1e12, nt=300)                       # to steady state
    z_ss = rn.get_segment_field(0, "z").copy()
    outlet = rn.edge_of(rn.mouth_segments()[0])[1]
    rn.set_node_field(outlet, "z_bl", -60.0)         # base-level drop
    proc.step(dt=1e11, nt=15)                        # partial: t << L/celerity
    drop = z_ss - rn.get_segment_field(0, "z")       # how much each node has fallen
    # signal enters from the outlet: the downstream end has dropped, the head has not
    assert drop[-1] > 10.0                            # downstream node responded
    assert drop[0] < 1.0                              # head not yet reached
    assert drop[-1] > 10.0 * drop[0]                  # monotone upstream propagation


def test_scheduler_one_liner_runs_bedrock():
    rn = _reach(Q=10.0)
    Scheduler(rn, [StreamPower(rn, K=K, m=M, U=U)]).run(nt=200, dt=1e12)
    z = rn.get_segment_field(0, "z")
    assert np.all(np.isfinite(z))
    assert np.all(np.diff(z) < 0)                     # descends downstream
    assert rn.graph.graph["t"] == pytest.approx(200 * 1e12)


def test_implicit_stability_huge_dt():
    rn = _reach(Q=10.0)
    StreamPower(rn, K=K, m=M, U=U).step(dt=1e18, nt=3)   # absurd dt
    assert np.all(np.isfinite(rn.get_segment_field(0, "z")))
