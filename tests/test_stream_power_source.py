"""
Distributed vertical forcing on the bedrock river.

``advection`` now takes the same ``source_rate`` argument ``diffusion`` does, and
``StreamPower`` builds it from ``uplift``/``source_sink`` through
``common.rates.broadcast_rate`` -- so a forcing specified once drives either
door. These tests pin that the forcing arrives, that it composes the way the
module docstring says, and that the scalar interface is untouched.
"""

import numpy as np
import pytest

from fluvtree import Scheduler
from fluvtree.processes import StreamPower, build_network
from fluvtree.solvers.advection import evolve_streampower_n1, sweep_order, incise_n1_step


_D = 1000.0
K, M = 1e-9, 0.5
U = 1e-3 / 3.15e7          # 1 mm/yr in m/s


def _reach(nx=20, Q=10.0):
    x = _D * np.arange(1, nx + 1)
    return build_network(
        [[]], [[]], [x], [np.zeros(nx)], [Q * np.ones(nx)], [np.ones(nx)],
        S0=0.01, x_bl=_D * (nx + 1), z_bl=0.0)


def _steady_slope(U_rate, Q=10.0):
    return U_rate / (K * Q**M)


# ------------------------------------------------- the scalar path is untouched

def test_scalar_U_is_unchanged_by_the_new_argument():
    a, b = _reach(), _reach()
    evolve_streampower_n1(a, K, M, U, dt=1e12, nt=200)
    # the same uplift given through source_rate instead
    evolve_streampower_n1(b, K, M, 0.0, dt=1e12, nt=200, source_rate=U)
    assert np.allclose(a.get_segment_field(0, "z"), b.get_segment_field(0, "z"),
                       rtol=0, atol=1e-9)


def test_the_three_contributions_sum():
    # U + source_rate must behave as one rate: half in each gives the same
    # profile as all of it in either.
    both, only_U = _reach(), _reach()
    evolve_streampower_n1(both, K, M, 0.5 * U, dt=1e12, nt=200,
                          source_rate=0.5 * U)
    evolve_streampower_n1(only_U, K, M, U, dt=1e12, nt=200)
    assert np.allclose(both.get_segment_field(0, "z"),
                       only_U.get_segment_field(0, "z"), rtol=0, atol=1e-9)


# ------------------------------------------- spatially variable uplift, at last

def test_uplift_stepped_across_a_fault_steps_the_steady_slope():
    # The capability the process did not have. At steady state
    # S = U/(K Q**m), so a step in U is a step in S in the same ratio.
    nx = 30
    rn = _reach(nx=nx)
    U_arr = np.full(nx, U)
    U_arr[: nx // 2] = 3.0 * U            # upstream block rises three times faster
    proc = StreamPower(rn, K=K, m=M, U=0.0, uplift=[U_arr])
    Scheduler(rn, [proc]).run(nt=400, dt=1e12)

    x = rn.get_segment_field(0, "x")
    z = rn.get_segment_field(0, "z")
    S = -np.diff(np.append(z, 0.0)) / np.diff(np.append(x, _D * (nx + 1)))
    assert np.allclose(S, U_arr / (K * 10.0**M), rtol=1e-5, atol=0)


def test_uplift_as_a_callable_of_x():
    # broadcast_rate's callable form, which is why it is worth reusing rather
    # than requiring a precomputed array.
    nx = 30
    rn = _reach(nx=nx)
    proc = StreamPower(rn, K=K, m=M, U=0.0,
                       uplift=lambda x: np.where(x < 15000.0, 3.0 * U, U))
    Scheduler(rn, [proc]).run(nt=400, dt=1e12)

    x = rn.get_segment_field(0, "x")
    z = rn.get_segment_field(0, "z")
    S = -np.diff(np.append(z, 0.0)) / np.diff(np.append(x, _D * (nx + 1)))
    expected = np.where(x < 15000.0, 3.0 * U, U) / (K * 10.0**M)
    assert np.allclose(S, expected, rtol=1e-5, atol=0)


def test_source_sink_sums_with_uplift():
    nx = 20
    rn_split, rn_whole = _reach(nx=nx), _reach(nx=nx)
    Scheduler(rn_split, [StreamPower(rn_split, K=K, m=M, U=0.0,
                                     uplift=0.5 * U, source_sink=0.5 * U)]
              ).run(nt=200, dt=1e12)
    Scheduler(rn_whole, [StreamPower(rn_whole, K=K, m=M, U=U)]
              ).run(nt=200, dt=1e12)
    assert np.allclose(rn_split.get_segment_field(0, "z"),
                       rn_whole.get_segment_field(0, "z"), rtol=0, atol=1e-9)


def test_no_forcing_means_no_source_rate_is_built():
    rn = _reach()
    proc = StreamPower(rn, K=K, m=M, U=U)
    assert proc._static_source_rate() is None        # nothing to broadcast


# ----------------------------------------------------- the network case

def test_distributed_uplift_on_a_network():
    # Y: two heads into a trunk. Uplift given per segment, the trunk slower.
    up, down = [[], [], [0, 1]], [[2], [2], []]
    x = [_D * np.arange(1, 5.0), _D * np.arange(1, 5.0), _D * np.arange(5, 9.0)]
    Q = [5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)]
    rn = build_network(up, down, x, [np.zeros(4)] * 3, Q, [np.ones(4)] * 3,
                       S0=0.01, x_bl=_D * 9, z_bl=0.0)
    rates = [2.0 * U * np.ones(4), 2.0 * U * np.ones(4), U * np.ones(4)]
    Scheduler(rn, [StreamPower(rn, K=K, m=M, U=0.0, uplift=rates)]
              ).run(nt=600, dt=1e12)

    for s, q, u in ((0, 5.0, 2.0 * U), (1, 5.0, 2.0 * U), (2, 10.0, U)):
        z = rn.get_segment_field(s, "z")
        xs = rn.get_segment_field(s, "x")
        S = -np.diff(z) / np.diff(xs)
        assert np.allclose(S, u / (K * q**M), rtol=1e-4, atol=0), f"segment {s}"
