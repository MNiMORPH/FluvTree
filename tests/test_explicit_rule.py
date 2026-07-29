"""
Tests for the door-2 explicit-rule process.

Validated against analytic rates whose forward-Euler trajectory is an exact closed
form, so the mechanism is checked without reference to any physics:
  - uniform uplift    dz/dt = U        => z_n = z0 + U*n*dt
  - linear decay      dz/dt = -k z      => z_n = z0*(1 - k dt)^n   (exact for Euler)
"""

import numpy as np
import pytest

from fluvtree import RiverNetwork, Scheduler
from fluvtree.processes import ExplicitRule


_CHAIN = dict(up=[[], [0]], down=[[1], []])


def _chain_with_z():
    rn = RiverNetwork.from_segment_lists(_CHAIN["up"], _CHAIN["down"])
    for s in rn.segment_ids:
        rn.set_segment_field(s, "z", np.array([10.0, 20.0, 30.0]))
    return rn


def test_uniform_uplift_exact():
    rn = _chain_with_z()
    U = 1e-3
    z0 = {s: rn.get_segment_field(s, "z").copy() for s in rn.segment_ids}

    def uplift(net, dt):
        return {s: U * np.ones_like(net.get_segment_field(s, "z"))
                for s in net.segment_ids}

    rule = ExplicitRule(rn, uplift)
    nt, dt = 7, 3.15e7
    for _ in range(nt):
        rule.step(dt)
    for s in rn.segment_ids:
        assert np.allclose(rn.get_segment_field(s, "z"), z0[s] + U * nt * dt)


def test_linear_decay_exact():
    rn = _chain_with_z()
    k = 0.1
    z0 = {s: rn.get_segment_field(s, "z").copy() for s in rn.segment_ids}

    def decay(net, dt):
        return {s: -k * net.get_segment_field(s, "z") for s in net.segment_ids}

    rule = ExplicitRule(rn, decay)
    dt, nt = 1.0, 5
    for _ in range(nt):
        rule.step(dt)
    for s in rn.segment_ids:
        assert np.allclose(rn.get_segment_field(s, "z"), z0[s] * (1 - k * dt) ** nt)


def test_omitted_segments_unchanged():
    rn = _chain_with_z()
    z_seg1 = rn.get_segment_field(1, "z").copy()

    def only_seg0(net, dt):
        return {0: np.ones_like(net.get_segment_field(0, "z"))}

    ExplicitRule(rn, only_seg0).step(1.0)
    assert np.allclose(rn.get_segment_field(1, "z"), z_seg1)  # untouched


def test_evolves_named_field():
    rn = _chain_with_z()
    for s in rn.segment_ids:
        rn.set_segment_field(s, "sed", np.zeros(3))

    def deposit(net, dt):
        return {s: np.ones(3) for s in net.segment_ids}

    ExplicitRule(rn, deposit, field="sed").step(2.0)
    for s in rn.segment_ids:
        assert np.allclose(rn.get_segment_field(s, "sed"), 2.0)
        assert np.allclose(rn.get_segment_field(s, "z"), [10.0, 20.0, 30.0])


def test_runs_under_scheduler():
    rn = _chain_with_z()
    U = 2e-3

    def uplift(net, dt):
        return {s: U * np.ones_like(net.get_segment_field(s, "z"))
                for s in net.segment_ids}

    z0 = {s: rn.get_segment_field(s, "z").copy() for s in rn.segment_ids}
    sched = Scheduler(rn, [ExplicitRule(rn, uplift)])
    nt, dt = 4, 100.0
    sched.run(nt, dt)
    for s in rn.segment_ids:
        assert np.allclose(rn.get_segment_field(s, "z"), z0[s] + U * nt * dt)
    assert rn.graph.graph["t"] == nt * dt
