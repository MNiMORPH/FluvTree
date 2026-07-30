"""
The FluvTree model: the single object that holds a river and runs its ruleset.

Pins the fluent front door (``m = ft.FluvTree(net); m.add(ft.GravelBed(D=...));
m.run(...)``) and that it is exactly the Scheduler path underneath -- deferred
binding changes ergonomics, not numerics.
"""

import numpy as np
import pytest

import fluvtree as ft
from fluvtree.processes import build_network


_D = 2000.0


def _confluence():
    up, down = [[], [], [0, 1]], [[2], [2], []]
    x = [_D * np.arange(1, 5.0), _D * np.arange(1, 5.0), _D * np.arange(5, 9.0)]
    Q = [5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)]
    B = [100.0 * np.ones(4)] * 3
    z = [0.015 * (_D * 9 - xi) for xi in x]
    return build_network(up, down, x, z, Q, B, S0=0.015, x_bl=_D * 9, z_bl=0.0)


def test_fluent_equals_scheduler_path():
    # model with a deferred-bound process ...
    net = _confluence()
    m = ft.FluvTree(net)
    assert m.add(ft.GravelBed(D=0.05)) is m           # add returns self (chains)
    assert m.processes[0].network is net              # bound to the model's network
    m.run(nt=200, dt=3.15e10)

    # ... is bit-for-bit the Scheduler path with an eagerly-bound process
    ref = _confluence()
    ft.Scheduler(ref, [ft.GravelBed(ref, D=0.05)]).run(nt=200, dt=3.15e10)
    for s in net.segment_ids:
        assert np.array_equal(net.get_segment_field(s, "z"),
                              ref.get_segment_field(s, "z"))
    assert m.t == pytest.approx(200 * 3.15e10)


def test_run_until_lands_on_target():
    m = ft.FluvTree(_confluence()).add(ft.GravelBed(D=0.05))
    m.run(until=6.3e12, dt=3.15e10)
    assert m.t == pytest.approx(6.3e12)
    # a dt that does not divide the target still lands exactly on it
    m2 = ft.FluvTree(_confluence()).add(ft.GravelBed(D=0.05))
    m2.run(until=1.0e12, dt=3.0e10)
    assert m2.t == pytest.approx(1.0e12)


def test_run_needs_exactly_one_of_nt_until():
    m = ft.FluvTree(_confluence()).add(ft.GravelBed(D=0.05))
    with pytest.raises(ValueError):
        m.run(dt=1.0e10)                       # neither
    with pytest.raises(ValueError):
        m.run(dt=1.0e10, nt=5, until=1.0e12)   # both


def test_unbound_process_steps_only_after_binding():
    proc = ft.GravelBed(D=0.05)                # no network -> unbound
    assert proc.network is None
    with pytest.raises(RuntimeError):
        proc.step(1.0)                         # clear error, not an AttributeError
    proc.bind(_confluence())
    proc.step(3.15e10)                         # now fine


def test_ordered_ruleset_runs_in_order():
    # uplift (a Rule) then gravel transport, composed on one model
    net = _confluence()

    def uplift(network, dt):
        return {s: np.ones(len(network.get_segment_field(s, "z")))
                * 1e-11 for s in network.segment_ids}

    m = ft.FluvTree(net).add(ft.Rule(rate=uplift)).add(ft.GravelBed(D=0.05))
    assert [type(p).__name__ for p in m.processes] == ["Rule", "GravelBed"]
    m.run(nt=50, dt=3.15e10)
    for s in net.segment_ids:
        assert np.all(np.isfinite(net.get_segment_field(s, "z")))
