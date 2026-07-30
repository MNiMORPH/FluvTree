"""Non-erodible bedrock: FixedBed returns its bed unchanged."""

import numpy as np

from fluvtree import Scheduler
from fluvtree.processes import FixedBed, build_network


def _network():
    D = 2000.0
    up, down = [[], [], [0, 1]], [[2], [2], []]
    x = [D * np.arange(1, 5.0), D * np.arange(1, 5.0), D * np.arange(5, 9.0)]
    Q = [5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)]
    B = [100.0 * np.ones(4)] * 3
    z = [np.linspace(60, 20, 4), np.linspace(60, 20, 4), np.linspace(20, 5, 4)]
    return build_network(up, down, x, z, Q, B, S0=0.015, x_bl=D * 9, z_bl=0.0)


def test_fixed_bed_returns_its_bed_unchanged():
    rn = _network()
    z0 = {s: rn.get_segment_field(s, "z").copy() for s in rn.segment_ids}
    Scheduler(rn, [FixedBed(rn)]).run(nt=5, dt=1.0e6)
    for s in rn.segment_ids:
        assert np.array_equal(rn.get_segment_field(s, "z"), z0[s])   # bit-for-bit unchanged


def test_fixed_bed_writes_nothing():
    rn = _network()
    assert FixedBed(rn).writes == ()
