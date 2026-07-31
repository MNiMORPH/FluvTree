"""
Tectonics: uplift / subsidence as a coupled source (the second common module).

Optional, off by default. Uplift is a source in the bed mass balance (not an
operator-split step), so it couples into the implicit solve exactly as GRLP's
``set_uplift_rate`` -- validated bit-for-bit -- and its sign/spatial handling are
pinned here.
"""

import numpy as np
import pytest

from fluvtree.processes import GravelBed, build_network
from fluvtree.common.tectonics import uplift_rate


_D = 2000.0
_U = 3.0e-11           # ~1 mm/yr


def _confluence(slope=0.015):
    up = [[], [], [], [], [0, 1], [4, 2], [5, 3]]
    down = [[4], [4], [5], [6], [5], [6], []]
    x = ([_D * np.arange(1, 5.0)] * 4
         + [_D * np.arange(5, 9.0), _D * np.arange(9, 13.0), _D * np.arange(13, 17.0)])
    Q = [5 * np.ones(4)] * 4 + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)]
    z = [slope * (_D * 17 - xi) for xi in x]
    B = [100.0 * np.ones(len(xi)) for xi in x]
    return build_network(up, down, x, z, B=B, Q=Q, S0=0.015, x_bl=_D * 17, z_bl=0.0)


def test_uplift_rate_scalar_per_segment_and_callable():
    rn = _confluence()
    segs = list(rn.segment_ids)
    # scalar -> uniform on every node
    u = uplift_rate(rn, _U)
    assert all(np.all(u[i] == _U) for i in range(len(segs)))
    assert u[0].shape == rn.get_segment_field(segs[0], "x").shape
    # per-segment sequence
    per = [float(i) for i in range(len(segs))]
    u2 = uplift_rate(rn, per)
    assert np.all(u2[3] == 3.0)
    # callable of x (spatial pattern)
    u3 = uplift_rate(rn, lambda x: 1e-11 * x)
    x0 = rn.get_segment_field(segs[0], "x")
    assert np.allclose(u3[0], 1e-11 * x0)


def test_uplift_matches_grlp():
    pytest.importorskip("grlp")
    from fluvtree.processes import GRLP

    def _cfg(g):
        g.set_time_integration(2)
        g.set_niter(3)
        for seg in g.segments:
            seg.set_uplift_rate(_U)

    rn_ext, rn_new = _confluence(), _confluence()
    GRLP(rn_ext, configure=_cfg).step(3.15e10, nt=6)
    GravelBed(rn_new, niter=3, uplift=_U).step(3.15e10, nt=6)
    for s in rn_new.segment_ids:
        assert np.allclose(rn_ext.get_segment_field(s, "z"),
                           rn_new.get_segment_field(s, "z"), rtol=0, atol=1e-9)


def test_uplift_raises_subsidence_lowers_and_off_by_default():
    # off by default
    rn_none, rn_default = _confluence(), _confluence()
    GravelBed(rn_none, niter=3, uplift=None).step(3.15e10, nt=6)
    GravelBed(rn_default, niter=3).step(3.15e10, nt=6)
    for s in rn_none.segment_ids:
        assert np.array_equal(rn_none.get_segment_field(s, "z"),
                              rn_default.get_segment_field(s, "z"))

    # +uplift raises the bed, -subsidence lowers it, relative to no forcing
    rn_up, rn_sub = _confluence(), _confluence()
    GravelBed(rn_up, niter=3, uplift=_U).step(3.15e10, nt=6)
    GravelBed(rn_sub, niter=3, uplift=-_U).step(3.15e10, nt=6)
    base = np.mean(rn_default.get_segment_field(0, "z"))
    assert np.mean(rn_up.get_segment_field(0, "z")) > base
    assert np.mean(rn_sub.get_segment_field(0, "z")) < base
