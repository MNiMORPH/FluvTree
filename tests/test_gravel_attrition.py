"""
Gravel attrition (Sternberg downstream fining) -- the first common module.

An optional, opt-in sink: the transported gravel abrades and fines downstream,
leaving the bedload. These tests pin (1) the ``compute_Q_s`` kernel against GRLP's
own network sediment-flux walk, (2) the full sink bit-for-bit against GRLP's
``update_gravel_loss`` run, and (3) that it is genuinely a sink and off by default.
"""

import numpy as np
import pytest

from fluvtree.processes import GravelBed, build_network
from fluvtree.closures.gravel import GravelClosure
from fluvtree.common.gravel_attrition import compute_Q_s


_D = 2000.0
_COEFF = 0.05          # fractional gravel-load loss per km


def _confluence(slope=0.03):
    up = [[], [], [], [], [0, 1], [4, 2], [5, 3]]
    down = [[4], [4], [5], [6], [5], [6], []]
    x = ([_D * np.arange(1, 5.0)] * 4
         + [_D * np.arange(5, 9.0), _D * np.arange(9, 13.0), _D * np.arange(13, 17.0)])
    Q = [5 * np.ones(4)] * 4 + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)]
    z = [slope * (_D * 17 - xi) for xi in x]
    B = [100.0 * np.ones(len(xi)) for xi in x]
    return build_network(up, down, x, z, B=B, Q=Q, S0=0.015, x_bl=_D * 17, z_bl=0.0)


def test_compute_Q_s_matches_grlp():
    pytest.importorskip("grlp")
    from fluvtree.processes import GRLP

    rn = _confluence()
    g = GRLP(rn)                            # internal grlp network, same z
    g.grlp_network.compute_Q_s()
    ours = compute_Q_s(rn, GravelClosure(),
                       [rn.get_segment_field(s, "z") for s in rn.segment_ids])
    for i, s in enumerate(rn.segment_ids):
        assert np.allclose(g.grlp_network.segments[i].Q_s, ours[i], rtol=0, atol=1e-12)


def test_attrition_matches_grlp_sternberg():
    pytest.importorskip("grlp")
    from fluvtree.processes import GRLP

    def _cfg(g):
        g.set_time_integration(2)          # BDF2 (default)
        g.set_niter(3)
        for seg in g.segments:
            seg.set_Sternberg_gravel_loss(_COEFF)

    rn_ext, rn_new = _confluence(), _confluence()
    GRLP(rn_ext, configure=_cfg).step(3.15e10, nt=6)
    GravelBed(rn_new, niter=3, gravel_attrition=_COEFF).step(3.15e10, nt=6)
    for s in rn_new.segment_ids:
        assert np.allclose(rn_ext.get_segment_field(s, "z"),
                           rn_new.get_segment_field(s, "z"), rtol=0, atol=1e-9)


def test_attrition_is_active_and_off_by_default():
    # off by default: gravel_attrition=None is identical to not passing it
    rn_none, rn_default = _confluence(), _confluence()
    GravelBed(rn_none, niter=3, gravel_attrition=None).step(3.15e10, nt=6)
    GravelBed(rn_default, niter=3).step(3.15e10, nt=6)
    for s in rn_none.segment_ids:
        assert np.array_equal(rn_none.get_segment_field(s, "z"),
                              rn_default.get_segment_field(s, "z"))

    # turning it on changes the profile (it is a real sink, not a no-op)
    rn_on = _confluence()
    GravelBed(rn_on, niter=3, gravel_attrition=_COEFF).step(3.15e10, nt=6)
    diff = max(float(np.max(np.abs(rn_on.get_segment_field(s, "z")
                                   - rn_default.get_segment_field(s, "z"))))
               for s in rn_on.segment_ids)
    assert diff > 1.0
