"""
Generic distributed source/sink -- the shared ``source_sink`` channel.

A signed per-node sediment source/sink rate [m/s] (GRLP's
``set_source_sink_distributed`` -- named ``source_sink`` here to avoid the
"suspended sediment discharge" reading of "ssd"). It is generic: the same channel
for gravel and sand (GRLP and SRLP contexts), and it composes additively with
uplift and attrition. Validated bit-for-bit against GRLP.
"""

import numpy as np
import pytest

from fluvtree.processes import GravelBed, SandBed, build_network


_D = 2000.0
_SS = 5.0e-11          # distributed sediment source rate [m/s]
_U = 3.0e-11


def _confluence(slope=0.015):
    up = [[], [], [], [], [0, 1], [4, 2], [5, 3]]
    down = [[4], [4], [5], [6], [5], [6], []]
    x = ([_D * np.arange(1, 5.0)] * 4
         + [_D * np.arange(5, 9.0), _D * np.arange(9, 13.0), _D * np.arange(13, 17.0)])
    Q = [5 * np.ones(4)] * 4 + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)]
    z = [slope * (_D * 17 - xi) for xi in x]
    B = [100.0 * np.ones(len(xi)) for xi in x]
    return build_network(up, down, x, z, B=B, Q=Q, S0=0.015, x_bl=_D * 17, z_bl=0.0)


def test_source_sink_matches_grlp():
    pytest.importorskip("grlp")
    from fluvtree.processes import GRLP

    def _cfg(g):
        g.set_time_integration(2)
        g.set_niter(3)
        for seg in g.segments:
            seg.set_source_sink_distributed(_SS)

    rn_ext, rn_new = _confluence(), _confluence()
    GRLP(rn_ext, configure=_cfg).step(3.15e10, nt=6)
    GravelBed(rn_new, niter=3, source_sink=_SS).step(3.15e10, nt=6)
    for s in rn_new.segment_ids:
        assert np.allclose(rn_ext.get_segment_field(s, "z"),
                           rn_new.get_segment_field(s, "z"), rtol=0, atol=1e-9)


def test_source_sink_composes_with_uplift():
    pytest.importorskip("grlp")
    from fluvtree.processes import GRLP

    def _cfg(g):
        g.set_time_integration(2)
        g.set_niter(3)
        for seg in g.segments:
            seg.set_source_sink_distributed(_SS)
            seg.set_uplift_rate(_U)

    rn_ext, rn_new = _confluence(), _confluence()
    GRLP(rn_ext, configure=_cfg).step(3.15e10, nt=6)
    GravelBed(rn_new, niter=3, source_sink=_SS, uplift=_U).step(3.15e10, nt=6)
    for s in rn_new.segment_ids:
        assert np.allclose(rn_ext.get_segment_field(s, "z"),
                           rn_new.get_segment_field(s, "z"), rtol=0, atol=1e-9)


def test_source_raises_sink_lowers_and_off_by_default():
    rn_none, rn_default = _confluence(), _confluence()
    GravelBed(rn_none, niter=3, source_sink=None).step(3.15e10, nt=6)
    GravelBed(rn_default, niter=3).step(3.15e10, nt=6)
    for s in rn_none.segment_ids:
        assert np.array_equal(rn_none.get_segment_field(s, "z"),
                              rn_default.get_segment_field(s, "z"))

    rn_src, rn_snk = _confluence(), _confluence()
    GravelBed(rn_src, niter=3, source_sink=_SS).step(3.15e10, nt=6)     # +source
    GravelBed(rn_snk, niter=3, source_sink=-_SS).step(3.15e10, nt=6)    # -sink
    base = np.mean(rn_default.get_segment_field(0, "z"))
    assert np.mean(rn_src.get_segment_field(0, "z")) > base
    assert np.mean(rn_snk.get_segment_field(0, "z")) < base


def test_source_sink_works_for_sand():
    # same generic channel serves the sand (SRLP) context; smoke: it runs and acts
    x = _D * np.arange(1, 13.0)
    z = 1.0e-3 * (_D * 13 - x)                     # sloped IC (sand singular at S=0)
    rn = build_network([[]], [[]], [x], [z], [10 * np.ones(12)], [200.0 * np.ones(12)],
                       S0=1.0e-3, x_bl=_D * 13, z_bl=0.0)
    rn_off = build_network([[]], [[]], [x], [z.copy()], [10 * np.ones(12)],
                           [200.0 * np.ones(12)], S0=1.0e-3, x_bl=_D * 13, z_bl=0.0)
    SandBed(rn, D=1e-3, n=0.02, tau_crit_bank=2, niter=3, source_sink=_SS).step(3e10, nt=6)
    SandBed(rn_off, D=1e-3, n=0.02, tau_crit_bank=2, niter=3).step(3e10, nt=6)
    assert np.all(np.isfinite(rn.get_segment_field(0, "z")))
    assert np.max(np.abs(rn.get_segment_field(0, "z")
                         - rn_off.get_segment_field(0, "z"))) > 0
