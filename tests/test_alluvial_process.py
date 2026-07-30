"""
The physics-named processes drive the in-tree engine through the graph.

  1. GravelBed (in-tree diffusion solver) reproduces GRLP (external grlp)
     to machine precision -- same physics, one engine now lives in FluvTree.
  2. SandBed reproduces stable SRLP end-to-end (graph -> process -> engine).
  3. The Scheduler one-liner runs a sand-bed network (the boxed flow, collapsed).
"""

import importlib.util
import os

import numpy as np
import pytest

from fluvtree import RiverNetwork, Scheduler
from fluvtree.processes import GravelBed, SandBed, build_network


_D = 2000.0

MULTI = dict(
    up=[[], [], [], [], [0, 1], [4, 2], [5, 3]],
    down=[[4], [4], [5], [6], [5], [6], []],
    x=[_D * np.arange(1, 5.0)] * 4
      + [_D * np.arange(5, 9.0), _D * np.arange(9, 13.0), _D * np.arange(13, 17.0)],
    Q=[5 * np.ones(4)] * 4 + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
    x_bl=_D * 17,
)


def _gravel_rn():
    return build_network(
        MULTI["up"], MULTI["down"], MULTI["x"],
        [np.zeros(len(xi)) for xi in MULTI["x"]], MULTI["Q"],
        [100.0 * np.ones(len(xi)) for xi in MULTI["x"]],
        S0=0.015, x_bl=MULTI["x_bl"], z_bl=0.0)


def test_gravel_process_matches_external_grlp_process():
    pytest.importorskip("grlp")
    from fluvtree.processes import GRLP

    # Compare like scheme to like: the in-tree solver is backward Euler (v1), so
    # configure external grlp for backward Euler with the same fixed Picard count.
    def _be(g):
        g.set_time_integration(1)
        g.set_niter(3)

    rn_ext, rn_in = _gravel_rn(), _gravel_rn()
    p_ext = GRLP(rn_ext, configure=_be)
    p_in = GravelBed(rn_in, niter=3)
    for _ in range(6):
        p_ext.step(3.15e10)
        p_in.step(3.15e10)
    for s in rn_in.segment_ids:
        assert np.allclose(rn_ext.get_segment_field(s, "z"),
                           rn_in.get_segment_field(s, "z"), rtol=0, atol=1e-9)


def _load_srlp():
    path = os.path.expanduser("~/models/SRLP/srlp/srlp.py")
    if not os.path.exists(path):
        pytest.skip("SRLP checkout not present")
    spec = importlib.util.spec_from_file_location("srlp", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_sand_process_reproduces_srlp():
    srlp = _load_srlp()
    S0, Bmax, P_xB = 0.003, 250.0, 0.2
    lp = srlp.LongProfile()
    lp.set_niter(3)
    lp.set_D(1e-3); lp.set_Mannings_roughness(0.02)
    lp.set_Darcy_Weisbach_friction(0.1); lp.set_tau_crit_bank(2)
    lp.basic_constants(); lp.sediment_lumped_constants(); lp.set_hydrologic_constants()
    lp.set_x(dx=500, nx=180, x0=10e3); lp.set_z(S0=-S0, z1=0); lp.set_A(k_xA=1.0)
    lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7 / 4.0 * 0.7)
    lp.set_B(k_xB=Bmax / np.max(lp.x ** P_xB), P_xB=P_xB); lp.set_z_bl(0)
    Qs0 = lp.k_Qs * lp.Q[0] * S0 ** (5 / 6.0)
    lp.set_Qs_input_upstream(Qs0); lp.set_uplift_rate(0)
    lp.evolve_threshold_width_river(20, 1e14)
    z_srlp, x, Q, B = lp.z.copy(), lp.x.copy(), lp.Q.copy(), lp.B.copy()

    x_bl = x[-1] + (x[-1] - x[-2])
    rn = build_network(
        [[]], [[]], [x], [S0 * (x_bl - x)], [Q], [B],   # sloped IC
        S0=S0, x_bl=x_bl, z_bl=0.0, Q_s_0=Qs0)          # flux boundary, as SRLP
    proc = SandBed(rn, D=1e-3, n=0.02, tau_crit_bank=2)
    proc.step(1e13, nt=200)

    relief = z_srlp.max() - z_srlp.min()
    assert np.max(np.abs(rn.get_segment_field(0, "z") - z_srlp)) / relief < 1e-3


def test_scheduler_one_liner_runs_sand_network():
    # the boxed flow, collapsed: build graph -> Scheduler([SandBed]).run()
    up, down = [[], [], [0, 1]], [[2], [2], []]
    x = [_D * np.arange(1, 5.0), _D * np.arange(1, 5.0), _D * np.arange(5, 9.0)]
    Q = [5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)]
    B = [200.0 * np.ones(4)] * 3
    S0, x_bl = 1e-3, _D * 9
    z = [S0 * (x_bl - xi) for xi in x]                  # sloped IC (sand)

    rn = build_network(up, down, x, z, Q, B, S0=S0, x_bl=x_bl, z_bl=0.0)
    Scheduler(rn, [SandBed(rn, D=1e-3, n=0.02, tau_crit_bank=2)]).run(nt=300, dt=3e10)

    for s in rn.segment_ids:
        zs = rn.get_segment_field(s, "z")
        assert np.all(np.isfinite(zs))
        assert np.all(np.diff(zs) < 0)                  # descends downstream
    assert rn.graph.graph["t"] == pytest.approx(300 * 3e10)
