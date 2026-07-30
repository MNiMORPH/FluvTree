"""
Benchmark the diffusion solver's sand path against stable SRLP.

Mirrors the gravel-vs-GRLP benchmark: drive the *same* single-channel setup as
SRLP's ``Simple_Qs`` example through (1) stable SRLP itself and (2) the
the diffusion solver :class:`DiffusionSolver` with a ``SandClosure`` on a
single-segment network, and require the steady profiles to agree. They agree to
~0.02% of relief -- the residual is the difference between SRLP's older
tridiagonal engine and the walking assembler, not physics; the sand transport law
is reproduced exactly (matching ``k_Qs`` and ``p = 5/6``).

Skipped if the SRLP checkout is not present (loaded from a sibling path, since
SRLP is not pip-installed).
"""

import importlib.util
import os

import numpy as np
import pytest

from fluvtree.processes import build_network
from fluvtree.solvers.diffusion import DiffusionSolver
from fluvtree.closures.sand import SandClosure


def _load_srlp():
    path = os.path.expanduser("~/models/SRLP/srlp/srlp.py")
    if not os.path.exists(path):
        pytest.skip("SRLP checkout not present")
    spec = importlib.util.spec_from_file_location("srlp", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_sand_engine_reproduces_srlp_steady_state():
    srlp = _load_srlp()
    S0, Bmax, P_xB = 0.003, 250.0, 0.2

    # --- SRLP oracle (its Simple_Qs single-channel example, run to steady) ---
    lp = srlp.LongProfile()
    lp.set_niter(3)
    lp.set_D(1e-3)
    lp.set_Mannings_roughness(0.02)
    lp.set_Darcy_Weisbach_friction(0.1)
    lp.set_tau_crit_bank(2)
    lp.basic_constants()
    lp.sediment_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=500, nx=180, x0=10e3)
    lp.set_z(S0=-S0, z1=0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7 / 4.0 * 0.7)
    lp.set_B(k_xB=Bmax / np.max(lp.x ** P_xB), P_xB=P_xB)
    lp.set_z_bl(0)
    Qs0 = lp.k_Qs * lp.Q[0] * S0 ** (5 / 6.0)
    lp.set_Qs_input_upstream(Qs0)
    lp.set_uplift_rate(0)
    lp.evolve_threshold_width_river(20, 1e14)
    z_srlp = lp.z.copy()
    x, Q, B = lp.x.copy(), lp.Q.copy(), lp.B.copy()

    # --- diffusion solver + SandClosure, SRLP's own geometry ---
    x_bl = x[-1] + (x[-1] - x[-2])
    closure = SandClosure(D=1e-3, n=0.02, tau_crit_bank=2)
    rn = build_network(
        [[]], [[]],
        x=[x.copy()], z=[S0 * (x_bl - x)],   # sloped IC (sand conductance ~ S**-1/6)
        Q=[Q.copy()], B=[B.copy()],
        S0=S0, x_bl=x_bl, z_bl=0.0, Q_s_0=Qs0)
    solver = DiffusionSolver(rn, closure)
    solver.evolve(nt=200, dt=1e13, tol=1e-4)
    z_new = rn.get_segment_field(0, "z")

    # exact transport law, near-exact profile (old vs walking assembler only)
    assert np.isclose(closure.k_Qs, lp.k_Qs, rtol=0, atol=0)
    relief = z_srlp.max() - z_srlp.min()
    assert np.max(np.abs(z_new - z_srlp)) / relief < 1e-3
