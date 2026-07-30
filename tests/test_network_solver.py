"""
The power-law nonlinear-diffusion solver, on the shared network.

``DiffusionSolver`` walks the canonical :class:`RiverNetwork` -- its own
segment fields and topology -- and a transport closure, building no parallel
engine object. These tests pin that it is correct two independent ways:

  1. **Bit-for-bit vs external GRLP** in backward-Euler mode (fixed Picard count),
     on a single chain, a 1-into-1 series, and a multi-tributary confluence: the
     lifted interior stencil, head/outlet boundaries, and conservative confluence
     junction cell reproduce GRLP to float round-off. Needs ``grlp`` (CI skips).
  2. **Analytic steady state**, no external dependency: a network whose channel
     heads' summed sediment supply equals the outlet discharge relaxes to a
     uniform transport-limited slope ``S = S0``. This is the scheme-independent
     physics fixed point and runs without ``grlp``.
"""

import numpy as np
import pytest

from fluvtree.processes import build_network
from fluvtree.solvers.diffusion import DiffusionSolver
from fluvtree.closures.gravel import GravelClosure
from fluvtree.closures.sand import SandClosure


_D = 2000.0
NITER = 3


def _be_config(g):
    """Configure external grlp for backward Euler with a fixed Picard count."""
    g.set_time_integration(1)
    g.set_niter(NITER)


def _chain():
    x = [_D * np.arange(1, 9.0)]
    z = [0.015 * (_D * 9 - x[0])]
    return build_network([[]], [[]], x, z, [10 * np.ones(8)], [100.0 * np.ones(8)],
                         S0=0.015, x_bl=_D * 9, z_bl=0.0)


def _series():
    up, down = [[], [0]], [[1], []]
    x = [_D * np.arange(1, 5.0), _D * np.arange(5, 9.0)]
    z = [0.015 * (_D * 9 - xi) for xi in x]
    Q = [10 * np.ones(4), 10 * np.ones(4)]
    B = [100.0 * np.ones(4), 100.0 * np.ones(4)]
    return build_network(up, down, x, z, Q, B, S0=0.015, x_bl=_D * 9, z_bl=0.0)


def _confluence():
    up = [[], [], [], [], [0, 1], [4, 2], [5, 3]]
    down = [[4], [4], [5], [6], [5], [6], []]
    x = ([_D * np.arange(1, 5.0)] * 4
         + [_D * np.arange(5, 9.0), _D * np.arange(9, 13.0), _D * np.arange(13, 17.0)])
    Q = [5 * np.ones(4)] * 4 + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)]
    z = [0.015 * (_D * 17 - xi) for xi in x]
    B = [100.0 * np.ones(len(xi)) for xi in x]
    return build_network(up, down, x, z, Q, B, S0=0.015, x_bl=_D * 17, z_bl=0.0)


@pytest.mark.parametrize("maker", [_chain, _series, _confluence],
                         ids=["chain", "series", "confluence"])
def test_backward_euler_matches_external_grlp(maker):
    pytest.importorskip("grlp")
    from fluvtree.processes import GRLP

    rn_ext, rn_new = maker(), maker()
    p_ext = GRLP(rn_ext, configure=_be_config)
    solver = DiffusionSolver(rn_new, GravelClosure())
    for _ in range(6):
        p_ext.step(3.15e10)
        solver.evolve(nt=1, dt=3.15e10, niter=NITER)
    for s in rn_new.segment_ids:
        assert np.allclose(rn_ext.get_segment_field(s, "z"),
                           rn_new.get_segment_field(s, "z"), rtol=0, atol=1e-9)


def test_steady_state_uniform_slope():
    # 4 heads (Q=5, S0=0.015) sum to the outlet discharge (20), so the
    # transport-limited steady state is a uniform slope S = S0 everywhere.
    rn = _confluence()
    solver = DiffusionSolver(rn, GravelClosure())
    solver.evolve(nt=400, dt=3.15e11, tol=1e-4)
    for s in rn.segment_ids:
        z = rn.get_segment_field(s, "z")
        x = rn.get_segment_field(s, "x")
        S = -np.diff(z) / np.diff(x)
        assert np.allclose(S, 0.015, rtol=2e-3)


def test_sand_steady_state_obeys_sand_flux_law():
    """The closure's exponent is really carried through the solver: a sand
    (``p = 5/6``) steady state obeys the sand flux law ``S = (Qs/(k_Qs Q))**(1/p)``
    -- and is far off the gravel exponent, so the solver is not secretly gravel."""
    x = _D * np.arange(1, 9.0)
    Q = np.linspace(30.0, 60.0, 8)
    B = 200.0 * np.ones(8)
    S0 = 1.0e-3
    sc = SandClosure(D=0.3e-3, n=0.03, tau_crit_bank=5.0)

    rn = build_network([[]], [[]], x=[x.copy()], z=[S0 * (_D * 9 - x)],  # sloped IC
                       Q=[Q.copy()], B=[B.copy()],
                       S0=S0, x_bl=_D * 9, z_bl=0.0)
    solver = DiffusionSolver(rn, sc)
    solver.evolve(nt=500, dt=3.0e10, tol=1e-4)

    z = rn.get_segment_field(0, "z")
    Sf = np.abs(np.diff(z) / np.diff(x))
    Qf = 0.5 * (Q[:-1] + Q[1:])
    supply = sc.k_Qs * Q[0] * S0 ** sc.p
    qs = sc.k_Qs * Qf * Sf ** sc.p                     # conserved during-flood flux
    assert np.max(np.abs(qs - supply)) / supply < 0.01  # coarse 8-node grid -> ~1%
    S_sand = (supply / (sc.k_Qs * Qf)) ** (1 / sc.p)
    S_gravel = (supply / (sc.k_Qs * Qf)) ** (6 / 7.0)
    assert np.max(np.abs(Sf - S_sand) / Sf) < 0.02
    assert np.max(np.abs(Sf - S_gravel) / Sf) > 1.0
