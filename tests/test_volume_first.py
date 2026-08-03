"""
Volume-first solve: the diffusion solver conserves stored sediment volume ``V`` and
recovers ``z`` through a valley-storage geometry.

Two things to pin:
- the **rectangular default reproduces the constant-B model** -- a constant-width
  ``ValleyGeometry`` must give the same profile as the default rectangular valley,
  bit-for-bit (the transform is an exact row-scaling there);
- a **z-varying valley genuinely engages** -- a valley that widens with height must
  change the transient (else the geometry is being silently ignored).

The bit-for-bit-vs-GRLP tests elsewhere already exercise the transform on the
constant-B path; these fix the geometry plumbing itself.
"""

import numpy as np

from fluvtree.processes import build_network
from fluvtree.closures.gravel import GravelClosure
from fluvtree.solvers.diffusion import DiffusionSolver
from fluvtree.valley import ValleyGeometry

_D = 2000.0


def _chain():
    x = _D * np.arange(1, 13.0)
    z = 0.03 * (_D * 13 - x)                    # steeper-than-steady transient IC
    Q = 10.0 * np.ones(12)
    B = 100.0 * np.ones(12)
    return build_network([[]], [[]], [x], [z], [Q], [B],
                         S0=0.015, x_bl=_D * 13, z_bl=0.0)


def test_constant_geometry_matches_rectangular_default():
    clo = GravelClosure(D=0.05)
    rn_def, rn_tab = _chain(), _chain()
    # explicit constant-B valley spanning the elevation range
    zlev = [np.array([-50.0, 2000.0])] * 12
    Blev = [np.array([100.0, 100.0])] * 12
    vg = ValleyGeometry(zlev, Blev, clo.lambda_p)
    DiffusionSolver(rn_def, clo).evolve(nt=6, dt=3.15e10, niter=3)
    DiffusionSolver(rn_tab, clo, geometry=vg).evolve(nt=6, dt=3.15e10, niter=3)
    assert np.allclose(rn_def.get_segment_field(0, "z"),
                       rn_tab.get_segment_field(0, "z"), rtol=0, atol=1e-9)


def test_widening_valley_changes_the_transient():
    clo = GravelClosure(D=0.05)
    rn_const, rn_wide = _chain(), _chain()
    zlev = [np.array([-50.0, 2000.0])] * 12
    B_const = [np.array([100.0, 100.0])] * 12
    B_wide = [np.array([100.0, 100.0 + 0.5 * (2000.0 + 50.0)])] * 12   # widens with z
    DiffusionSolver(rn_const, clo,
                    geometry=ValleyGeometry(zlev, B_const, clo.lambda_p)
                    ).evolve(nt=6, dt=3.15e10, niter=3)
    DiffusionSolver(rn_wide, clo,
                    geometry=ValleyGeometry(zlev, B_wide, clo.lambda_p)
                    ).evolve(nt=6, dt=3.15e10, niter=3)
    z_const = rn_const.get_segment_field(0, "z")
    z_wide = rn_wide.get_segment_field(0, "z")
    assert np.all(np.isfinite(z_wide))
    assert np.max(np.abs(z_wide - z_const)) > 0.1     # the geometry actually engages
