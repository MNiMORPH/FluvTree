"""
Solvers: the numerical engines, named by their mathematical form.

- ``diffusion`` -- power-law nonlinear diffusion (transport-limited long profile);
  the implicit ``DiffusionSolver`` on the shared network, parameterized by a
  transport closure.
- ``advection`` -- power-law nonlinear advection (stream-power incision); the
  outlet->upstream sweep. ``AdvectionSolver`` wraps it; the ``n = 1`` linear rung
  is implemented.

See docs/DESIGN-structure-and-naming.md.
"""

from fluvtree.solvers.diffusion import DiffusionSolver
from fluvtree.solvers.advection import (
    AdvectionSolver, sweep_order, incise_n1_step, evolve_streampower_n1,
)

__all__ = ["DiffusionSolver", "AdvectionSolver",
           "sweep_order", "incise_n1_step", "evolve_streampower_n1"]
