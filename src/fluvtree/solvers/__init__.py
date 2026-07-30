"""
Solvers: the numerical engines, named by their mathematical form.

- ``diffusion`` -- power-law nonlinear diffusion (transport-limited long profile);
  the implicit ``DiffusionSolver`` on the shared network, parameterized by a
  transport closure.

See docs/DESIGN-structure-and-naming.md.
"""

from fluvtree.solvers.diffusion import DiffusionSolver

__all__ = ["DiffusionSolver"]
