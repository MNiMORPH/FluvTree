"""
Transport closures: the constitutive laws that close the alluvial governing
equation. ``base`` holds the ``TransportClosure`` interface (the plug socket the
diffusion solver reads); ``gravel`` and ``sand`` are the concrete laws.

Closures are per-solver-form (see docs/DESIGN-structure-and-naming.md): these are
the closure family for the power-law nonlinear-diffusion solver.
"""

from fluvtree.closures.base import TransportClosure
from fluvtree.closures.gravel import GravelClosure
from fluvtree.closures.sand import SandClosure

__all__ = ["TransportClosure", "GravelClosure", "SandClosure"]
