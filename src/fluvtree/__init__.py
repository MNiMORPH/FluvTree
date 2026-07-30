"""
FluvTree: a directed convergent graph framework for river-network evolution.

``FluvTree`` holds what is; its methods modify what is. The canonical state of a
river network is a single directed graph -- the :class:`RiverNetwork` (edges =
segments/reaches, nodes = junctions) -- and physical models are *processes* that
read named fields off it and write them back.

This is the single public front door: ``import fluvtree as ft`` exposes the
structure (:class:`RiverNetwork`), the physics you attach (``GravelBed``,
``SandBed``, ``StreamPower``, ``FixedBed``, ``Rule``), the external cross-check
adapter (``GRLP``), and the scheduler that runs an ordered ruleset. The internal
engines live in the ``solvers`` / ``closures`` subpackages.

See docs/DESIGN-structure-and-naming.md and docs/fluvtree-engine-architecture.md.
"""

from fluvtree.network import RiverNetwork
from fluvtree.model import FluvTree
from fluvtree.scheduler import Scheduler
from fluvtree.processes import (
    GravelBed, SandBed, StreamPower, FixedBed, Rule, GRLP,
    DiffusionProcess, Process, build_network,
)

__all__ = ["FluvTree", "RiverNetwork", "Scheduler",
           "GravelBed", "SandBed", "StreamPower", "FixedBed", "Rule", "GRLP",
           "DiffusionProcess", "Process", "build_network"]
