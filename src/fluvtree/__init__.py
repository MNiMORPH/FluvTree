"""
FluvTree: a directed convergent graph representation for river-network evolution.

FluvTree holds what is; its methods modify what is. The canonical state of a river
network is a single directed graph (edges = segments/reaches, nodes = junctions);
physical models are *processes* that read named fields off the graph and write them
back. See docs/fluvtree-engine-architecture.md.
"""

from .network import RiverNetwork
from .scheduler import Scheduler

__all__ = ["RiverNetwork", "Scheduler"]
