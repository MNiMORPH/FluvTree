"""
FluvTree processes: physical models that read named fields off the canonical
network graph and write them back. A process owns no state; the graph does.
"""

from fluvtree.processes.grlp import GRLP, build_grlp_network
from fluvtree.processes.alluvial import DiffusionProcess, GravelBed, SandBed
from fluvtree.processes.rule import Rule
from fluvtree.processes.transport_limited import TransportLimitedRate, default_k_Qs

build_network = build_grlp_network

__all__ = ["DiffusionProcess", "GravelBed", "SandBed",
           "Rule", "GRLP", "build_grlp_network", "build_network",
           "TransportLimitedRate", "default_k_Qs"]
