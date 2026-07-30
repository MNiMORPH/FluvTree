"""
FluvTree processes: physical models that read named fields off the canonical
network graph and write them back. A process owns no state; the graph does.
"""

from fluvtree.processes.grlp import GRLP, build_grlp_network
from fluvtree.processes.rule import Rule
from fluvtree.processes.transport_limited import TransportLimitedRate, default_k_Qs

# Physics-neutral name for the graph builder.
build_network = build_grlp_network

__all__ = ["GRLP", "build_grlp_network", "build_network",
           "Rule", "TransportLimitedRate", "default_k_Qs"]
