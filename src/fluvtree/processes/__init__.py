"""
FluvTree processes: physical models that read named fields off the canonical
network graph and write them back. A process owns no state; the graph does.

See docs/fluvtree-engine-architecture.md.
"""

from .grlp_process import GRLPProcess, build_grlp_network
from .explicit_rule import ExplicitRule

__all__ = ["GRLPProcess", "build_grlp_network", "ExplicitRule"]
