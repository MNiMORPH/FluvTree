"""
FluvTree processes: physical models that read named fields off the canonical
network graph and write them back. A process owns no state; the graph does.

Named by the specific law/mechanism, not the domain family (GravelBed/SandBed are
transport-closure presets; StreamPower is a specific incision law; FixedBed a null
behavior; Rule an arbitrary supplied dz/dt). See docs/DESIGN-structure-and-naming.md.
"""

from fluvtree.processes.base import Process
from fluvtree.processes.grlp import GRLP, build_grlp_network
from fluvtree.processes.alluvial import DiffusionProcess, GravelBed, SandBed
from fluvtree.processes.fixed_bed import FixedBed
from fluvtree.processes.stream_power import StreamPower
from fluvtree.processes.rule import Rule
from fluvtree.processes.transport_limited import TransportLimitedRate, default_k_Qs

# Physics-neutral name for the graph builder (it just stamps x/z/Q/B + boundaries
# onto a RiverNetwork; not GRLP-specific).
build_network = build_grlp_network

__all__ = ["Process", "DiffusionProcess", "GravelBed", "SandBed",
           "StreamPower", "FixedBed", "Rule",
           "GRLP", "build_grlp_network", "build_network",
           "TransportLimitedRate", "default_k_Qs"]
