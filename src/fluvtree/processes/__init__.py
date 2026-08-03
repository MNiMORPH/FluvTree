"""
FluvTree processes: physical models that read named fields off the canonical
network graph and write them back. A process owns no state; the graph does.

Named by the specific law/mechanism, not the domain family (GravelBed/SandBed are
transport-closure presets; StreamPower is a specific incision law; FixedBed a null
behavior; Rule an arbitrary supplied dz/dt). See docs/DESIGN-structure-and-naming.md.
"""

from fluvtree.network import RiverNetwork
from fluvtree.processes.base import Process
from fluvtree.processes.grlp import GRLP
from fluvtree.processes.alluvial import DiffusionProcess, GravelBed, SandBed
from fluvtree.processes.fixed_bed import FixedBed
from fluvtree.processes.stream_power import StreamPower
from fluvtree.processes.rule import Rule
from fluvtree.processes.transport_limited import TransportLimitedRate, default_k_Qs

# Convenience re-export of the explicit network constructor. Its canonical home is
# RiverNetwork.from_arrays (it stamps x/z/Q/B + boundaries onto a network; it is not
# a process). Kept here because much of the suite builds test networks with it.
build_network = RiverNetwork.from_arrays

__all__ = ["Process", "DiffusionProcess", "GravelBed", "SandBed",
           "StreamPower", "FixedBed", "Rule",
           "GRLP", "build_network",
           "TransportLimitedRate", "default_k_Qs"]
