"""
FluvTree generate: *synthesis* of networks and inputs from parameters.

One of the three ways an input network arrives (see docs/DESIGN-structure-and-naming.md):

- **synthesis** -- invented from parameters (this package): arbitrary/random network
  topologies, idealized profiles;
- *explicit* -- caller-provided values -> ``RiverNetwork.from_segment_lists`` /
  ``RiverNetwork.from_arrays``;
- *loading* -- from data (DEM, measured) -> a future ``io`` layer.

``ShreveRandomNetwork`` generates a random drainage topology (Shreve 1974) and
``random_network`` returns it as a :class:`~fluvtree.network.RiverNetwork`
(topology; attach fields with ``from_arrays`` or the planned field generators).
"""

from fluvtree.generate.topology import ShreveRandomNetwork, random_network

__all__ = ["ShreveRandomNetwork", "random_network"]
