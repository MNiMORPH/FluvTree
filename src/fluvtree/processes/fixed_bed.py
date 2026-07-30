"""
Bedrock-river processes.

Two ends of a spectrum:
  - **erodible** bedrock -- stream-power incision, ``dz/dt = U - K A^m S^n``
    (advective, its own core; future work);
  - **non-erodible** bedrock -- :class:`FixedBed`, the bed returned
    unchanged.

``FixedBed`` is the rigid endpoint: no erosion, no deposition, no sediment
moved -- **conservative by construction**. Note that in the current whole-network
solver it is a *marker* of intent (it declares which reaches are bedrock), not yet
an *enforcer*: making a reach truly non-erodible while a co-running alluvial process
evolves the network needs the solver to hold those reaches. That enforcement --
together with a smooth alluvial->bedrock taper and the sediment mass accounting it
requires -- is deferred; this rigid bed is the anchor that taper will taper toward.
(Overwriting ``z`` back to a fixed value each step would enforce it, but that is
nonconservative -- it fabricates/destroys the sediment the erosion moved -- which is
exactly why this process moves nothing.)
"""


class FixedBed(object):
    """
    Non-erodible bedrock: returns its bed unchanged (no erosion possible).

    Parameters
    ----------
    network : RiverNetwork
        The canonical network.
    segments : iterable, optional
        The bedrock reaches (default: the whole network).
    """

    reads = ()
    writes = ()

    def __init__(self, network, segments=None):
        self.network = network
        self.segments = (list(network.segment_ids) if segments is None
                         else list(segments))

    def step(self, dt):
        """No erosion possible: the bed is returned unchanged."""
        pass
