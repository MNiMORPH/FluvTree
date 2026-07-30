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

from fluvtree.processes.base import Process


class FixedBed(Process):
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

    def __init__(self, network=None, segments=None):
        self._segments_arg = segments
        self.segments = None
        super().__init__(network)

    def _on_bind(self):
        self.segments = (list(self.network.segment_ids)
                         if self._segments_arg is None else list(self._segments_arg))

    def step(self, dt, nt=1):
        """No erosion possible: the bed is returned unchanged."""
        pass
