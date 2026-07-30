"""
Stream-power incision as a FluvTree process (erodible bedrock, n = 1).

The erodible end of the bedrock family (the counterpart to
:class:`FixedBed`). Drives ``fluvtree.solvers.advection``'s implicit
outlet->upstream sweep on the canonical graph -- no engine object to build, since
the sweep mutates the graph's ``z`` arrays directly.

"Stream power" here is detachment-limited bedrock incision, ``dz/dt = U - K Q**m S``
-- structurally different from the alluvial (transport-limited) processes: advection,
not diffusion.
"""

from fluvtree.processes.base import Process
from fluvtree.solvers.advection import sweep_order, incise_n1_step


class StreamPower(Process):
    """
    n = 1 stream-power bedrock incision.

    Parameters
    ----------
    network : RiverNetwork, optional
        Canonical network carrying ``x, z, Q`` (edges) and base level ``x_bl``,
        ``z_bl`` (outlet node). Omit to bind later via a :class:`FluvTree` model or
        :class:`Scheduler`.
    K : float
        Erodibility.
    m : float, optional
        Discharge exponent (default 0.5).
    U : float, optional
        Uplift rate (scalar, default 0.0). Per-segment uplift via a graph field is
        a later extension.
    """

    reads = ("x", "z", "Q", "x_bl", "z_bl")
    writes = ("z",)

    def __init__(self, network=None, K=None, m=0.5, U=0.0):
        self.K = K
        self.m = m
        self.U = U
        self._order = None
        super().__init__(network)

    def _on_bind(self):
        self._order = sweep_order(self.network)     # fixed topology: computed once

    def step(self, dt, nt=1):
        """Advance ``nt`` implicit steps of ``dt`` (in place on the graph's z)."""
        self._require_bound()
        for _ in range(int(nt)):
            incise_n1_step(self.network, self._order, self.K, self.m, self.U, dt)
