"""
Stream-power incision as a FluvTree process (erodible bedrock, n = 1).

The erodible end of the bedrock family (the counterpart to
:class:`FixedBed`). Drives ``fluvtree.solvers.advection``'s implicit
outlet->upstream sweep on the canonical graph -- no engine object to build, since
the sweep mutates the graph's ``z`` arrays directly.

"Stream power" here is detachment-limited bedrock incision, ``dz/dt = U - K Q**m S``
-- structurally different from the alluvial (transport-limited) processes: advection,
not diffusion.

Vertical forcing is specified the same way as on the alluvial side -- ``uplift``
and ``source_sink``, each a scalar, a per-segment value, or a callable of ``x``,
normalized by :func:`fluvtree.common.rates.broadcast_rate` -- so one forcing
statement drives either door. The caveat is physical rather than mechanical: a
*sediment* source or sink has no meaning for detachment-limited incision, which
carries no sediment budget, so ``source_sink`` here should only ever carry
vertical bed motion.
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
        Uplift rate as a bare scalar (default 0.0), the original interface. Kept
        for compatibility; ``uplift`` is the general form and the two sum.
    uplift : scalar, per-segment, or callable, optional
        Vertical rock motion [m/s], positive up. Accepts whatever
        :func:`fluvtree.common.rates.broadcast_rate` accepts, so a spatial
        pattern can be given as ``uplift=lambda x: ...`` without precomputing.
    source_sink : scalar, per-segment, or callable, optional
        A further vertical rate [m/s], the same shape as ``uplift``. See the
        module docstring on why only vertical motion belongs here.
    """

    reads = ("x", "z", "Q", "x_bl", "z_bl")
    writes = ("z",)

    def __init__(self, network=None, K=None, m=0.5, U=0.0, uplift=None,
                 source_sink=None):
        self.K = K
        self.m = m
        self.U = U
        self._uplift = uplift
        self._source_sink = source_sink
        self._order = None
        super().__init__(network)

    def _on_bind(self):
        self._order = sweep_order(self.network)     # fixed topology: computed once

    def _static_source_rate(self):
        """Sum the static per-node vertical forcings into one ``source_rate`` for
        the solver, or ``None`` if neither is set. Mirrors
        :meth:`fluvtree.processes.alluvial.Alluvial._static_source_rate`."""
        from fluvtree.common.rates import broadcast_rate
        nseg = len(self.network.segment_ids)
        total = None
        for spec in (self._uplift, self._source_sink):
            if spec is None:
                continue
            contrib = broadcast_rate(self.network, spec)
            if total is None:
                total = [c.copy() for c in contrib]
            else:
                total = [total[si] + contrib[si] for si in range(nseg)]
        return total

    def step(self, dt, nt=1):
        """Advance ``nt`` implicit steps of ``dt`` (in place on the graph's z)."""
        self._require_bound()
        source_rate = self._static_source_rate()
        for _ in range(int(nt)):
            incise_n1_step(self.network, self._order, self.K, self.m, self.U, dt,
                           source_rate=source_rate)
