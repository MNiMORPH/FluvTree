"""
Door 2: an arbitrary explicit rule as a process.

``Rule`` forward-integrates a user-supplied rate over the canonical graph.
This is the general escape hatch: any evolution law expressible as a rate of change
of a per-segment field can be run as a FluvTree process, with no implicit machinery
and no ties to GRLP.

The rate callable has the signature ``rate(network, dt) -> {seg_id: dfield_dt}``,
returning one interior array per segment it wants to change (segments it omits are
left unchanged). ``step(dt)`` advances the field by explicit (forward-Euler)
integration: ``field <- field + dt * rate``.

The tradeoff, as the design note flags: generality is bought with the stability
limit of explicit stepping (the rate's own CFL condition). Door 1 (the implicit
GRLP process) keeps large-``dt`` stability; door 2 does not.
"""

import numpy as np

from fluvtree.processes.base import Process


class Rule(Process):
    """
    Evolve a per-segment field by explicitly integrating a user rate.

    Parameters
    ----------
    network : RiverNetwork, optional
        The canonical network the rule reads and writes. Omit to bind later via a
        :class:`FluvTree` model or :class:`Scheduler`.
    rate : callable
        ``rate(network, dt) -> {seg_id: dfield_dt_array}``. Each returned array is
        the time derivative of ``field`` on that segment's interior nodes.
    field : str, optional
        The per-segment field to evolve (default ``"z"``).
    """

    def __init__(self, network=None, rate=None, field="z"):
        self.rate = rate
        self.field = field
        self.writes = (field,)
        super().__init__(network)

    def step(self, dt, nt=1):
        """Advance ``field`` by one forward-Euler step of ``dt`` [s]."""
        self._require_bound()
        rates = self.rate(self.network, dt)
        for seg, dfield_dt in rates.items():
            value = self.network.get_segment_field(seg, self.field)
            self.network.set_segment_field(
                seg, self.field, value + dt * np.asarray(dfield_dt, dtype=float))
