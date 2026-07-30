"""
The FluvTree model: holds a river and the ordered ruleset that evolves it.

*FluvTree holds what is; its methods modify what is.* A :class:`FluvTree` wraps a
:class:`~fluvtree.network.RiverNetwork` -- the structure and its variables -- and
owns the **processes** that modify those variables, together with the scheduler
that runs them and the clock. It is the single object a user instantiates, attaches
physics to, and runs::

    import fluvtree as ft
    net   = ft.RiverNetwork.from_segment_lists(up, down)   # the structure + state
    model = ft.FluvTree(net)                               # the model
    model.add(ft.GravelBed(D=0.05))                        # attach physics
    model.run(until=3.0e12, dt=3.0e10)                     # run it

The :class:`~fluvtree.scheduler.Scheduler` is the machinery inside :meth:`run`; use
it directly for advanced ordered-ruleset control. Processes attached here are bound
to this model's network automatically (see :class:`~fluvtree.processes.base.Process`).
"""

from fluvtree.scheduler import Scheduler


class FluvTree(object):
    """
    A river-network model: a :class:`RiverNetwork` plus its ordered ruleset.

    Parameters
    ----------
    network : RiverNetwork
        The canonical structure and state the model evolves.
    processes : iterable, optional
        An initial ordered ruleset; each is attached (and bound) as by :meth:`add`.
        More can be added later.
    """

    def __init__(self, network, processes=None):
        self.network = network
        self._scheduler = Scheduler(network)
        if processes is not None:
            for p in processes:
                self.add(p)

    # ------------------------------------------------------------------ #
    # State
    # ------------------------------------------------------------------ #

    @property
    def processes(self):
        """The ordered ruleset (list of bound processes), in run order."""
        return self._scheduler.processes

    @property
    def t(self):
        """Canonical model time [s] -- it lives on the graph, so it *is* state."""
        return float(self.network.graph.graph.get("t", 0.0))

    # ------------------------------------------------------------------ #
    # Building the ruleset
    # ------------------------------------------------------------------ #

    def add(self, process):
        """
        Attach a process to the ruleset, after those already added.

        If the process is unbound (constructed with just its physics, e.g.
        ``ft.GravelBed(D=0.05)``), it is bound to this model's network here. Returns
        ``self`` so attachments chain: ``model.add(a).add(b)``.
        """
        self._scheduler.add(process)
        return self

    # ------------------------------------------------------------------ #
    # Running
    # ------------------------------------------------------------------ #

    def run(self, dt, nt=None, until=None):
        """
        Advance the model, running the whole ruleset each step.

        Give exactly one of:

        - ``nt`` -- take ``nt`` steps of size ``dt``;
        - ``until`` -- step by ``dt`` until the canonical time reaches ``until``
          (the final step is shortened so time lands exactly on the target).

        Returns the new canonical time :attr:`t`.
        """
        if (nt is None) == (until is None):
            raise ValueError("give exactly one of nt or until")
        if nt is not None:
            self._scheduler.run(int(nt), dt)
        else:
            target = float(until)
            eps = 1.0e-9 * max(abs(target), 1.0)
            while self.t < target - eps:
                self._scheduler.step(min(dt, target - self.t))
        return self.t

    # ------------------------------------------------------------------ #
    # Views
    # ------------------------------------------------------------------ #

    def plot(self, ax=None, **kwargs):
        """
        Plot the current long profile (elevation vs downstream distance).

        Delegates to :func:`fluvtree.plot.long_profile`, imported lazily so bare
        ``import fluvtree`` pulls in no matplotlib. Returns the matplotlib ``Axes``.
        """
        from fluvtree.plot import long_profile
        return long_profile(self.network, ax=ax, **kwargs)
