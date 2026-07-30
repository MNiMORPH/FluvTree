"""
The process base class: the schedulable-operator protocol.

A *process* reads named fields off the :class:`~fluvtree.network.RiverNetwork`,
computes, and writes fields back. It owns no canonical state -- the graph does --
but it may cache network-derived setup (an index map, a solver, a sweep order),
which is why it *binds* to a network.

**Binding is deferred.** Construct a process with just its physics
(``ft.GravelBed(D=0.05)``) and let a :class:`~fluvtree.model.FluvTree` model (or a
:class:`~fluvtree.scheduler.Scheduler`) bind it to their network when it is added.
Passing a network at construction (``GravelBed(net, D=0.05)``) binds eagerly --
both forms work. Fixed topology makes any network-derived cache valid for the whole
run, so ``_on_bind`` is the right place to build it.
"""


class Process(object):
    """Base class for a schedulable process on a :class:`RiverNetwork`.

    Subclasses declare ``reads``/``writes``, build any network-derived cache in
    :meth:`_on_bind`, and implement :meth:`step`. Store physics parameters in
    ``__init__`` *before* calling ``super().__init__(network)``, so ``_on_bind`` can
    use them when a network is passed eagerly.
    """

    #: named graph fields this process consumes (for scheduling / validation)
    reads = ()
    #: named graph fields this process produces
    writes = ()

    def __init__(self, network=None):
        self.network = None
        if network is not None:
            self.bind(network)

    def bind(self, network):
        """Bind to a canonical network and build any network-derived cache.
        Returns ``self`` so callers can chain. Re-binding rebuilds the cache."""
        self.network = network
        self._on_bind()
        return self

    def _on_bind(self):
        """Hook: build network-derived state (a solver, a sweep order, an engine).
        Default: nothing -- a process that only reads/writes fields needs no cache."""
        pass

    def _require_bound(self):
        if self.network is None:
            raise RuntimeError(
                "process is not bound to a network; add it to a FluvTree or "
                "Scheduler (or pass a network at construction) before running")

    def step(self, dt, nt=1):
        """Advance the process ``nt`` steps of ``dt`` [s], in place on the graph."""
        raise NotImplementedError
