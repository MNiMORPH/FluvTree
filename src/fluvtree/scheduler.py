"""
The scheduler: run an ordered ruleset over the canonical network, each time step.

This is the generalization of the hand-written GRLP<->TerraPIN coupling loop (see
docs/terrapin-grlp-coupling.md): the *sequence of processes is data, not code*. A
``Scheduler`` holds the canonical :class:`~fluvtree.network.RiverNetwork` and an
ordered list of **processes** (the ruleset). Each step, it runs every process's
``step(dt)`` in order -- each reads and writes the shared graph, so process *k*
sees what processes ``0..k-1`` just did -- then advances time.

Order is the modeller's to set (e.g. a lateral process before the vertical one).
Processes declare ``reads``/``writes`` (see the process interface); validating an
ordering against those is a deferred hook, not done here yet.

Canonical time lives on the graph (``G.graph["t"]``), so "what is" includes when.
"""


class Scheduler(object):
    """
    Advance a river network by running an ordered ruleset of processes.

    Parameters
    ----------
    network : RiverNetwork
        The canonical state every process reads and writes.
    processes : iterable, optional
        The ordered ruleset. Each process must implement ``step(dt)`` and, if it
        binds a network, bind *this* one. Appendable later with :meth:`add`.
    """

    def __init__(self, network, processes=None):
        self.network = network
        self.processes = []
        self.t = 0.0
        network.graph.graph.setdefault("t", self.t)
        if processes is not None:
            for p in processes:
                self.add(p)

    def add(self, process):
        """Append a process to the ruleset (runs after those already added)."""
        bound = getattr(process, "network", None)
        if bound is not None and bound is not self.network:
            raise ValueError(
                "process is bound to a different RiverNetwork than the scheduler")
        self.processes.append(process)
        return process

    def step(self, dt):
        """Run every process once, in order, then advance canonical time by ``dt``."""
        for process in self.processes:
            process.step(dt)
        self.t += dt
        self.network.graph.graph["t"] = self.t
        return self.t

    def run(self, nt, dt):
        """Advance ``nt`` steps of ``dt`` [s]; return the new canonical time."""
        for _ in range(int(nt)):
            self.step(dt)
        return self.t
