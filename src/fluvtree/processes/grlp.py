"""
The GRLP long-profile process: door 1 (transport-limited) on the substrate.

A thin adapter that reuses GRLP's solver numerics unchanged. It reads the
canonical state off the FluvTree graph, drives ``grlp.Network``, and writes the
updated bed elevation ``z`` back to the graph. Under FluvTree's canonical-state
design the graph is the truth; the internal ``grlp.Network`` is a computational
cache, built **once** because the topology is fixed.

Reads ``x, z, Q, B`` (per-segment interior arrays on edges) and the boundary
fields ``S0`` (per channel-head node), ``x_bl``/``z_bl`` (outlet node), and
optional network-level ``Q_s_0`` (graph attribute). Writes ``z``.

State layout note: GRLP has no shared junction degree of freedom -- each segment
owns *all* its nodes, and a confluence is a special row in GRLP's matrix coupling
segment-owned nodes (``grlp@366fb3e``, ``Network.compute_Q_s``). So for this
process the graph edge carries the **full** GRLP segment array and the junction
nodes are topological markers, not shared ``z`` state. (A process that does use
per-node degrees of freedom -- e.g. one cross-section per node -- interprets the
nodes differently; that is a per-process choice, not a universal convention.)

Dependency direction (see the design note): FluvTree depends on published GRLP,
never the reverse.
"""

import numpy as np

from fluvtree.network import RiverNetwork


def build_grlp_network(upstream_segment_IDs, downstream_segment_IDs,
                       x, z, Q, B, S0, x_bl, z_bl, Q_s_0=None):
    """
    Build a :class:`RiverNetwork` carrying GRLP-consumable state.

    ``x, z, Q, B`` are per-segment interior arrays (edge fields). ``S0`` is a
    scalar boundary slope applied at every channel head, or one value per head in
    ascending head-segment-ID order (matching GRLP's channel-head ordering).
    ``x_bl``/``z_bl`` place base level at the outlet; ``Q_s_0`` is an optional
    network-level upstream sediment supply.
    """
    rn = RiverNetwork.from_segment_lists(upstream_segment_IDs,
                                         downstream_segment_IDs)
    segs = rn.segment_ids
    for i, s in enumerate(segs):
        # copy, so the graph owns independent arrays: guards against a caller
        # passing aliased arrays (e.g. ``[np.zeros(n)] * k``, which repeats one
        # object) -- otherwise the reaches would share and clobber one another.
        rn.set_segment_field(s, "x", np.array(x[i], dtype=float))
        rn.set_segment_field(s, "z", np.array(z[i], dtype=float))
        rn.set_segment_field(s, "Q", np.array(Q[i], dtype=float))
        rn.set_segment_field(s, "B", np.array(B[i], dtype=float))
    heads = sorted(rn.head_segments())
    try:
        iter(S0)
        S0_per_head = list(S0)
    except TypeError:
        S0_per_head = [S0] * len(heads)
    for head, _S0 in zip(heads, S0_per_head):
        rn.set_node_field(rn.edge_of(head)[0], "S0", _S0)
    outlet = rn.edge_of(rn.mouth_segments()[0])[1]
    rn.set_node_field(outlet, "x_bl", x_bl)
    rn.set_node_field(outlet, "z_bl", z_bl)
    if Q_s_0 is not None:
        rn.graph.graph["Q_s_0"] = Q_s_0
    return rn


class GRLP(object):
    """
    GRLP long-profile evolution as a FluvTree process.

    Parameters
    ----------
    network : RiverNetwork
        The canonical network. Must carry the fields listed in the module
        docstring (see :func:`build_grlp_network`).
    configure : callable, optional
        ``configure(grlp_network)`` called once after the internal
        ``grlp.Network`` is built, to set the integration scheme (e.g.
        ``set_time_integration`` / ``set_iteration_tolerance`` / ``set_niter``).
        Default: leave GRLP's shipped defaults (BDF2, iterate-to-convergence).
    """

    reads = ("x", "z", "Q", "B", "S0", "x_bl", "z_bl")
    writes = ("z",)

    def __init__(self, network, configure=None):
        import grlp

        self.network = network
        segs = list(network.segment_ids)
        self._segs = segs

        up = [sorted(network.upstream_segments(s)) for s in segs]
        down = []
        for s in segs:
            d = network.downstream_segment(s)
            down.append([] if d is None else [d])
        x = [network.get_segment_field(s, "x") for s in segs]
        z = [network.get_segment_field(s, "z").copy() for s in segs]
        Q = [network.get_segment_field(s, "Q") for s in segs]
        B = [network.get_segment_field(s, "B") for s in segs]

        heads = sorted(network.head_segments())
        S0 = [network.get_node_field(network.edge_of(h)[0], "S0") for h in heads]
        outlet = network.edge_of(network.mouth_segments()[0])[1]

        gnet = grlp.Network()
        gnet.initialize(
            x_bl=network.get_node_field(outlet, "x_bl"),
            z_bl=network.get_node_field(outlet, "z_bl"),
            S0=S0,
            Q_s_0=network.graph.graph.get("Q_s_0"),
            upstream_segment_IDs=up,
            downstream_segment_IDs=down,
            x=x, z=z, Q=Q, B=B,
        )
        if configure is not None:
            configure(gnet)
        gnet.get_z_lengths()
        self.grlp_network = gnet

    def _pull_z(self):
        """Refresh the internal solver's ``z`` from the canonical graph."""
        for i, s in enumerate(self._segs):
            self.grlp_network.segments[i].z = \
                self.network.get_segment_field(s, "z").copy()

    def _push_z(self):
        """Write the solver's updated ``z`` back to the canonical graph."""
        for i, s in enumerate(self._segs):
            self.network.set_segment_field(
                s, "z", self.grlp_network.segments[i].z.copy())

    def step(self, dt, nt=1):
        """
        Advance the long profile by ``nt`` steps of ``dt`` [s].

        Pulls the current ``z`` off the graph first (so an upstream process this
        step is seen), evolves, and writes ``z`` back. Only ``z`` is re-pulled;
        re-pulling ``B``/``Q``/``ssd`` is the hook for a lateral process (e.g.
        TerraPIN) and is deferred until one exists.
        """
        self._pull_z()
        self.grlp_network.evolve_threshold_width_river_network(nt=nt, dt=dt)
        self._push_z()
