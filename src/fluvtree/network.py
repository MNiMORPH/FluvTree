"""
The FluvTree network substrate: the topology, and where variables are held on it.

Extracted (not grafted) from GRLP's physics-free network machinery,
``grlp@366fb3e`` (v2.1.0); its history stays in GRLP. See
docs/fluvtree-engine-architecture.md for the architecture and the reasons the
history is not transplanted.

The network generates the *structure* of a river as a directed convergent graph
and **holds** the variables attached to it -- it does not **own** them. Other
modules add, remove, and modify those variables; the network is only where they
live, keyed by the structure:

- **edges are segments (reaches)** and hold along-reach *interior* arrays
  (``x, z, Q, ...``): the reach's interior nodes;
- **nodes are junctions** (channel heads, confluences, the outlet) and hold the
  scalar values at that point.

A reach's two endpoints live on the graph nodes it connects, held once and shared
by every process that touches that junction (the "reach-assembly convention":
``hstack(upstream node, edge interior, downstream node)``).

This module knows nothing about sediment transport and owns none of the variables
it holds. Traversal and attach/read/remove access only; a *process* (e.g. the GRLP
adapter) supplies the physics and owns what it places here.
"""

import networkx as nx
import numpy as np


class RiverNetwork(object):
    """
    A directed convergent graph of a river network: the structure that variables
    are attached to. It holds them, keyed by segment and node -- it does not own
    them; the modules that attach a variable own it.

    Node identities are opaque hashables (any networkx node key). Segment
    (reach) identities are integers carried on the edges as the ``seg`` attribute,
    so a segment is addressable independently of which two nodes it happens to
    connect -- topology is fixed, but segment ids never depend on node numbering.
    """

    def __init__(self, graph=None):
        self.graph = nx.DiGraph() if graph is None else graph
        self._rebuild_edge_index()

    # ------------------------------------------------------------------ #
    # Construction
    # ------------------------------------------------------------------ #

    @classmethod
    def from_segment_lists(cls, upstream_segment_IDs, downstream_segment_IDs):
        """
        Build the topology from per-segment upstream/downstream ID lists -- the
        way a GRLP network is specified (see grlp.Network.initialize).

        Segment ``s`` becomes one downstream-directed edge. Its upstream node is a
        channel head (``("source", s)``) if it has no upstream segments, else the
        confluence it flows out of (``("jcn", s)``); its downstream node is the
        outlet (``("outlet",)``) if it has no downstream segment, else the
        confluence named for the segment it flows into (``("jcn", d)``). A
        confluence node is named by the single segment leaving it, so every
        tributary entering it and that outgoing segment share the node.
        """
        n = len(upstream_segment_IDs)
        if len(downstream_segment_IDs) != n:
            raise ValueError("upstream/downstream lists must have equal length")
        G = nx.DiGraph()
        for s in range(n):
            up = upstream_segment_IDs[s]
            down = downstream_segment_IDs[s]
            u = ("jcn", s) if up else ("source", s)
            v = ("jcn", down[0]) if down else ("outlet",)
            G.add_edge(u, v, seg=s)
        return cls(G)

    @classmethod
    def from_arrays(cls, upstream_segment_IDs, downstream_segment_IDs,
                    x, z, Q, B, S0, x_bl, z_bl, Q_s_0=None):
        """
        Construct a network from explicit per-segment arrays and boundary values.

        The topology comes from the upstream/downstream ID lists (as in
        :meth:`from_segment_lists`); ``x, z, Q, B`` are per-segment interior arrays
        stamped onto the edges, and ``S0`` (a scalar applied at every channel head,
        or one value per head in ascending head-segment-ID order), ``x_bl``/``z_bl``
        (outlet base level), and optional network-level ``Q_s_0`` set the boundaries.

        This is the *explicit* constructor: the caller provides the values. Synthetic
        networks belong in ``fluvtree.generate``, DEM/measured loading in an ``io``
        layer -- see docs/DESIGN-structure-and-naming.md.
        """
        rn = cls.from_segment_lists(upstream_segment_IDs, downstream_segment_IDs)
        for i, s in enumerate(rn.segment_ids):
            # copy, so the graph holds independent arrays: guards against a caller
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

    # ------------------------------------------------------------------ #
    # Topology
    # ------------------------------------------------------------------ #

    def _rebuild_edge_index(self):
        """Cache seg id -> (u, v). Fixed topology: built once, valid for the run."""
        self._edge_of = {
            d["seg"]: (u, v) for u, v, d in self.graph.edges(data=True)
        }

    @property
    def segment_ids(self):
        """All segment (reach) ids, sorted."""
        return sorted(self._edge_of)

    def edge_of(self, seg):
        """The ``(upstream_node, downstream_node)`` pair carrying segment ``seg``."""
        return self._edge_of[seg]

    def head_segments(self):
        """Segments whose upstream node is a channel head (no inflow)."""
        return [s for s, (u, v) in self._edge_of.items()
                if self.graph.in_degree(u) == 0]

    def mouth_segments(self):
        """Segments whose downstream node is an outlet (no outflow)."""
        return [s for s, (u, v) in self._edge_of.items()
                if self.graph.out_degree(v) == 0]

    def upstream_segments(self, seg):
        """Segments immediately upstream of ``seg`` (its tributaries)."""
        u, _ = self._edge_of[seg]
        return [d["seg"] for _, _, d in self.graph.in_edges(u, data=True)]

    def downstream_segment(self, seg):
        """The single segment immediately downstream, or ``None`` at the outlet.

        Convergent tree: a segment has at most one downstream neighbour.
        """
        _, v = self._edge_of[seg]
        outs = [d["seg"] for _, _, d in self.graph.out_edges(v, data=True)]
        if len(outs) > 1:
            raise ValueError(
                "non-convergent topology: segment %r has %d downstream reaches"
                % (seg, len(outs)))
        return outs[0] if outs else None

    def walk_downstream(self, seg):
        """All segments from ``seg`` to the outlet, inclusive, in flow order."""
        out = [seg]
        nxt = self.downstream_segment(seg)
        while nxt is not None:
            out.append(nxt)
            nxt = self.downstream_segment(nxt)
        return out

    def walk_upstream(self, seg):
        """All segments upstream of (and including) ``seg``, depth-first."""
        out = []

        def _recurse(s):
            out.append(s)
            for up in self.upstream_segments(s):
                _recurse(up)

        _recurse(seg)
        return out

    # ------------------------------------------------------------------ #
    # Variable access (attach / read / remove -- the network holds, not owns)
    # ------------------------------------------------------------------ #

    def set_segment_field(self, seg, name, value):
        """Attach or modify the along-reach interior array ``name`` on segment
        ``seg``. The caller owns the variable; the network only holds it."""
        u, v = self._edge_of[seg]
        self.graph.edges[u, v][name] = value

    def get_segment_field(self, seg, name):
        """Read the along-reach interior array ``name`` from segment ``seg``."""
        u, v = self._edge_of[seg]
        return self.graph.edges[u, v][name]

    def remove_segment_field(self, seg, name):
        """Remove the along-reach interior array ``name`` from segment ``seg``."""
        u, v = self._edge_of[seg]
        del self.graph.edges[u, v][name]

    def set_node_field(self, node, name, value):
        """Attach or modify the scalar ``name`` on junction ``node``. The caller
        owns the variable; the network only holds it."""
        self.graph.nodes[node][name] = value

    def get_node_field(self, node, name):
        """Read a scalar ``name`` from junction ``node``."""
        return self.graph.nodes[node][name]

    def remove_node_field(self, node, name):
        """Remove the scalar ``name`` from junction ``node``."""
        del self.graph.nodes[node][name]

    def reach_profile(self, seg, name):
        """
        Assemble the full profile of field ``name`` over reach ``seg``:
        ``hstack(upstream node, edge interior, downstream node)`` -- the
        reach-assembly convention. Endpoint values come from the shared nodes.
        """
        u, v = self._edge_of[seg]
        return np.hstack((
            [self.graph.nodes[u][name]],
            self.graph.edges[u, v][name],
            [self.graph.nodes[v][name]],
        ))
