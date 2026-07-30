"""
n = 1 stream-power incision on a FluvTree network: the implicit outlet->upstream
sweep (Braun & Willett, 2013).

Implicit backward Euler for ``dz/dt = U - K Q**m S`` (``S`` = downstream slope,
discharge-based stream power). At a node ``i`` with downstream receiver ``r``:

    (z_i - z_i_old)/dt = U - K Q_i**m (z_i - z_r)/dx        (dx = x_r - x_i > 0)

which, using the *already-updated* receiver value ``z_r`` (the outlet->upstream
sweep), is closed-form and linear:

    z_i = (z_i_old + dt*U + c*z_r) / (1 + c),     c = dt*K*Q_i**m / dx

The receiver is the next node downstream within a reach; for a reach's last node it
is the first node of the downstream reach, or base level ``(x_bl, z_bl)`` at the
outlet. Because information flows upstream, sweeping reaches downstream-first and
nodes last-to-first means every receiver is solved before the node that needs it --
so it is exact in one pass, and unconditionally stable at any ``dt``.

State lives on the graph: ``x, z, Q`` per reach (edges), base level on the outlet
node. ``z`` arrays are mutated in place, so the graph is updated directly.
"""

import numpy as np


def sweep_order(network):
    """
    Reaches in downstream-first (outlet -> upstream) order -- the valid sweep
    order, since ``walk_upstream`` visits a reach before its upstream neighbours.
    Topology is fixed, so this is computed once and reused.
    """
    mouth = network.mouth_segments()[0]
    return network.walk_upstream(mouth)


def incise_n1_step(network, order, K, m, U, dt):
    """
    One implicit ``n = 1`` stream-power step over the network, in place on ``z``.

    ``K`` erodibility, ``m`` discharge exponent, ``U`` uplift rate (scalar), ``dt``
    time step. ``order`` from :func:`sweep_order`.
    """
    outlet_node = network.edge_of(network.mouth_segments()[0])[1]
    z_bl = network.get_node_field(outlet_node, "z_bl")
    x_bl = network.get_node_field(outlet_node, "x_bl")

    # start-of-step elevations (the time term), snapshot before any node moves
    z_old = {s: network.get_segment_field(s, "z").copy() for s in order}

    for s in order:
        x = network.get_segment_field(s, "x")
        z = network.get_segment_field(s, "z")          # mutated in place == graph
        Q = network.get_segment_field(s, "Q")
        L = len(z)

        # receiver of this reach's downstream-most node: the downstream reach's
        # first node (already updated), or base level at the outlet
        d = network.downstream_segment(s)
        if d is None:
            z_end, x_end = z_bl, x_bl
        else:
            z_end = network.get_segment_field(d, "z")[0]
            x_end = network.get_segment_field(d, "x")[0]

        for i in range(L - 1, -1, -1):
            if i == L - 1:
                z_r, x_r = z_end, x_end
            else:
                z_r, x_r = z[i + 1], x[i + 1]
            c = dt * K * Q[i] ** m / (x_r - x[i])
            z[i] = (z_old[s][i] + dt * U + c * z_r) / (1.0 + c)


def evolve_streampower_n1(network, K, m, U, dt, nt=1):
    """Advance ``nt`` implicit ``n = 1`` stream-power steps of ``dt`` in place."""
    order = sweep_order(network)
    for _ in range(int(nt)):
        incise_n1_step(network, order, K, m, U, dt)


class AdvectionSolver(object):
    """
    Power-law nonlinear-advection solver (stream-power incision) on a
    :class:`RiverNetwork`, the object-oriented sibling of :class:`DiffusionSolver`.

    Implements the ``n = 1`` linear rung of ``dz/dt = U - K Q**m S**n`` via the
    implicit outlet->upstream sweep (celerity ``~ S**(n-1)`` is slope-independent at
    ``n = 1``, so the sweep is exact in one pass). ``n != 1`` (a Newton step per
    node) is a later addition. ``K`` erodibility, ``m`` discharge exponent, ``U``
    uplift rate (scalar). Topology is fixed, so the sweep order is cached once.
    """

    def __init__(self, network, K, m=0.5, U=0.0):
        self.network = network
        self.K = K
        self.m = m
        self.U = U
        self._order = sweep_order(network)

    def evolve(self, nt, dt):
        """Advance ``nt`` implicit ``n = 1`` steps of ``dt`` [s], in place on the graph."""
        for _ in range(int(nt)):
            incise_n1_step(self.network, self._order, self.K, self.m, self.U, dt)
