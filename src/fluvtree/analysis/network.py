"""
Across-the-network (spatial) analysis: quantities derived along the river network.

Pure numpy; returns data, not figures. The sediment discharge ``Q_s`` (the same
topology walk GRLP's ``Network.compute_Q_s`` uses) and the slope-area relationship
live here; the sediment budget and network morphometry are the planned additions
(see docs/GRLP-parity-and-gaps.md). ``fluvtree.plot`` draws what these return.
"""

import numpy as np


def compute_Q_s(network, closure, z, intermittency=1.0, sinuosity=1.0):
    """
    Sediment discharge ``Q_s = k_Qs * I * Q * S**p`` at each node of the network.

    Slope ``S`` is a central difference on the ghost-padded profile (the boundary
    slope ``S0`` at a channel head, base level at the outlet, the neighbour's
    endpoint across an interior junction); at a confluence head the tributaries
    give one profile each and ``Q_s`` is their average -- the same walk GRLP's
    ``Network.compute_Q_s`` uses.

    Parameters
    ----------
    network : RiverNetwork
    closure : TransportClosure
        Supplies ``k_Qs`` and ``p``.
    z : list of ndarray
        Current per-segment elevation, in ``network.segment_ids`` order.
    intermittency, sinuosity : float, optional

    Returns
    -------
    list of ndarray
        Per-segment ``Q_s`` (signed positive downhill), in segment-id order.
    """
    segs = list(network.segment_ids)
    zmap = {s: np.asarray(z[i], float) for i, s in enumerate(segs)}
    k_Qs, p = closure.k_Qs, closure.p
    out = []
    for s in segs:
        zs = zmap[s]
        xs = np.asarray(network.get_segment_field(s, "x"), float)
        Qseg = np.asarray(network.get_segment_field(s, "Q"), float)

        # downstream ghost: base level at the outlet, else the downstream reach's
        # first node
        d = network.downstream_segment(s)
        if d is None:
            outlet = network.edge_of(s)[1]
            z_down = network.get_node_field(outlet, "z_bl")
            try:
                x_down = network.get_node_field(outlet, "x_bl")
            except KeyError:
                x_down = 2 * xs[-1] - xs[-2]
        else:
            z_down = zmap[d][0]
            x_down = network.get_segment_field(d, "x")[0]

        # upstream ghost(s): S0 at a channel head, else each tributary's last node
        z_up, x_up = [], []
        ups = network.upstream_segments(s)
        if not ups:
            head = network.edge_of(s)[0]
            S0 = network.get_node_field(head, "S0")
            _xg = 2 * xs[0] - xs[1]
            x_up.append(_xg)
            z_up.append(zs[0] + S0 * (xs[0] - _xg))
        else:
            for u in ups:
                z_up.append(zmap[u][-1])
                x_up.append(float(network.get_segment_field(u, "x")[-1]))

        Q_s = []
        for _zu, _xu in zip(z_up, x_up):
            _z = np.hstack((_zu, zs, z_down))
            _x = np.hstack((_xu, xs, x_down))
            _dx = _x[2:] - _x[:-2]
            S = np.abs((_z[2:] - _z[:-2]) / _dx) / sinuosity
            Q_s.append(-np.sign(_z[2:] - _z[:-2])
                       * k_Qs * intermittency * Qseg * S ** p)
        out.append(np.mean(Q_s, axis=0))
    return out


def slope_area(network, against="Q"):
    """
    The slope-area relationship, as data (``fluvtree.plot.slope_area`` draws it).

    Bed slope ``|dz/dx|`` at each interior interval of every reach, against an
    area-like abscissa (``against``, a per-segment field: default discharge ``"Q"``,
    ``"A"`` for drainage area) at the interval midpoint. Returns ``(abscissa,
    slope)`` as flat arrays over all reaches, positives only, so a log-log plot or
    power-law fit can consume them directly.
    """
    slopes, absc = [], []
    for s in network.segment_ids:
        x = np.asarray(network.get_segment_field(s, "x"), float)
        z = np.asarray(network.get_segment_field(s, "z"), float)
        a = np.asarray(network.get_segment_field(s, against), float)
        slopes.append(np.abs(np.diff(z) / np.diff(x)))
        absc.append(0.5 * (a[:-1] + a[1:]))
    S = np.concatenate(slopes)
    A = np.concatenate(absc)
    keep = (S > 0) & (A > 0)                 # log-log needs positives
    return A[keep], S[keep]
