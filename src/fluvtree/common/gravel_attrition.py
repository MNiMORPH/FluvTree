"""
Sternberg downstream fining of the transported gravel load (optional).

As gravel is carried it abrades -- grains lose mass and fine downstream (Sternberg
1875). The abraded material leaves the bedload as wash load, so it is a *sink* in
the sediment mass balance. This is GRLP's ``update_gravel_loss`` /
``gravel_fractional_loss_per_km``, lifted into a shared, opt-in module: it depends
only on the gravel load and its flux, so it applies to any river carrying gravel --
gravel-bedded or bedrock-floored. It computes **nothing about the bed type**;
bedrock abrasion (tools wearing the rock) is a separate concern.

The fining sink is ``-(coeff/1000) * Q_s / ((1 - lambda_p) * B)`` in bed-lowering
rate [m/s], with ``coeff`` the fractional load loss per km and ``Q_s`` the local
sediment discharge. Because ``Q_s`` depends on the evolving bed, the sink is
relinearized every Picard iterate (fed to the diffusion solver as a dynamic
source), exactly as in GRLP.
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


def fining_rate(network, closure, z, coeff, intermittency=1.0, sinuosity=1.0):
    """
    Sternberg fining sink [m/s] per segment: ``-(coeff/1000) Q_s / ((1-lambda_p) B)``.

    ``coeff`` is the gravel fractional load loss per km. Negative (a sink: gravel
    leaves the bedload as it abrades to wash load). Sign follows GRLP exactly.
    """
    Q_s = compute_Q_s(network, closure, z, intermittency, sinuosity)
    lam = closure.lambda_p
    out = []
    for i, s in enumerate(network.segment_ids):
        B = np.asarray(network.get_segment_field(s, "B"), float)
        out.append(-coeff / 1000.0 * Q_s[i] / ((1.0 - lam) * B))
    return out


def dynamic_source(network, closure, coeff, intermittency=1.0, sinuosity=1.0):
    """
    Build the per-Picard-iterate dynamic-source callback the diffusion solver takes.

    Returns ``f(z) -> list of ndarray`` (the fining sink [m/s] for the current
    iterate ``z``), so the sink relinearizes with the evolving bed exactly as
    GRLP's ``update_gravel_loss`` does inside its Picard loop.
    """
    def _source(z):
        return fining_rate(network, closure, z, coeff, intermittency, sinuosity)
    return _source
