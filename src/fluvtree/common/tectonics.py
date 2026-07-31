"""
Tectonic forcing of the bed (optional).

Vertical motion for starters -- **uplift and subsidence** -- as a source in the bed
mass balance ``dz/dt = transport + U`` (GRLP's ``set_uplift_rate``). Positive is
uplift (the bed rises relative to base level, equivalently base-level fall);
negative is subsidence. It is a *source*, not a separate operator-split step, so it
couples into the implicit solve exactly as in GRLP.

Off by default (a common, opt-in module). Non-vertical motions -- horizontal
advection of the profile, fault offset, tilting -- are a later extension; the
interface below (scalar / per-segment / spatial rate) is the vertical case, which
those would sit alongside.
"""

import numpy as np


def uplift_rate(network, rate):
    """
    Per-segment vertical source rate [m/s] from an uplift/subsidence specification.

    Parameters
    ----------
    network : RiverNetwork
    rate : float, sequence, or callable
        The uplift rate (positive) or subsidence (negative). One of:
        a scalar (uniform); one value/array per segment (in ``segment_ids`` order);
        or a callable ``rate(x) -> array`` for a spatial pattern along each reach.

    Returns
    -------
    list of ndarray
        The source rate on each segment's nodes, in segment-id order -- pass to the
        diffusion solver's ``source_rate`` (or a process's ``uplift``).
    """
    out = []
    for i, s in enumerate(network.segment_ids):
        x = np.asarray(network.get_segment_field(s, "x"), float)
        if callable(rate):
            out.append(np.broadcast_to(np.asarray(rate(x), float), x.shape).copy())
        elif np.isscalar(rate):
            out.append(np.full(x.shape, float(rate)))
        else:
            out.append(np.broadcast_to(np.asarray(rate[i], float), x.shape).copy())
    return out
