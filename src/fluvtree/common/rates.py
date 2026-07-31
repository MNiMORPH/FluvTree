"""
Shared helper: broadcast a rate specification onto the network's nodes.

Several optional forcings -- tectonic uplift, a distributed sediment source/sink --
are just a per-node rate [m/s] added to the bed mass balance. They accept the same
flexible specification (uniform scalar, per-segment, or a spatial function of ``x``),
normalized here to one array per segment so the solver receives a plain
``source_rate``.
"""

import numpy as np


def broadcast_rate(network, spec):
    """
    Normalize a rate ``spec`` to a per-segment list of node arrays [in ``x`` order].

    ``spec`` is one of: a scalar (uniform on every node); one value/array per segment
    (in ``segment_ids`` order); or a callable ``spec(x) -> array`` giving a spatial
    pattern along each reach.
    """
    out = []
    for i, s in enumerate(network.segment_ids):
        x = np.asarray(network.get_segment_field(s, "x"), float)
        if callable(spec):
            out.append(np.broadcast_to(np.asarray(spec(x), float), x.shape).copy())
        elif np.isscalar(spec):
            out.append(np.full(x.shape, float(spec)))
        else:
            out.append(np.broadcast_to(np.asarray(spec[i], float), x.shape).copy())
    return out
