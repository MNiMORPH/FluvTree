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

from fluvtree.common.rates import broadcast_rate


def uplift_rate(network, rate):
    """
    Per-segment vertical source rate [m/s] from an uplift/subsidence specification.

    ``rate`` is positive for uplift, negative for subsidence, and may be a scalar
    (uniform), one value/array per segment, or a callable ``rate(x)`` for a spatial
    pattern (see :func:`fluvtree.common.rates.broadcast_rate`). Returns one array per
    segment (segment-id order) -- pass to the diffusion solver's ``source_rate`` or a
    process's ``uplift``.
    """
    return broadcast_rate(network, rate)
