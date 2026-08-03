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

# Q_s is a general across-the-network diagnostic, not attrition-specific; it lives
# in fluvtree.analysis and is imported here (the fining sink is built from it).
from fluvtree.analysis.network import compute_Q_s


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
