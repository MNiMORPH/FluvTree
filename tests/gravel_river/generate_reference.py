"""
Generate the gravel-river golden reference from stable GRLP.

Fixed, re-runnable code (run it to regenerate ``reference.npz``): builds each
topology in :mod:`topologies` as a ``grlp.Network``, evolves it to steady state
with the shipped BDF2 + iterate-to-convergence solver, and saves the steady-state
bed elevation and during-flood sediment discharge per segment. These GRLP steady
states are the authority the door-2 rule must reproduce.

    python tests/gravel_river/generate_reference.py

Provenance: grlp v2.1.0 (grlp@366fb3e).
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import grlp
from topologies import SPECS, widths, B, S0, Z_BL

# Large steps drive these small domains to machine-precision steady state.
NT = 1000
DT = 3.0e10


def _steady_grlp(spec):
    n = len(spec["x"])
    n_heads = sum(1 for i in range(n) if len(spec["up"][i]) == 0)
    net = grlp.Network()
    net.initialize(
        x_bl=spec["x_bl"], z_bl=Z_BL, S0=[S0] * n_heads, Q_s_0=None,
        upstream_segment_IDs=spec["up"], downstream_segment_IDs=spec["down"],
        x=[xi.copy() for xi in spec["x"]],
        z=[np.zeros(len(xi)) for xi in spec["x"]],
        Q=[qi.copy() for qi in spec["Q"]],
        B=widths(spec),
    )
    net.set_time_integration(2)
    net.set_iteration_tolerance(1.0e-4)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=NT, dt=DT)
    return net.list_of_LongProfile_objects


def _during_flood_Qs(lp):
    """Sediment discharge on interior faces (GRLP's network-correctness metric)."""
    S = np.abs(np.diff(lp.z) / np.diff(lp.x))
    Q_mid = 0.5 * (lp.Q[:-1] + lp.Q[1:])
    return lp.k_Qs * Q_mid * S ** (7 / 6.0)


def main():
    out = {}
    for name, spec in SPECS.items():
        segs = _steady_grlp(spec)
        out["%s__z" % name] = np.hstack([lp.z for lp in segs])
        out["%s__lengths" % name] = np.array([len(lp.z) for lp in segs])
        out["%s__qs" % name] = np.hstack([_during_flood_Qs(lp) for lp in segs])
    path = os.path.join(os.path.dirname(__file__), "reference.npz")
    np.savez(path, **out)
    print("wrote %s (%d topologies)" % (path, len(SPECS)))


if __name__ == "__main__":
    main()
