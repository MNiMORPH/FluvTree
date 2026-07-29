"""
Shared gravel-river test topologies (the reference "test set").

Each spec fully specifies a fixed-topology network: per-segment node positions
``x`` and discharge ``Q``, the upstream/downstream segment-ID lists, and the
boundary conditions (a head supply slope ``S0``, base level ``(x_bl, z_bl)``, and a
uniform valley width ``B``). Both the golden-reference generator and the door-2
reproduction test import this module, so the physics is pinned in one place.

Steady states here are solver-independent (constant sediment flux => constant
slope), so they are the natural golden benchmark: GRLP's implicit solver and the
door-2 finite-volume rule must reach the same profile.
"""

import numpy as np

_D = 2000.0

B = 100.0
S0 = 0.015
Z_BL = 0.0


def _x(a, b):
    """Node positions from ``a*_D`` to ``(b-1)*_D`` inclusive, spacing ``_D``."""
    return _D * np.arange(a, b, dtype=float)


SPECS = {
    "single_segment": dict(
        up=[[]], down=[[]],
        x=[_x(1, 11)],
        Q=[10.0 * np.ones(10)],
        x_bl=_D * 11,
    ),
    "symmetric_confluence": dict(
        up=[[], [], [0, 1]], down=[[2], [2], []],
        x=[_x(1, 5), _x(1, 5), _x(5, 9)],
        Q=[5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)],
        x_bl=_D * 9,
    ),
    "asymmetric_Q_confluence": dict(
        up=[[], [], [0, 1]], down=[[2], [2], []],
        x=[_x(1, 5), _x(1, 5), _x(5, 9)],
        Q=[5 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
        x_bl=_D * 9,
    ),
    "multi_level": dict(
        up=[[], [], [], [], [0, 1], [4, 2], [5, 3]],
        down=[[4], [4], [5], [6], [5], [6], []],
        x=[_x(1, 5), _x(1, 5), _x(1, 5), _x(1, 5),
           _x(5, 9), _x(9, 13), _x(13, 17)],
        Q=[5 * np.ones(4)] * 4
          + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
        x_bl=_D * 17,
    ),
}


def widths(spec):
    """Per-segment width arrays (uniform ``B``) for a spec."""
    return [B * np.ones(len(xi)) for xi in spec["x"]]
