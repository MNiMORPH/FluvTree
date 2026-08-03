"""
FluvTree analysis: derive quantities and metrics from a network's state.

Pure numpy, never matplotlib -- analysis returns *data* (arrays, scalars); the
separate ``fluvtree.plot`` layer consumes it and draws (one-way: plot -> analysis).
Keeping analysis matplotlib-free is what lets ``import fluvtree`` stay light.

Two axes:

- ``analysis.network`` -- *across the network* (spatial): sediment discharge
  ``Q_s``, the slope-area relationship, long-profile assembly. Planned home for the
  sediment budget and network morphometry (see docs/GRLP-parity-and-gaps.md).
- ``analysis.dynamics`` -- *in dynamics* (temporal): how the state evolves in time
  (planned; e.g. dz/dt, response/equilibration times, transfer functions).
"""

from fluvtree.analysis.network import compute_Q_s, slope_area

__all__ = ["compute_Q_s", "slope_area"]
