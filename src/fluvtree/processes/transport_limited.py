"""
A transport-limited long-profile rule for door 2 -- GRLP's physics, self-contained.

This is a demonstration that the explicit-rule hook (:class:`Rule`) can
carry real geomorphic physics and reproduce the GRLP test set, with **no GRLP
dependency**: the flux law and boundary conditions are re-expressed as a
finite-volume Exner update on the FluvTree substrate. Because a transport-limited
steady state is solver-independent (constant sediment flux => constant slope), this
explicit rule relaxes to the *same* steady state GRLP's implicit solver reaches.

Sediment flux on a face: ``Q_s = k_Qs * I * Q_face * S**(7/6)`` (Meyer-Peter-Muller
bedload, threshold-width channel), the same closure GRLP uses. Exner mass balance
at a node: ``dz/dt = -(1 / ((1 - lambda_p) * B)) * dQ_s/dx``. Boundaries: a
prescribed supply slope ``S0`` at each channel head, base level ``(x_bl, z_bl)`` at
the outlet, and at a confluence the shared downstream node receives the summed flux
of its tributaries (flux conservation => downstream ``Q_s`` = sum of head supplies,
GRLP's during-flood conservation).

The default ``k_Qs`` is computed from the same physical constants as GRLP, so the
reproduction is numerical, not merely qualitative.

Explicit, so CFL-bound (the door-2 tradeoff): fine on the small test domains; not a
substitute for door 1's implicit large-``dt`` stability.
"""

import numpy as np


def default_k_Qs(lambda_p=0.35, rho_s=2650.0, rho=1000.0, g=9.805,
                 epsilon=0.2, tau_star_c=0.0495, phi=3.97):
    """
    The lumped bedload coefficient ``k_Qs`` from GRLP's default constants
    (``Segment.basic_constants`` + ``bedload_lumped_constants``), so a rule using
    it reproduces GRLP numerically.
    """
    k_qs = phi * ((rho_s - rho) / rho) ** 0.5 * g ** 0.5 \
        * epsilon ** 1.5 * tau_star_c ** 1.5
    k_b = 0.17 * g ** (-0.5) * ((rho_s - rho) / rho) ** (-5 / 3.0) \
        * (1 + epsilon) ** (-5 / 3.0) * tau_star_c ** (-5 / 3.0)
    return k_b * k_qs


class TransportLimitedRate(object):
    """
    A finite-volume transport-limited Exner rate, callable as ``rate(network, dt)``
    for :class:`~fluvtree.processes.Rule`.

    Reads per-segment ``x, z, Q, B`` (edges) and the boundary fields ``S0`` (head
    nodes) and ``x_bl``/``z_bl`` (outlet node) off the graph -- the layout that
    :func:`~fluvtree.processes.build_grlp_network` stamps. Holds the transport
    constants as parameters.

    Parameters
    ----------
    k_Qs : float, optional
        Lumped bedload coefficient (default: :func:`default_k_Qs`, GRLP's value).
    intermittency : float, optional
        Flood intermittency ``I`` (default 1.0).
    lambda_p : float, optional
        Sediment porosity (default 0.35, GRLP's value).
    """

    def __init__(self, k_Qs=None, intermittency=1.0, lambda_p=0.35):
        self.k_Qs = default_k_Qs() if k_Qs is None else k_Qs
        self.intermittency = intermittency
        self.lambda_p = lambda_p

    def _face_Qs(self, Q_face, S):
        """Signed sediment flux through a face of slope ``S`` (downhill positive)."""
        return np.sign(S) * self.k_Qs * self.intermittency \
            * Q_face * np.abs(S) ** (7 / 6.0)

    def __call__(self, network, dt):
        segs = network.segment_ids
        z = {s: network.get_segment_field(s, "z") for s in segs}
        x = {s: network.get_segment_field(s, "x") for s in segs}
        Q = {s: network.get_segment_field(s, "Q") for s in segs}
        B = {s: network.get_segment_field(s, "B") for s in segs}

        rates = {}
        for s in segs:
            xs, zs, Qs_, Bs = x[s], z[s], Q[s], B[s]
            L = len(zs)
            flux = np.zeros(L + 1)  # face fluxes: flux[i] enters node i, flux[i+1] leaves

            # --- upstream face of node 0: head supply, or summed tributaries ---
            ups = network.upstream_segments(s)
            if not ups:
                head_node = network.edge_of(s)[0]
                S0 = network.get_node_field(head_node, "S0")
                flux[0] = self.k_Qs * self.intermittency * Qs_[0] * S0 ** (7 / 6.0)
            else:
                total = 0.0
                for u in ups:
                    dist = xs[0] - x[u][-1]
                    total += self._face_Qs(Q[u][-1], (z[u][-1] - zs[0]) / dist)
                flux[0] = total

            # --- interior faces ---
            S_int = (zs[:-1] - zs[1:]) / (xs[1:] - xs[:-1])
            Q_int = 0.5 * (Qs_[:-1] + Qs_[1:])
            flux[1:L] = self._face_Qs(Q_int, S_int)

            # --- downstream face of the last node: base level, or the confluence ---
            d = network.downstream_segment(s)
            if d is None:
                outlet_node = network.edge_of(s)[1]
                x_bl = network.get_node_field(outlet_node, "x_bl")
                z_bl = network.get_node_field(outlet_node, "z_bl")
                flux[L] = self._face_Qs(Qs_[-1], (zs[-1] - z_bl) / (x_bl - xs[-1]))
            else:
                dist = x[d][0] - xs[-1]
                flux[L] = self._face_Qs(Qs_[-1], (zs[-1] - z[d][0]) / dist)

            # --- node control-volume lengths (half-cells; ends use the end cell) ---
            Lv = np.empty(L)
            Lv[1:-1] = (xs[2:] - xs[:-2]) / 2.0
            Lv[0] = xs[1] - xs[0]
            Lv[-1] = xs[-1] - xs[-2]

            rates[s] = (flux[:-1] - flux[1:]) / ((1 - self.lambda_p) * Bs * Lv)

        return rates
