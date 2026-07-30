"""
Sand-bed transport closure (SRLP) for the diffusion solver.

A thin, external concrete closure: it subclasses ``fluvtree.closures.TransportClosure``
and supplies the sand-bed transport law and hydraulic geometry. It depends on
the closure interface, never the reverse -- adding a transport law is a new
module like this one, not a core edit.

Naming: "sand" is conventional shorthand, **not** a grain-size switch. The
discriminators are (1) **cohesive** (mud) banks -- the width is set by the critical
stress to erode them (``tau_crit_bank``), not by a grain-motion threshold -- and
(2) the effect of **form drag** (bedforms partitioning the total boundary stress),
which lowers the transport slope exponent to ``p = 5/6``. Grain size ``D`` only
resolves the diagnostic width and depth. The common name is kept because the field
uses it and it is the usual sand-bed river combination; this docstring defines what
it means.

Provenance: SRLP (MNiMORPH/SRLP) -- sand-bed transport, Manning friction.
"""

from fluvtree.closures.base import TransportClosure


class SandClosure(TransportClosure):
    """
    Sand-bed closure (SRLP): **cohesive** (mud) banks plus the effect of **form
    drag** on transport. The width is set by the critical stress to erode cohesive
    mud banks (``tau_crit_bank``); form drag lowers the transport slope exponent to
    5/6.

    ``Q_s proportional to Q * S**(5/6)``. Requires grain size ``D``, Manning's
    roughness ``n``, and the critical bank stress ``tau_crit_bank`` (mud banks).
    """

    p = 5 / 6.0
    conductance_exponent = -1 / 6.0   # SRLP's literal (== p - 1, declared exactly)

    def __init__(self, D, n, tau_crit_bank, lambda_p=0.35, rho_s=2650.0,
                 rho=1000.0, g=9.805, epsilon=0.2):
        self.D = D
        self.n = n
        self.tau_crit_bank = tau_crit_bank
        self.lambda_p = lambda_p
        self.rho_s = rho_s
        self.rho = rho
        self.g = g
        self.epsilon = epsilon

    @property
    def R(self):
        """Submerged specific gravity ``(rho_s - rho) / rho``."""
        return (self.rho_s - self.rho) / self.rho

    @property
    def k_Qs(self):
        # SRLP Segment.sediment_lumped_constants
        return ((0.05 / self.n) * 1.0 / (self.g ** 0.5 * self.R ** 2)
                * ((1 + self.epsilon) * self.tau_crit_bank
                   / (self.rho * self.g)) ** (7 / 6.0) * (1.0 / self.D))

    def channel_width(self, Q, S):
        # SRLP compute_channel_width (paper draft Eq. 11)
        return (self.n * (self.rho * self.g) ** (5 / 3.0)
                * ((1 + self.epsilon) * self.tau_crit_bank) ** (-5 / 3.0)
                * Q * S ** (7 / 6.0))

    def channel_depth(self, S):
        # SRLP compute_flow_depth (paper draft Eq. 7)
        return (1 + self.epsilon) * self.tau_crit_bank / (self.rho * self.g * S)
