"""
Gravel-bed transport closure (GRLP) for the diffusion solver.

A thin, external concrete closure: it subclasses ``fluvtree.closures.TransportClosure``
and supplies the gravel-bed transport law and hydraulic geometry. It depends on
the closure interface, never the reverse -- adding a transport law is a new
module like this one, not a core edit.

Naming: "gravel" is conventional shorthand, **not** a grain-size switch. The
discriminators are (1) **noncohesive**, grain-threshold banks (which set the
channel width) and (2) the **absence of form drag** (which sets the slope exponent
``p = 7/6``). Grain size ``D`` only resolves the diagnostic width and depth, not
the profile evolution. The common name is kept because the field uses it and it is
the usual gravel-bed river combination; this docstring defines what it means.

Provenance: GRLP v2.1.0 (grlp@366fb3e) -- Wong & Parker MPM, threshold width.
"""

from fluvtree.closures.base import TransportClosure


class GravelClosure(TransportClosure):
    """
    Gravel-bed closure (GRLP): **noncohesive** (grain-threshold) banks, **no form
    drag**. Wong & Parker MPM bedload under the threshold-width assumption -- the
    channel widens until the bed sits at ``(1+epsilon)`` times the grain-motion
    threshold.

    ``Q_s proportional to Q * S**(7/6)``. Grain size ``D`` is needed only to
    resolve width and depth, not to evolve the profile.
    """

    p = 7 / 6.0
    conductance_exponent = 1 / 6.0    # GRLP's literal (== p - 1, declared exactly)

    def __init__(self, D=None, lambda_p=0.35, rho_s=2650.0, rho=1000.0,
                 g=9.805, epsilon=0.2, tau_star_c=0.0495, phi=3.97):
        self.D = D
        self.lambda_p = lambda_p
        self.rho_s = rho_s
        self.rho = rho
        self.g = g
        self.epsilon = epsilon
        self.tau_star_c = tau_star_c
        self.phi = phi

    @property
    def R(self):
        """Submerged specific gravity ``(rho_s - rho) / rho``."""
        return (self.rho_s - self.rho) / self.rho

    @property
    def k_Qs(self):
        # GRLP Segment.bedload_lumped_constants: k_Qs = k_b * k_qs
        k_qs = (self.phi * self.R ** 0.5 * self.g ** 0.5
                * self.epsilon ** 1.5 * self.tau_star_c ** 1.5)
        k_b = (0.17 * self.g ** (-0.5) * self.R ** (-5 / 3.0)
               * (1 + self.epsilon) ** (-5 / 3.0) * self.tau_star_c ** (-5 / 3.0))
        return k_b * k_qs

    def _require_D(self):
        if self.D is None:
            raise ValueError("grain size D is required to resolve width/depth")

    def channel_width(self, Q, S):
        self._require_D()
        coeff = 0.17 / (self.g ** 0.5 * self.R ** (5 / 3.0)
                        * (1 + self.epsilon) ** (5 / 3.0)
                        * self.tau_star_c ** (5 / 3.0))
        return coeff * Q * S ** (7 / 6.0) / self.D ** 1.5

    def channel_depth(self, S):
        self._require_D()
        return self.R * (1 + self.epsilon) * self.tau_star_c * self.D / S
