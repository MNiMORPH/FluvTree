"""
The transport-closure interface: the template solver's plug socket.

An alluvial long-profile model is the shared :class:`DiffusionSolver` plus
a small, swappable *closure* -- the sediment-transport law and the channel
hydraulic geometry. Gravel-bed (GRLP) and sand-bed (SRLP) rivers differ *only* in
this closure; the solver (network assembly, implicit time stepping, Exner mass
balance) is identical. This module is that seam: the abstract ``TransportClosure``
contract only. The concrete gravel and sand closures live *outside* the diffusion solver,
in the thin ``gravel`` and ``sand`` modules that depend on this interface -- adding
a transport law is a new module, never a core edit.

Every closure shares one flux form:

    Q_s = k_Qs * I * Q * S**p        (signed by the slope)

and supplies three things: the lumped coefficient ``k_Qs``, the slope exponent
``p`` (gravel 7/6, sand 5/6), and the width/depth relations. The exponent also
threads into the *implicit* solver, whose sediment-flux conductance goes as
``S**(p - 1)`` (gravel ``S**(1/6)``, sand ``S**(-1/6)``); the solver reads ``p``
and this ``conductance_exponent`` off the closure rather than hard-coding them.
"""

import numpy as np


class TransportClosure(object):
    """
    Base class for an alluvial-river transport + hydraulic-geometry closure.

    Subclasses set the slope exponent ``p`` and provide ``k_Qs`` and the
    ``channel_width`` / ``channel_depth`` relations. The solver consumes ``k_Qs``,
    ``p`` (and ``conductance_exponent`` = ``p - 1``), and ``lambda_p``.
    """

    #: sediment-flux slope exponent (``Q_s`` proportional to ``S**p``)
    p = None

    #: sediment porosity, used by the solver's Exner balance ``(1 - lambda_p)``
    lambda_p = 0.35

    @property
    def k_Qs(self):
        """Lumped sediment-transport coefficient (model-specific)."""
        raise NotImplementedError

    #: Slope exponent of the implicit sediment-flux conductance (``= p - 1``).
    #: Declared as an explicit literal by each closure (gravel 1/6, sand -1/6) so
    #: the solver matches the source model *bit-for-bit* -- ``p - 1`` computed in
    #: float differs from the literal in the last bits. The base value is the
    #: computed fallback for custom closures.
    @property
    def conductance_exponent(self):
        """Slope exponent of the implicit sediment-flux conductance (``p - 1``)."""
        return self.p - 1.0

    def sediment_discharge(self, Q, S, intermittency=1.0):
        """
        During-flood sediment discharge ``Q_s = k_Qs * I * Q * S**p``, signed by
        the slope (positive downhill). Shared functional form; ``k_Qs`` and ``p``
        carry the model.
        """
        S = np.asarray(S, dtype=float)
        return (np.sign(S) * self.k_Qs * intermittency
                * np.asarray(Q, dtype=float) * np.abs(S) ** self.p)

    def channel_width(self, Q, S):
        """Channel width ``b`` (model-specific)."""
        raise NotImplementedError

    def channel_depth(self, S):
        """Channel flow depth ``h`` (model-specific)."""
        raise NotImplementedError
