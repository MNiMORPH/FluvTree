"""
Physics-named long-profile processes on the diffusion solver.

``DiffusionProcess`` drives the :class:`DiffusionSolver`
directly on the canonical :class:`RiverNetwork`: the solver walks the graph's own
``x, z, Q, B`` fields, topology, and boundaries -- there is no separate engine
object to build and no ``z`` to pull/push, because the solver reads and writes the
graph in place. The transport **closure** selects the physics; ``GravelBed``
and ``SandBed`` are thin presets.

"Gravel"/"sand" are conventional shorthand -- the real discriminators are bank
stability (noncohesive grain-threshold vs cohesive mud banks) and form drag; see
the naming note in ``fluvtree.closures``.

This is the in-tree counterpart to :class:`GRLP`, which wraps the external
published ``grlp`` for cross-validation. Here the solver lives inside FluvTree and
a closure picks the model, so gravel and sand are the *same* process.

Time integration is **BDF2** (second-order, as in GRLP), self-started with one
lower-order step. The nonlinear conductance is relinearized by Picard iteration: by
default iterate-to-convergence (``tol``) -- GRLP's default; pass ``niter`` for a
fixed count instead. Note sand still needs a sloped initial ``z`` (its conductance
is singular at ``S = 0``; see the parked issue).
"""

from fluvtree.processes.base import Process
from fluvtree.solvers.diffusion import DiffusionSolver
from fluvtree.closures.gravel import GravelClosure
from fluvtree.closures.sand import SandClosure


class DiffusionProcess(Process):
    """
    Long-profile evolution with the diffusion solver and a transport closure.

    Parameters
    ----------
    network : RiverNetwork, optional
        Canonical network carrying ``x, z, Q, B`` (edges) and ``S0`` (head nodes),
        ``x_bl``/``z_bl`` (outlet node); optional ``Q_s_0`` (graph attribute). Omit
        to bind later via a :class:`FluvTree` model or :class:`Scheduler`.
    closure : TransportClosure
        The physics (e.g. ``GravelClosure`` or ``SandClosure``).
    tol : float, optional
        Picard convergence tolerance (default ``1e-4``): each step iterates the
        nonlinear conductance until ``max|z_k - z_{k-1}| < tol``. Ignored when
        ``niter`` is given.
    niter : int, optional
        Fixed number of Picard iterations per step instead of iterating to ``tol``
        (matches GRLP's ``set_niter``; useful for bit-for-bit comparison).
    gravel_attrition : float, optional
        Turn on Sternberg downstream fining of the gravel load: the fractional load
        loss per km (``fluvtree.common.gravel_attrition``). Off by default. The
        fining sink relinearizes each Picard iterate.
    """

    reads = ("x", "z", "Q", "B", "S0", "x_bl", "z_bl")
    writes = ("z",)

    def __init__(self, network=None, closure=None, tol=1.0e-4, niter=None,
                 gravel_attrition=None):
        self.closure = closure
        self._tol = None if niter is not None else tol
        self._niter = 3 if niter is None else int(niter)
        self._gravel_attrition = gravel_attrition
        self.solver = None
        super().__init__(network)

    def _on_bind(self):
        self.solver = DiffusionSolver(self.network, self.closure)

    def step(self, dt, nt=1):
        """Advance the long profile ``nt`` steps of ``dt`` [s], in place on the graph."""
        self._require_bound()
        dyn = None
        if self._gravel_attrition is not None:
            from fluvtree.common.gravel_attrition import dynamic_source
            dyn = dynamic_source(self.network, self.closure, self._gravel_attrition,
                                 intermittency=self.solver.intermittency,
                                 sinuosity=self.solver.sinuosity)
        self.solver.evolve(nt=nt, dt=dt, niter=self._niter, tol=self._tol,
                           dynamic_source=dyn)


class GravelBed(DiffusionProcess):
    """Gravel-bed long-profile process (noncohesive banks, no form drag).

    Thin preset over :class:`DiffusionProcess` with a :class:`GravelClosure`. Grain
    size ``D`` is optional (needed only to resolve width/depth, not to evolve).
    ``gravel_attrition`` (fractional load loss per km) turns on Sternberg fining."""

    def __init__(self, network=None, D=None, tol=1.0e-4, niter=None,
                 gravel_attrition=None, **closure_kwargs):
        super().__init__(network, GravelClosure(D=D, **closure_kwargs),
                         tol=tol, niter=niter, gravel_attrition=gravel_attrition)


class SandBed(DiffusionProcess):
    """Sand-bed long-profile process (cohesive mud banks + form drag).

    Thin preset over :class:`DiffusionProcess` with a :class:`SandClosure`. Requires
    grain size ``D``, Manning's ``n``, and the critical bank stress
    ``tau_crit_bank``. Use a sloped initial ``z`` (sand is singular at ``S=0``)."""

    def __init__(self, network=None, D=None, n=None, tau_crit_bank=None,
                 tol=1.0e-4, niter=None, **closure_kwargs):
        super().__init__(
            network,
            SandClosure(D=D, n=n, tau_crit_bank=tau_crit_bank, **closure_kwargs),
            tol=tol, niter=niter)
