"""
Power-law nonlinear-diffusion solver, on the shared network.

The transport-limited long-profile PDE class

    dz/dt = (1 / ((1 - lambda_p) * B)) * dQ_s/dx ,   Q_s = k_Qs * I * Q * S**p

is a nonlinear diffusion of bed elevation ``z`` (the flux is a power law of the
slope ``S``). This module solves it implicitly on a :class:`RiverNetwork`. The
exponent ``p`` and coefficient ``k_Qs`` come from a *closure* (gravel ``p=7/6``,
sand ``p=5/6``); the solver hard-codes no transport law -- it is the template.

**It lives on the network.** The assembly walks the shared ``RiverNetwork``'s own
segment fields (``x, z, Q, B``), topology (upstream/downstream segments), and node
boundaries (head ``S0``, outlet ``z_bl``/``x_bl``), and reads the physics off the
closure per segment. There is no parallel ``Segment``/``Network`` object graph:
the only state the solver holds is the unknown vector ``z`` it is solving for, as
any linear solver must. A river's per-segment arrays carry all their own nodes
(the confluence value is duplicated across the mainstem's first node and each
tributary's last node, separated by a finite face); junction *nodes* on the graph
are topological markers holding boundary scalars, not shared ``z`` degrees of
freedom. That is the alluvial state layout (see ``processes/grlp_process`` note),
not a duplicate network.

The assembly math -- the ``p``-parameterized interior stencil, the conservative
confluence junction cell (flux coefficient, not Jacobian), and the head/outlet
boundaries -- is lifted from GRLP's validated network solver (``grlp@366fb3e``)
and reproduces it bit-for-bit. Only the data *access* is rewritten onto the shared
network.

Time integration is **second-order BDF2** (as in GRLP), self-started with one
lower-order step; ``B`` is constant. Picard iteration relinearizes the nonlinear
conductance on the current iterate while the right-hand-side history is frozen at
the step's start.

Not yet ported from GRLP (both additive; neither affects the constant-``B`` BDF2
solve): the Sternberg gravel-abrasion / downstream-fining sink, and the
volume-first transform for spatially/temporally varying ``B`` (dynamic valley
width).
"""

import warnings

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


class DiffusionSolver(object):
    """
    Implicit power-law nonlinear-diffusion solver on a :class:`RiverNetwork`.

    Parameters
    ----------
    network : RiverNetwork
        Canonical network carrying per-segment arrays ``x, z, Q, B`` (full node
        arrays, GRLP layout) and node boundaries ``S0`` (channel heads),
        ``x_bl``/``z_bl`` (outlet); optional network-level ``Q_s_0`` (graph
        attribute) sets the upstream sediment supply instead of ``S0``.
    closure : TransportClosure
        The transport physics: ``p``, ``conductance_exponent`` (``= p - 1``),
        ``k_Qs``, ``lambda_p``.
    intermittency : float, optional
        Flood intermittency ``I`` in ``Q_s = k_Qs I Q S**p`` (default 1).
    sinuosity : float, optional
        Channel sinuosity (default 1); enters ``C0`` as ``sinuosity**-p`` and the
        head-supply conversion.

    Topology is fixed, so all index maps, node classifications, and confluence
    land areas are built once at construction. Only ``z`` changes during a solve.
    """

    def __init__(self, network, closure, intermittency=1.0, sinuosity=1.0):
        self.network = network
        self.closure = closure
        self.p = float(closure.p)
        self.p_conductance = float(closure.conductance_exponent)
        self.lambda_p = float(closure.lambda_p)
        self.intermittency = float(intermittency)
        self.sinuosity = float(sinuosity)

        segs = list(network.segment_ids)
        self._segs = segs
        self._pos = {s: i for i, s in enumerate(segs)}   # seg id -> position

        # Constant-in-time reach state: references into the network's own arrays.
        self.x = [np.asarray(network.get_segment_field(s, "x"), float) for s in segs]
        self.Q = [np.asarray(network.get_segment_field(s, "Q"), float) for s in segs]
        self.B = [np.asarray(network.get_segment_field(s, "B"), float) for s in segs]
        self._lengths = [len(xi) for xi in self.x]

        # Global node ordering: segment by segment, as in GRLP's assembler.
        self._starts = np.cumsum([0] + self._lengths)[:-1]
        self._n = int(np.sum(self._lengths))

        # Topology, resolved to positions in the segment list.
        self._up = [sorted(self._pos[u] for u in network.upstream_segments(s))
                    for s in segs]
        self._down = []
        for s in segs:
            d = network.downstream_segment(s)
            self._down.append(None if d is None else self._pos[d])

        # Boundaries. Head slope S0 per channel head; base level at the outlet.
        # A network-level Q_s_0 overrides S0 (converted per head via the closure).
        Q_s_0 = network.graph.graph.get("Q_s_0")
        self._S0 = [None] * len(segs)
        for s in segs:
            i = self._pos[s]
            if not self._up[i]:                     # channel head
                u_node = network.edge_of(s)[0]
                if Q_s_0 is not None:
                    self._S0[i] = self._S0_from_Qs0(Q_s_0, self.Q[i][0])
                else:
                    self._S0[i] = float(network.get_node_field(u_node, "S0"))
        outlet_seg = network.mouth_segments()[0]
        outlet_node = network.edge_of(outlet_seg)[1]
        self._mouth = self._pos[outlet_seg]
        self.z_bl = float(network.get_node_field(outlet_node, "z_bl"))
        try:
            self.x_bl = float(network.get_node_field(outlet_node, "x_bl"))
        except KeyError:
            self.x_bl = None

        self._check_confluence_sizes()
        self._land_area = self._compute_land_areas()

        # C0 = k_Qs I / ((1 - lambda_p) sinuosity**p) * dt, per unit dt (times dt
        # in the assembler). Uniform closure -> one scalar for the whole network.
        self._C0_per_dt = (closure.k_Qs * self.intermittency
                           / ((1.0 - self.lambda_p) * self.sinuosity ** self.p))

        # Two-level BDF2 history, persisted across evolve() calls so the scheme
        # engages even when the scheduler advances one step at a time. None until
        # the solver has taken its self-starting first step.
        self._zold2 = None      # z^{n-1}
        self._dt_prev = None    # the previous step's dt (for variable-step weights)

    # ------------------------------------------------------------------ #
    # One-time setup
    # ------------------------------------------------------------------ #

    def _S0_from_Qs0(self, Q_s_0, Q0):
        """Head boundary slope equivalent to an upstream sediment supply ``Q_s_0``
        (inverse of ``Q_s = k_Qs I Q S**p`` at the head), sinuosity external."""
        return (np.sign(Q0) * self.sinuosity
                * (np.abs(Q_s_0) / (self.closure.k_Qs * self.intermittency
                                    * np.abs(Q0))) ** (1.0 / self.p))

    def _check_confluence_sizes(self):
        """The three-node junction cell reaches to the confluence's second
        interior node and each tributary's second-to-last node, so adjacent
        segments must be long enough. Fail clearly, not with an IndexError."""
        for i, s in enumerate(self._segs):
            if len(self._up[i]) > 1:
                if self._lengths[i] < 3:
                    raise ValueError(
                        "confluence segment %r needs >= 3 nodes (has %d)"
                        % (s, self._lengths[i]))
                for j in self._up[i]:
                    if self._lengths[j] < 2:
                        raise ValueError(
                            "tributary segment %r into confluence %r needs >= 2 "
                            "nodes (has %d)"
                            % (self._segs[j], s, self._lengths[j]))

    def _compute_land_areas(self):
        """Confluence land area per segment: the sum over tributaries of
        (half the gap to the confluence) times the tributary's end width, plus
        (half the first mainstem cell) times the mainstem's first width -- GRLP's
        ``compute_land_areas_around_confluences``. ``None`` where not a confluence."""
        land = [None] * len(self._segs)
        for i in range(len(self._segs)):
            if len(self._up[i]) <= 1:
                continue
            above = 0.0
            for j in self._up[i]:
                half_dx = (self.x[i][0] - self.x[j][-1]) / 2.0
                above += half_dx * self.B[j][-1]
            half_dx_below = (self.x[i][1] - self.x[i][0]) / 2.0
            below = half_dx_below * self.B[i][0]
            land[i] = above + below
        return land

    # ------------------------------------------------------------------ #
    # Assembly (lifted math, network data access)
    # ------------------------------------------------------------------ #

    def _face_conductance(self, z_up, z_down, Q, x_up, x_down, C0):
        """Shared junction-face sediment-flux coefficient: ``conductance * dz =
        k_Qs I Q (S/sinuosity)**p`` at Picard convergence. This is the flux
        coefficient itself, NOT the linearized (p) Jacobian dQ_s/dS -- using the
        Jacobian under-applies distributed sources by exactly ``1/p`` at every
        confluence (the 6/7 fix)."""
        L_face = x_down - x_up
        return (C0 * Q * (np.abs(z_up - z_down) / L_face) ** self.p_conductance
                / L_face)

    def _assemble(self, Z, dt, C0, src, time_diag):
        """Build the global sparse LHS and RHS for one Picard iterate.

        ``Z`` is the current-iterate elevation (list of per-segment arrays); ``src``
        is the frozen right-hand-side history (the BDF2 elevation history plus
        source terms) and ``time_diag`` is the time-derivative diagonal (3/2 for
        uniform BDF2, the variable-step value otherwise, 1 for the startup step).
        Constant B, so the volume-first transform is the identity and is omitted."""
        starts, lengths = self._starts, self._lengths
        n = self._n
        p, twop = self.p, 2.0 * self.p
        rows, cols, vals = [], [], []
        RHS = np.zeros(n)
        for si in range(len(self._segs)):
            offset = starts[si]
            L = lengths[si]
            z, x, Q, B = Z[si], self.x[si], self.Q[si], self.B[si]
            up_ids, down_id = self._up[si], self._down[si]
            for i in range(L):
                g = offset + i
                is_confluence = (i == 0 and len(up_ids) > 1)
                down_is_confluence = (
                    i == L - 1 and down_id is not None
                    and len(self._up[down_id]) > 1)
                up_is_confluence = (i == 1 and len(up_ids) > 1)

                # ===== multi-tributary junction: shared-flux three-node cell ====
                if is_confluence:
                    A = self._land_area[si]
                    cond_down = self._face_conductance(
                        z[0], z[1], 0.5 * (Q[0] + Q[1]), x[0], x[1], C0)
                    cond_sum = cond_down
                    rows.append(g); cols.append(g + 1)
                    vals.append(-cond_down / A)
                    for j in up_ids:
                        upg = starts[j] + lengths[j] - 1
                        cond_up = self._face_conductance(
                            Z[j][-1], z[0], self.Q[j][-1],
                            self.x[j][-1], x[0], C0)
                        cond_sum += cond_up
                        rows.append(g); cols.append(upg)
                        vals.append(-cond_up / A)
                    rows.append(g); cols.append(g)
                    vals.append(time_diag + cond_sum / A)
                    RHS[g] = src[si][i]
                    continue
                if down_is_confluence:
                    dj = down_id
                    downg = starts[dj]
                    land = B[-1] * 0.5 * ((x[-1] - x[-2])
                                          + (self.x[dj][0] - x[-1]))
                    cond_downseg = self._face_conductance(
                        z[-1], Z[dj][0], Q[-1], x[-1], self.x[dj][0], C0)
                    cond_up = self._face_conductance(
                        z[-2], z[-1], 0.5 * (Q[-2] + Q[-1]), x[-2], x[-1], C0)
                    rows.append(g); cols.append(g - 1)
                    vals.append(-cond_up / land)
                    rows.append(g); cols.append(downg)
                    vals.append(-cond_downseg / land)
                    rows.append(g); cols.append(g)
                    vals.append(time_diag + (cond_up + cond_downseg) / land)
                    RHS[g] = src[si][i]
                    continue
                if up_is_confluence:
                    land = B[1] * 0.5 * ((x[1] - x[0]) + (x[2] - x[1]))
                    cond_up = self._face_conductance(
                        z[0], z[1], 0.5 * (Q[0] + Q[1]), x[0], x[1], C0)
                    cond_down = self._face_conductance(
                        z[1], z[2], 0.5 * (Q[1] + Q[2]), x[1], x[2], C0)
                    rows.append(g); cols.append(g - 1)
                    vals.append(-cond_up / land)
                    rows.append(g); cols.append(g + 1)
                    vals.append(-cond_down / land)
                    rows.append(g); cols.append(g)
                    vals.append(time_diag + (cond_up + cond_down) / land)
                    RHS[g] = src[si][i]
                    continue

                # ----- upstream neighbour (or head ghost) -----
                if i > 0:
                    up_g = g - 1
                    z_up, x_up, Q_up = z[i - 1], x[i - 1], Q[i - 1]
                    is_head = False
                elif len(up_ids) == 0:
                    is_head = True
                    up_g = None
                    x_up = 2 * x[0] - x[1]
                    z_up = z[0] + self._S0[si] * (x[0] - x_up)
                    Q_up = 2 * Q[0] - Q[1]
                else:                               # single upstream segment
                    is_head = False
                    j = up_ids[0]
                    up_g = starts[j] + lengths[j] - 1
                    z_up, x_up, Q_up = Z[j][-1], self.x[j][-1], self.Q[j][-1]

                # ----- downstream neighbour (or outlet ghost) -----
                if i < L - 1:
                    down_g = g + 1
                    z_down, x_down, Q_down = z[i + 1], x[i + 1], Q[i + 1]
                    is_outlet = False
                elif down_id is None:
                    is_outlet = True
                    down_g = None
                    x_down = self.x_bl if self.x_bl is not None \
                        else 2 * x[-1] - x[-2]
                    z_down = self.z_bl
                    Q_down = 2 * Q[-1] - Q[-2]
                else:
                    is_outlet = False
                    dj = down_id
                    down_g = starts[dj]
                    z_down, x_down, Q_down = Z[dj][0], self.x[dj][0], self.Q[dj][0]

                # ----- interior stencil -----
                dx_up = x[i] - x_up
                dx_down = x_down - x[i]
                dx_2c = x_down - x_up
                dQ_2c = Q_down - Q_up
                S = np.abs(z_down - z_up) / dx_2c
                C1 = C0 * S ** self.p_conductance * Q[i] / B[i]
                center = -C1 / dx_2c * (twop * (-1 / dx_up - 1 / dx_down)) + time_diag
                left = -C1 / dx_2c * (twop / dx_up - dQ_2c / Q[i] / dx_2c)
                right = -C1 / dx_2c * (twop / dx_down + dQ_2c / Q[i] / dx_2c)
                rhs_g = src[si][i]
                if is_head:                         # sediment-flux Neumann
                    right = -C1 / dx_2c * twop * (1 / dx_up + 1 / dx_down)
                    rhs_g += self._S0[si] * C1 * (twop / dx_up
                                                  - dQ_2c / Q[i] / dx_2c)
                if is_outlet:                       # base-level Dirichlet
                    rhs_g += z_down * C1 / dx_2c * (
                        twop * (1 / (x[-1] - x[-2]) + 1 / (x_down - x[-1])) / 2.0
                        + dQ_2c / Q[i] / dx_2c)
                rows.append(g); cols.append(g)
                vals.append(center)
                if up_g is not None:
                    rows.append(g); cols.append(up_g)
                    vals.append(left)
                if down_g is not None:
                    rows.append(g); cols.append(down_g)
                    vals.append(right)
                RHS[g] = rhs_g
        LHS = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
        return LHS, RHS

    # ------------------------------------------------------------------ #
    # Time stepping
    # ------------------------------------------------------------------ #

    def _pull_z(self):
        """Read the current elevation off the network into working arrays."""
        return [np.asarray(self.network.get_segment_field(s, "z"), float).copy()
                for s in self._segs]

    def _push_z(self, Z):
        """Write the solved elevation back onto the network."""
        for si, s in enumerate(self._segs):
            self.network.set_segment_field(s, "z", Z[si].copy())

    def evolve(self, nt, dt, niter=3, tol=None, max_iter=100):
        """
        Advance ``nt`` BDF2 steps of ``dt`` [s], writing ``z`` back to the network.

        Time integration is **second-order BDF2** (as in GRLP), self-started with a
        single lower-order step the first time the solver runs -- a two-step method
        has no prior state to reach back to. The two-level history persists across
        ``evolve`` calls, so BDF2 engages from the second step onward even when the
        scheduler advances one step at a time; variable steps are exact (the weights
        use ``omega = dt / dt_prev``).

        Picard iteration relinearizes the nonlinear conductance on the current
        iterate while the RHS history is frozen at the step's start. With
        ``tol=None`` (default) it runs a fixed ``niter`` iterations (bit-for-bit
        comparable to GRLP's ``set_niter``); with ``tol`` set it iterates to
        ``max|z_k - z_{k-1}| < tol`` (``max_iter`` is the safety cap).
        """
        nseg = len(self._segs)
        C0 = self._C0_per_dt * dt
        # per-node source (m/s) -> elevation increment; optional "source" field.
        src_rate = []
        for s in self._segs:
            try:
                r = np.asarray(self.network.get_segment_field(s, "source"), float)
            except KeyError:
                r = 0.0
            src_rate.append(r)
        Z = self._pull_z()
        for _ in range(int(nt)):
            zold = [zi.copy() for zi in Z]              # z^n, frozen for this step
            if self._zold2 is not None:                 # BDF2: three-level history
                w = dt / self._dt_prev
                b, c = 1.0 + w, w * w / (1.0 + w)
                time_diag = (1.0 + 2.0 * w) / (1.0 + w)
                z_rhs = [b * zold[si] - c * self._zold2[si] for si in range(nseg)]
            else:                                       # self-start: one Euler step
                time_diag = 1.0
                z_rhs = zold
            src = [z_rhs[si] + src_rate[si] * dt for si in range(nseg)]
            converged = tol is None
            cap = int(niter) if tol is None else int(max_iter)
            change = 0.0
            for _k in range(cap):
                LHS, RHS = self._assemble(Z, dt, C0, src, time_diag)
                out = spsolve(LHS, RHS)
                change = 0.0
                for si in range(nseg):
                    znew = out[self._starts[si]:self._starts[si] + self._lengths[si]]
                    if tol is not None:
                        change = max(change, float(np.max(np.abs(znew - Z[si]))))
                    Z[si] = znew
                if tol is not None and change < tol:
                    converged = True
                    break
            if not converged:
                warnings.warn(
                    "Picard did not converge to tol=%g in %d iterations "
                    "(last change %g m); result may be under-converged."
                    % (tol, max_iter, change), RuntimeWarning)
            self._zold2 = zold                          # z^{n-1} for the next step
            self._dt_prev = dt
        self._push_z(Z)
        # advance the network clock, as the scheduler expects
        self.network.graph.graph["t"] = \
            self.network.graph.graph.get("t", 0.0) + int(nt) * dt
