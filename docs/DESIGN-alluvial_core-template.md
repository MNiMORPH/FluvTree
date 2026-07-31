# Design: the power-law nonlinear-diffusion solver (`fluvtree.solvers.diffusion`)

Status: **built** (branch `rebuild/alluvial-template`). Purpose-built (not GRLP
trimmed); GRLP's validated numerics were lifted *surgically and additively* into a
clean, minimal structure. `grlp`/`srlp` stay **external** validation oracles, never
vendored.

Outcome: `fluvtree/solvers/diffusion.py` (`DiffusionSolver`) walks the shared
`RiverNetwork` directly and reproduces external `grlp` bit-for-bit in backward
Euler (~1e-13) on a chain, a 1-into-1 series, and a multi-tributary confluence, and
matches `srlp`'s steady state; `fluvtree/closures/base.py` is the
`TransportClosure` interface only, with concrete closures in `fluvtree/closures/gravel.py` / `sand.py`. The 2258-line vendored engine is gone; `fluvtree.solvers.diffusion` depends on
`numpy` + `scipy.sparse` alone. It integrates in backward Euler with constant `B`.
**Fidelity note (2026-07-31): this is currently an *unfaithful* port.** GRLP's
default is BDF2 + iterate-to-convergence; BDF2 and the Sternberg gravel-abrasion
sink are not ported and were mis-framed here as "later add-ons" — they are
GRLP-parity gaps to close, not a designed scope (see "GRLP-parity gaps" below).
The rest of this note is the as-built design.

## What it is

A general implicit solver for the transport-limited long-profile PDE class -- a
**power-law nonlinear diffusion** -- on the shared `network`:

    dz/dt = (1 / ((1 - lambda_p) * B)) * dQ_s/dx ,   Q_s = k_Qs * I * Q * S^p

`S = downstream slope`. Because `Q_s` is a power law of slope, this is a nonlinear
diffusion of `z`; the exponent `p` and coefficient `k_Qs` come from a **closure**.
Gravel (`p = 7/6`) and sand (`p = 5/6`) are two closures; the solver hard-codes
none of it. This is the "template": one solver, closures select the model.

## What it depends on / provides

- **In:** the shared `network` (reaches, topology, traversal); a `TransportClosure`
  (already exists: `k_Qs`, `p`, conductance exponent); boundaries (head supply
  `S0`/`Q_s_0`, base level).
- **Out:** the evolved `z`. No own network, no user-facing wrapper.
- **Deps:** `numpy`, `scipy.sparse` only. (No matplotlib / scipy.optimize /
  scipy.stats / networkx -- those served GRLP's periphery.)

## Structure (target ~600-800 lines total)

`fluvtree.solvers.diffusion` is the **template** and knows nothing about gravel vs sand. The
concrete closures live *outside* it, in thin modules that call it (`gravel`,
`sand` -- the grlp/srlp thin shims); adding a transport law is a new module, never
a core edit.

- `closures.py` -- **only the `TransportClosure` interface** (the base class /
  contract: `p`, `k_Qs`, conductance exponent). The template's plug socket.
  `GravelClosure` / `SandClosure` do **not** live here -- they are thin external
  modules (`gravel`, `sand`) that subclass this interface and depend on
  `fluvtree.solvers.diffusion`.
- `solver.py` -- the implicit assembly + Picard, on the shared `network`:
  - assemble the global sparse system per Picard iterate (the `p`-parameterized
    interior stencil + the conservative confluence junction cell);
  - solve (`scipy.sparse` spsolve);
  - Picard-iterate to convergence (nonlinear coefficient relinearizes on the
    current iterate).
- `engine.py` (thin) -- orchestration: step / evolve `nt` steps, time integration
  (backward Euler; BDF2 is GRLP's default and is not yet ported -- a fidelity gap,
  see "GRLP-parity gaps"), reading reach state from the shared `network` and
  writing `z` back.

## Surgically lifted from GRLP (the hard-won, validated pieces)

These are lifted *verbatim-in-spirit* into the clean structure and validated bit-
for-bit against `grlp`:

1. **The conservative confluence junction cell** -- the flux-coefficient (not
   Jacobian) three-node coupling that conserves sources across confluences
   (`grlp@366fb3e` `solver.py`; the 6/7 fix).
2. **The interior implicit stencil**, `p`-parameterized: conductance `~ S^(p-1)`,
   the `2p` Jacobian factor, the discharge-gradient term.
3. **Picard iteration** (freeze RHS at start-of-step `z`, relinearize the
   coefficient each iterate) + the iterate-to-convergence control.
4. **Boundary handling** -- head sediment supply (`S0` or `Q_s_0`, via `1/p`) and
   the base-level Dirichlet outlet.

## Written fresh (not lifted)

- **The solver lives *on* the shared `RiverNetwork` directly -- NO parallel reach
  objects, NO own `Network`/`Segment`, NO own graph.** Everywhere GRLP's solver
  reads `seg.z` / `net.segments` / its index maps, the new solver reads off the
  shared `RiverNetwork` (`get_segment_field`, `upstream_segments`,
  `downstream_segment`, node boundary fields) and gets physics (`k_Qs`, `p`,
  conductance exponent, `C0`) from the **closure**, applied per node. This is the
  crux: lift the assembly *math*, rewrite the data *access* to the one network.
  (A "lightweight reach" object list is still a duplicate network -- avoid it.)
- The step/evolve orchestration (thin, purpose-built).

## Dropped entirely (GRLP periphery)

`LongProfile` wrapper; `analytical_threshold_width*`; `slope_area`; plotting;
`compute_channel_width`/`compute_flow_depth` (these belong to the *closures* as
diagnostics if wanted); the dozens of `set_*` user setters; the redundant
`networkx` graph.

## GRLP-parity gaps and genuinely-later items

**GRLP-parity gaps (flagged 2026-07-31) — NOT agreed deferrals. The lift dropped
these and this note mis-framed them as "later"; they are fidelity gaps to close for
a faithful port:**

- **BDF2** — GRLP's *default* time integration; restore as FluvTree's default.
  Constant-`B` BDF2 needs only the three-level time term (`time_diag = 3/2`, RHS
  history `2 zⁿ − ½ zⁿ⁻¹`) + one extra history level; no volume-first required.
- **Sternberg gravel abrasion / downstream fining**, and the full source-term set
  (`ssd`, uplift `U`) wired through the processes.

**Genuinely later (separate work, or the modeller's call):**

- The volume-first transform for spatially/temporally varying `B` (transient valley
  width — Andy's separate valley-realism task; add when a process needs dynamic width).
- Adaptive time stepping — Andy found it a convenience, not a speed-up (a "fail" for
  GRLP); skip unless asked.
- The hydraulic-radius sand `S->0` regularization (parked issue) and its taper.

## Port spec (execution guide for the solver lift)

Lift the assembly *math* from the vendored GRLP solver (`grlp@366fb3e`,
`grlp/solver.py`), written fresh onto the shared `network` + closure. Scope actually
lifted: **backward Euler, constant B** (so the volume-first transform is identity,
`J = 1`) — a *subset* of GRLP. BDF2 (GRLP's default) and the abrasion sink were NOT
lifted; those are the GRLP-parity gaps above, not a planned phasing.

Per-reach inputs needed by the solver (build from `network` reach fields + closure):
`x, z, Q, B`; `C0 = closure.k_Qs * I * dt / ((1 - lambda_p) * sinuosity**p)`;
`p`, `p_conductance`; topology (up/down); boundary (`S0`, `z_bl`, ghosts).

The three assembly rules to lift (formulas verified in the vendored `solver.py`):
- **Interior node** i (neighbours give `x_up,z_up,Q_up` and `x_down,z_down,Q_down`):
  `dx_up=x[i]-x_up; dx_down=x_down-x[i]; dx_2c=x_down-x_up; dQ_2c=Q_down-Q_up`;
  `S=|z_down-z_up|/dx_2c`; `C1=C0 * S**p_conductance * Q[i]/B[i]`;
  `center = -C1/dx_2c * (2p*(-1/dx_up - 1/dx_down)) + 1` (time_diag=1 for BE);
  `left  = -C1/dx_2c * (2p/dx_up  - dQ_2c/Q[i]/dx_2c)`;
  `right = -C1/dx_2c * (2p/dx_down + dQ_2c/Q[i]/dx_2c)`;  `RHS = zold[i] + src[i]`.
- **Head (Neumann via S0):** `right = -C1/dx_2c * 2p*(1/dx_up+1/dx_down)`;
  `RHS += S0 * C1 * (2p/dx_up - dQ_2c/Q[i]/dx_2c)`.
- **Outlet (Dirichlet via z_bl):** `RHS += z_bl * C1/dx_2c * (2p*(1/(x[-1]-x[-2]) +
  1/(x_down-x[-1]))/2 + dQ_2c/Q[i]/dx_2c)`.
- **Confluence junction cell** (multi-tributary head node): the conservative
  three-node cell using the *flux coefficient* (not the Jacobian) --
  `conductance = C0 * Q * (|z_up - z_down|/L_face)**p_conductance / L_face`,
  summed over tributaries and divided by `land_area_around_confluence`. Lift
  verbatim-in-spirit from vendored `solver.py` (the 6/7 fix). Needs the confluence
  land area; compute it from reach widths/spacings as GRLP does.

`2p` here is `2*closure.p`; for gravel `2*(7/6)=7/3` reproduces GRLP bit-for-bit.
Picard: freeze `RHS` at start-of-step `z`; recompute `C1` (via `S`) on the current
iterate each pass; iterate to `max|z_k - z_{k-1}| < tol`. Solve the sparse system
each pass (`scipy.sparse.linalg.spsolve`).

Validate against `grlp` in **euler** scheme (backward Euler + fixed 3 Picard;
`tests/network_helpers.apply_scheme("euler")`) so the time integration matches.
Steady state is scheme-independent, so `grlp`(bdf2) and `srlp` steady states also
work as oracles.

## Validation (the guardrail, at every step)

Reuse the existing benchmarks with `grlp`/`srlp` as **external** oracles:
- **gravel == GRLP** to machine precision on the network topologies;
- **sand == SRLP** steady state (single-channel oracle).
Build additively; each lifted piece is landed only once its oracle test is green.
