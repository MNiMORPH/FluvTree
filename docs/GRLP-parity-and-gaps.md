# GRLP → FluvTree parity: what did *not* make it across

**Status: audit, 2026-08-01.** A full method-by-method comparison of the GRLP
package (`grlp/grlp.py` — `Segment`/`LongProfile`/`Network` — plus `grlp/solver.py`
and `grlp/build_synthetic_network.py`) against FluvTree's `src/fluvtree/`.

**Read the scope first, or the list misleads.** FluvTree deliberately ported the
GRLP **solver engine** — evolve a bed profile `z` on a river network,
transport-limited, with sources — and *generalized* it (closure-driven exponent →
gravel *and* sand; canonical `RiverNetwork`; process/scheduler/model layers). By
that scope the engine is a faithful, bit-for-bit port. Most of what is listed
below as "not ported" is **GRLP's gravel-threshold-width *science* layer**
(analytical solutions, linear-response transfer functions, network morphometry) —
capability that GRLP has and FluvTree does not, but that was never inside the
engine port's frame and may or may not belong in FluvTree's core. This document
is a *factual inventory of absence*, **not** a to-do list: scope decisions are
Andy's.

The earlier in-session solver audit (transcript row 3694) was complete *for
`solver.py`* — it flagged every engine gap below. Its unstated limit was that it
never enumerated `grlp.py`'s `Segment`/`LongProfile`/`Network` science surface;
that surface is the bulk of this document.

---

## A. Solver engine (`grlp/solver.py`) — genuine engine gaps

All four were flagged in the row-3694 audit.

1. **Volume-first transform** — *the big one.* GRLP's `assemble` row-scales the
   elevation system by the storage Jacobian `J = dV/dz`, carries the BDF2 history
   in volume space (`Vhist`), and adds the nonlinear-map linearization correction
   `Vcorr` (`solver.py:70–138, 325–336`). FluvTree's `_assemble` is the **z-based
   rectangular collapse** — mathematically identical for constant `B` (hence
   bit-for-bit validation), but it omits the transform entirely
   (`diffusion.py:206`). This is the foundation GRLP built for dynamic `B`
   (issue #19); without it, FluvTree cannot do a widening/narrowing valley.
   *Requires, together with it:* the geometry primitives in §B.

2. **Adaptive time-stepping** — `evolve_adaptive` + `_trial_step` (step-doubling,
   I-controller, reject/retry) (`solver.py:423–555`) and its control setter
   `set_adaptive_timestep`. Absent in FluvTree (`evolve` is fixed-`dt` only). Note
   Andy has characterized adaptive stepping as a convenience, not a speed-up.

3. **`Q_ghost_upstream` / `Q_ghost_downstream` boundary overrides** — GRLP lets an
   explicit ghost discharge be set at head/outlet (`solver.py:249–252, 283–286`);
   FluvTree always linearly extrapolates `Q` at boundaries. Minor.

4. **`dz_dt` per-step diagnostic** — GRLP stamps `seg.dz_dt = (z − zold)/dt` each
   step (`solver.py:420, 521`); FluvTree does not expose it. Minor. (Note this is
   also the hook GRLP's dynamic-`B` design updates `B` from, between steps.)

---

## B. Segment geometry & hydraulics (`grlp.py` `Segment`)

- **`valley_width(z)` / `storage_jacobian(z)` / `storage_volume(z)`**
  (`grlp.py:352–379`) — the valley-storage geometry primitives the volume-first
  transform (§A1) consumes. **Not ported** (FluvTree closures carry no storage
  geometry). Rectangular default only, implicit.
- **`set_B` as a *dynamic* width** — spatially varying constant-in-time `B` **is**
  supported (FluvTree reads a per-node `B` field). Time/`z`-varying `B` is **not**
  (that is §A1). Listed here because "dynamic B" is often thought of as "call
  `set_B`," which it is not.
- **`compute_channel_width` / `compute_flow_depth`** — **reimplemented**, not
  absent: FluvTree closures expose `channel_width(Q, S)` / `channel_depth(S)`
  (gravel and sand forms). Equivalent capability, different home.
- **Segment setters** (`set_uplift_rate`, `set_source_sink_distributed`,
  `set_Sternberg_gravel_loss`, `set_Qs_input_upstream`, `set_S0`, `set_bl`/
  `set_z_bl`/`set_x_bl`, `set_intermittency`) — **functionally covered** via
  network fields, graph attributes, and process parameters (`uplift=`,
  `source_sink=`, `gravel_attrition=`, the `Q_s_0` graph attribute, `S0`/`z_bl`/
  `x_bl` node fields). Different API, same physics; all validated bit-for-bit.

---

## C. Analytical solutions (`grlp.py` `LongProfile`) — NOT ported

Independent, closed-/semi-analytical steady-state profiles — the *external*
validation baselines (FluvTree has so far validated against GRLP's **numerics**,
not against these).

- `analytical_threshold_width`
- `analytical_threshold_width_perturbation`
- `analytical_threshold_width_uplift`

Gravel-threshold-width–specific. **Highest near-term value of the un-ported
surface** — they are the principled, GRLP-independent check on the engine.

---

## D. Linear-response / transfer-function diagnostics (`grlp.py` `LongProfile`) — NOT ported

Wickert & Schildgen (2019) spectral / linear-systems theory for the gravel
long-profile response:

- `compute_diffusivity`
- `compute_equilibration_time`
- `compute_e_folding_time`
- `compute_wavenumber`
- `compute_series_coefficient`
- `compute_z_series_terms` / `compute_Qs_series_terms`
- `compute_z_gain` / `compute_Qs_gain`
- `compute_z_lag` / `compute_Qs_lag`

Gravel-specific analytical theory; none present in FluvTree.

---

## E. Network morphometry / geomorphometry (`grlp.py` `Network`) — NOT ported

Quantitative drainage-network structure metrics:

- `compute_absolute_lengths`
- `compute_topological_lengths`
- `compute_topological_widths`
- `compute_strahler_orders`
- `compute_horton_ratios`
- `compute_tokunaga_metrics`
- `compute_jarvis_E`
- `compute_mean_diffusivity`
- `compute_network_properties`
- `find_hack_parameters` / `find_hack_parameters_non_dim`

None present in FluvTree.

---

## F. Synthetic / random network generation (`grlp/build_synthetic_network.py`) — NOT ported

Random-topology test-network builders:

- `generate_random_network`, `Shreve_Random_Network`
- `generate_x_domain`, `generate_discharges`, `generate_ssds`,
  `generate_variable_widths`, `generate_zs`
- topology recursion helpers (`_downstream_IDs`, `_upstream_IDs`, …)

FluvTree builds networks only by **explicit construction**
(`RiverNetwork.from_segment_lists`, plus the `build_grlp_network` test helper). No
random/synthetic generation.

---

## G. Network construction internals — mostly reimplemented; one open item

GRLP's `Network` construction/topology methods (`build_graph`, `initialize`,
`find_downstream_IDs`/`find_upstream_IDs`, `create_list_of_channel_head/
mouth_segment_IDs`, `compute_land_areas_around_confluences`, `build_ID_list`,
`get_z_lengths`) are **reimplemented** in FluvTree's `RiverNetwork` + the solver
(`from_segment_lists`, `upstream_segments`/`downstream_segment`,
`head_segments`/`mouth_segments`, `_compute_land_areas`). Equivalent capability.

- **`update_Q` — discharge accumulation from the network — NOT clearly ported.**
  GRLP can accumulate `Q` downstream through the topology; FluvTree currently
  **reads `Q` as a given per-node field** and does not derive it from drainage
  area / accumulation. This overlaps the parked **A → Q** item in the
  DEM→RiverNetwork pipeline work. Flagging it as a real gap, not an internal.

---

## H. For contrast — what *was* ported (so this list isn't read as mostly-absent)

The engine core, all validated bit-for-bit against GRLP (and, where applicable,
SRLP):

- Interior implicit stencil — **ported + generalized** (closure-driven exponent →
  gravel *and* sand)
- Confluence junction cell, all three cases (flux-coefficient / "6-7" fix)
- Head Neumann (`S0` / `Q_s_0`) + outlet Dirichlet (`z_bl`) boundaries
- Picard nonlinear iteration — fixed `niter` **and** iterate-to-tolerance
- **BDF2** (second-order, variable-step, self-started) — the shipped default
- `compute_Q_s` sediment-flux network walk (`fluvtree.common.gravel_attrition`)
- **Sternberg gravel attrition / downstream fining** (`gravel_attrition`)
- **Tectonic uplift / subsidence** as a coupled source (`tectonics`)
- **Generic distributed source/sink** (`source_sink`)
- `C0`, sinuosity, intermittency (parameterized)
- Long-profile / slope–area / planform plotting (reimplemented in `fluvtree.plot`)

---

## Summary of the true gaps, by kind

| Kind | Items | Likely scope call |
|---|---|---|
| **Engine — dynamic B** | volume-first transform (§A1) + geometry primitives (§B) | in scope; the valley-realism foundation |
| **Engine — convenience** | adaptive stepping, `Q_ghost`, `dz_dt` (§A2–4) | optional; adaptive de-prioritized by Andy |
| **Validation baselines** | 3 analytical solutions (§C) | high value; independent engine check |
| **Gravel linear theory** | ~10 transfer-function diagnostics (§D) | GRLP-context; maybe not core |
| **Network morphometry** | ~11 metrics (§E) | GRLP-context; maybe not core |
| **Synthetic networks** | random-network generators (§F) | test/util; maybe not core |
| **Discharge accumulation** | `update_Q` / A→Q (§G) | real gap; ties to DEM pipeline |
