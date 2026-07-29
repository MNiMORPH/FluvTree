# FluvTree engine architecture — canonical state + processes

Status: **design settled (option A); to be built.**
Author: Andy Wickert + Claude Code. Date: 2026-07-29.

## Guiding principle

> **FluvTree holds what is. Its methods modify what is.**

The network's state is a single canonical object — a directed convergent graph.
Every physical model (GRLP, TerraPIN, stream-power incision, a user rule) is a
**process**: a method that reads named fields off the graph, computes, and writes
named fields back. No process owns state; the graph owns state, and processes
operate on it. This is what makes an *arbitrary* ruleset possible — adding a new
process is writing one class, not wiring new pairwise sync between models.

## The topology is fixed

Topology does not change through time (no capture, avulsion, or divide migration).
This is the enabling assumption: every process's mapping from graph structure to
its own working indices is built **once** and stays valid for the whole run. GRLP's
global sparse-matrix index map, TerraPIN's one-cross-section-per-node array — all
of it can be cached. Nothing renumbers.

## Canonical state: what lives where

The graph is a `networkx.DiGraph`, edges point downstream.

- **Edges are segments (reaches).** Each carries the along-reach interior arrays:
  `x, y, z, Q`, and whatever a process needs (`B`, sediment, grain size, …). These
  are the *interior* nodes of the reach; the reach's two endpoints live on the
  graph nodes it connects.
- **Nodes are junctions** (channel heads, confluences, the outlet). Each carries
  the scalar state at that point: `x, y, z, Q, s`, and the boundary role
  (source / confluence / outlet).

**Reach-assembly convention.** A full profile over a reach is
`hstack(upstream node, edge interior array, downstream node)` — the pattern
`setupDomain.py` already uses. Endpoint values are held **once**, on the node, and
shared by every edge and every process that touches that junction. That shared node
`z` is exactly how two processes (GRLP's confluence cell and TerraPIN's per-node
cross-section) stay consistent for free.

## The process interface

A process is an operator over the graph:

```python
class Process:
    reads  = (...)   # named fields it consumes   (for scheduling / validation)
    writes = (...)   # named fields it produces
    def step(self, G, dt):
        ...          # read G's fields, compute, write G's fields back
```

Processes may hold **parameters and cached index maps** (fixed topology makes the
cache safe), but never the canonical physical state — that is always read fresh
from `G` and written back to `G`.

## The scheduler

An ordered list of processes, run each `dt`. This is the generalization of the
hand-written GRLP↔TerraPIN loop in `terrapin-grlp-coupling.md` — that loop *is* a
two-process schedule; here the sequence is data, not code:

```python
for _ in range(nt):
    for proc in ruleset:      # ordered list of Process objects
        proc.step(G, dt)
    t += dt
```

Order matters and is the modeller's to set (lateral before vertical, etc.). Because
processes declare `reads`/`writes`, the scheduler can later validate an ordering or
warn on a read-before-write — not needed for v1, but the declaration earns its place.

## The three doors, as processes

The earlier "three doors" were never three variants of one solver; under (A) they
are three peer processes on the shared graph:

| Door | Process | Numerics | Reuses |
|---|---|---|---|
| 1 — transport-law variants | GRLP-family | GRLP's **implicit** tridiagonal solve, `_face_conductance` swapped | GRLP matrix machinery + graph |
| 2 — arbitrary user rule | `ExplicitRule` | forward-integrate a supplied `dz/dt` | graph substrate only (CFL-bound) |
| 3 — stream-power incision | `StreamPower` | its **own** explicit/upwind scheme | graph substrate only |
| (lateral) | TerraPIN | event-driven per-node cross-section | graph substrate + `terrapin` |

The trap to avoid: do **not** force door 3 into GRLP's diffusion solver. Each
process owns its own numerics; the engine only owns the graph and the schedule.
GRLP's implicit large-`dt` stability is a property of the GRLP process, not of the
engine, and it survives only for the diffusion family (door 1).

**Three tiers, not two.** The extraction seam is not a clean "topology = generic,
everything else = GRLP." There is a middle tier: some of the *river machinery* —
the implicit tridiagonal assembly, the confluence junction cell, the Picard/BDF2
time integration, the ghost/boundary handling — is reusable by **any
transport-limited rule**, not just gravel. Only the flux law itself
(`_face_conductance` + `C0`, the `k_Qs · Q · S^{7/6}` closure) is GRLP-specific.
So the layering is: (1) **generic substrate** — topology + traversal + state,
no physics; (2) **transport-limited engine** — the implicit machinery, parameterized
by a supplied flux closure; (3) **GRLP closure** — the gravel flux law that plugs
into tier 2. Door-1 rules reuse tier 2 with a different closure; doors 2–3 sit only
on tier 1. Tier 2 is worth factoring out of GRLP eventually, but tier 1 is the
first and cleanest extraction.

## GRLP as a process under (A)

`grlp.Network` becomes an **adapter**, not a state owner:

1. On construction, read topology from `G`, build the global index map **once**.
2. Each `step(G, dt)`: borrow `z, Q, B` arrays from `G`'s edges/nodes, assemble and
   Picard-solve exactly as today (`solver.evolve` for one step), write `z` back to
   `G`.

GRLP's numerics are untouched. What changes is ownership: `z` no longer lives on
`Segment`; it lives on the graph and is borrowed for the solve. The fixed-topology
guarantee is what lets the index map be cached across steps.

## Sediment double-counting (carried over from the TerraPIN note)

Still holds and is now a scheduling property: GRLP's vertical `dz` already carries
the vertical sediment budget; TerraPIN's *lateral* `sediment_out` is the only new
contribution feeding GRLP's `ssd`. Read `sediment_out` after a **lateral** op,
discard it after the vertical `sweep`. See `terrapin-grlp-coupling.md`.

## Extraction & provenance

The generic substrate (tier 1) is **extracted, not grafted.** GRLP's network code
never lived in its own file — `Network` is a class inside a single large
`grlp/grlp.py`, its commits interleaved with physics — so a faithful per-file
history transplant is not available, and `git filter-repo` would only import the
whole file's history (physics and all) as noise. Therefore:

- **GRLP stays the living record** of that code's evolution; its history remains in
  GRLP where the file lives.
- **FluvTree gets a fresh extraction commit** with explicit provenance: extracted
  from `grlp@366fb3e` (v2.1.0). Record the source SHA in the commit and here.
- Hard-won rationale (why the junction cell uses the flux coefficient not the
  Jacobian; why de-pad) is carried as **code comments**, not as a git graft.

### Dependency direction — the guardrail

GRLP is published (RTD, PyPI, Zenodo DOIs). FluvTree is early. **Published GRLP must
not depend on prototype FluvTree** — that would chain a stable, citable release to a
moving target. So:

- The substrate's **home is FluvTree**; GRLP-standalone keeps its own `Network`
  as-is and shipped.
- The **GRLP process inside FluvTree is a thin adapter** that reuses GRLP's solver
  numerics on FluvTree's substrate.
- Re-homing GRLP-standalone onto the shared substrate is a **future** option, gated
  on FluvTree stabilizing — not a now-decision. One implementation is the goal;
  not coupling the library to the prototype is the constraint on getting there.

### Benchmark back to stable GRLP

The GRLP adapter must reproduce standalone GRLP. **Stable GRLP here (`grlp@366fb3e`,
v2.1.0) is ground truth:** evolve the same network from identical initial conditions
with both, and assert the elevation profiles match to tolerance. This is the
regression that lets the adapter refactor proceed safely.

## Open items

- **State schema.** Pin the exact field names and where each lives (edge vs node),
  and the reach-assembly convention, as a small written contract all processes obey.
- **GRLP adapter.** The focused refactor: `Network` reads/writes graph state instead
  of owning it. Numerics unchanged; ownership moved.
- **Scheduler v1.** Ordered list + `step(G, dt)`; `reads`/`writes` declared but
  validation deferred.
- **Units reconciliation** (TerraPIN areas ↔ GRLP `ssd`) — still open from the
  coupling note.
- **Whose `x,y`?** Node coordinates duplicated onto edge endpoints in
  `setupDomain.py`; under (A) they are canonical on the node — confirm no process
  needs a private copy.
