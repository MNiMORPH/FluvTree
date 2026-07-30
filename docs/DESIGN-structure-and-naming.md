# FluvTree — package structure & naming spec

Status: **locked (2026-07-30).** The reference the `src/`-layout rebuild is built
from and checked against. Names may still be refined later, but this is the agreed
target we build to now.

## Model & structure — the two nouns

- **`FluvTree` — the model.** Holds a `RiverNetwork` plus the ordered processes,
  the scheduler, and the clock; its methods modify the variables sitting on the
  structure. Embodies the guiding principle: *FluvTree holds what is; its methods
  modify what is.*
- **`RiverNetwork` — the structure.** The directed convergent graph and the
  variables on it (topology + state). Structural methods only (traversal, field
  get/set); no physics.

```python
import fluvtree as ft
net   = ft.RiverNetwork.from_…(…)     # the structure + its variables
model = ft.FluvTree(net)              # the model: holds structure, owns the methods
model.add(ft.GravelBed(D=0.05))       # attach a process that modifies the variables
model.run(until=…, dt=…)
```

## Scope

Convergent river networks only — a deliberate invariant. Distributary, cyclic
(braided), dynamic-topology, and step-backwater generalizations are analysed and
deferred (see the frontiers note / memory).

## Layout — single `fluvtree` namespace under `src/`

```
src/fluvtree/
  __init__.py          front door — import fluvtree as ft
  network.py           RiverNetwork
  model.py             FluvTree
  scheduler.py         Scheduler (the machinery inside FluvTree.run)
  solvers/
    diffusion.py       DiffusionSolver   — power-law nonlinear diffusion
    advection.py       AdvectionSolver   — power-law nonlinear advection (n=1 rung built)
  closures/
    base.py            TransportClosure
    gravel.py          GravelClosure
    sand.py            SandClosure
  processes/           GravelBed, SandBed, StreamPower, FixedBed, Rule, GRLP
scripts/               GIS utilities (setupDomain, drainageNetworkGRASS)
tests/  docs/  pyproject.toml
```

## Three-layer taxonomy

- **closure — the physics.** A constitutive law that closes an under-determined
  governing equation. Constitutive-closure sense (a flux/material law, like an
  equation of state or a friction law) — *not* moment-hierarchy truncation.
  Stateless, numerics-free, graph-free.
- **solver — the numerics.** Discretization + time integration + solve for a
  governing *form* on the network (assembly, the conservative junction cell,
  implicit matrix + Picard, or the explicit sweep). Law-agnostic; parameterized by
  a closure.
- **process — the operator.** A schedulable unit that composes a solver + a
  closure onto the `RiverNetwork`, declares `reads`/`writes`, and exposes
  `step(dt)`. Named by the specific physics it embodies.

## Naming rules (locked)

1. **Solvers are named by their math**, not their domain — `diffusion`,
   `advection`. Symmetry: flux `Q_s ∝ Q·S^p` → diffusivity `∝ S^(p-1)`; erosion
   `E ∝ Q^m·S^n` → celerity `∝ S^(n-1)`. Domain words (alluvial, bedrock) live in
   the process layer and docs — never as a solver or process class name.
2. **Closures are per-solver-form** — each solver form carries a closure family:
   transport closures → diffusion; incision (stream-power) closures → advection;
   friction closures → the (future) flow solver. There is **no single "bedrock
   closure"**: a genuinely different form (saltation–abrasion) is a *different
   solver* with its own family.
3. **Processes are named by the specific law/mechanism**, not the domain family:
   `GravelBed`, `SandBed` (specific transport closures), `StreamPower` (a specific
   incision law), `FixedBed` (a specific null behavior), `Rule` (an arbitrary
   supplied `dz/dt`). Future siblings — `SaltationAbrasion`, `Plucking` — slot in
   without renaming anything. "Bedrock" would drift as those arrive; the specific
   law does not.
4. **Register: name for use, document for precision.** Module paths and class
   names stay clean (`solvers.diffusion`, `DiffusionSolver`); the full precision
   ("power-law nonlinear diffusion; flux ∝ Q·S^p") lives in the first docstring
   line. Promote precision into the identifier only when a sibling forces
   disambiguation (e.g. a future *linear* hillslope diffusion).

## Locked names

| role | name |
|---|---|
| package / alias | `fluvtree` / `ft` |
| model | `FluvTree` |
| structure | `RiverNetwork` |
| solvers | `DiffusionSolver`, `AdvectionSolver` |
| closure interface | `TransportClosure` |
| transport closures | `GravelClosure`, `SandClosure` |
| processes | `GravelBed`, `SandBed`, `StreamPower`, `FixedBed`, `Rule` |
| external cross-check | `GRLP` — adapter to published `grlp`, for validation |
| scheduler | `Scheduler` |

Not yet locked: a generic bring-your-own-closure process (would be named by its
*solver form*, per rule 3), and incision closures becoming first-class once a
second bedrock law lands.

## Realized by

The born-in-`src` history rewrite: replay the substrate commits under
`src/fluvtree/` (content, messages, author, dates preserved — paths/imports only),
add the template diffusion solver + transport closures + the advection solver,
unify the front door, move GIS scripts to `scripts/`, switch `pyproject` to
src-layout. Substrate replay authorized.
