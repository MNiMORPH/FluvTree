# Enhancement (draft): general mechanistic bedrock-erosion module (abrasion + plucking)

Status: **parked enhancement** — to be filed as a GitHub issue (label: enhancement).
Scope: a physically-based bedrock-erosion module for `fluvtree.solvers.advection` / the process
layer, as a mechanistic alternative to the phenomenological stream-power law
(`StreamPower`), combining **abrasion** and **plucking**.

## Motivation

`StreamPower` (`dz/dt = U - K Q^m S^n`) is phenomenological: `K` lumps all
the mechanics. A mechanistic module resolves the actual incision processes, and --
importantly -- it is the **physically-grounded form of the alluvial<->bedrock
transition** we already want (see the sand `S->0` issue's taper idea and
`FixedBedProcess`). The **cover effect** *is* that transition, done from physics
instead of by a hand-tuned spline:

- sediment shields bedrock, so erosion depends on the **cover fraction ~ Q_s/Q_t**
  (a smooth, bounded, monotone transition variable -- the "spline" for free);
- plucking/abrasion convert **bedrock -> sediment**, which *is* the mass
  conservation the transition needs (rock removed becomes load, tracked, not
  fabricated).

So this enhancement **subsumes** the alluvial<->bedrock taper: `FixedBedProcess`
(fully covered / infinitely strong) and detachment-limited incision (starved bed)
both fall out as limits of one law.

## Physics to implement

- **Abrasion** (Sklar & Dietrich saltation-abrasion): erosion by impacting bedload
  "tools", reduced by bed "cover".
  `E_abrasion ~ (tools) x (1 - cover) / (rock resistance)`, with
  `tools = Q_s` (sediment supply), `cover ~ Q_s / Q_t` (capacity fraction),
  `rock resistance ~ sigma_T^2` (tensile strength). Produces the **humped**
  tools-and-cover response (max erosion at intermediate supply). Abrasion yields
  fines (wash load).
- **Plucking**: threshold/hydraulic removal of jointed rock blocks (lift/drag +
  impacts exceeding a block-stability threshold set by joint spacing / block size).
  Dominant in fractured rock; the main **sediment producer**. Formulation is less
  unified than abrasion -- threshold- and block-size-based; couple to
  weathering/fracture state if desired.
- (Optional) **macroabrasion**: fragmentation of the bed into pluckable blocks --
  Chatanantavet & Parker combine all three.

## Architecture / fit in FluvTree

This is inherently a **coupled** process (bedrock erosion <- sediment flux <-
alluvial transport, with sediment fed back) -- the first genuine multi-process
coupling, and exactly what the shared-graph + scheduler exist for.

- **New graph state: two layers** -- bedrock-surface elevation and an alluvial
  sediment-mantle thickness per reach (the SPACE / Shobe et al. structure; the bed
  is bedrock where the mantle is thin/absent, alluvial where it is thick).
- **Reuse `fluvtree.solvers.diffusion`** for `Q_s` and transport capacity `Q_t` (already in the
  closures) -- no new implicit engine; the erosion rate is **local and algebraic**
  (an explicit source/sink on the bedrock surface, closer in shape to `ExplicitRule`
  than to a PDE solver).
- **Process shape**: sits in the scheduler between the alluvial transport (produces
  `Q_s`, `Q_t`) and the sediment budget; reads capacity + cover, lowers bedrock,
  updates the mantle, and hands plucked/abraded rock to `Q_s`. Mass-conserving.
- **Limits (acceptance)**: `-> FixedBedProcess` (covered / strong / no tools) and
  `-> detachment-limited` (starved) recover as end members of the same law.

## Staged implementation (each step testable)

1. **Two-layer state** on the graph: bedrock surface + sediment-mantle thickness
   (fields + accessors). Entrainment/deposition bookkeeping.
2. **Abrasion "tools and cover"** process (`E ~ Q_s (1 - Q_s/Q_t) / sigma^2`),
   reusing `fluvtree.solvers.diffusion` for `Q_s`/`Q_t`. Validate the **humped** supply response.
3. **Plucking** (threshold + block size) and the **mass hand-off** to `Q_s`.
4. **Limit checks**: recover `FixedBedProcess` and detachment-limited behavior as
   end members; conserve the bedrock+sediment budget.

## Open questions / choices

- Plucking law: threshold formulation, block size / joint spacing as inputs;
  whether to couple to a fracture/weathering state.
- Cover function: linear `(1 - Q_s/Q_t)` vs Turowski-style probabilistic /
  intermittent cover (Lague 2010).
- Fate of abrasion fines: leave the reach as wash load, or track as suspended load?
- Two-layer numerics: explicit entrainment/deposition vs a small implicit coupling
  to the alluvial Exner; stability of the coupled step.
- Rock-strength / abrasion-coefficient parameterization and units.

## Relationships

- **Subsumes** the alluvial<->bedrock smooth-taper + mass-conservation idea
  (cover effect = that transition, physical).
- **Builds on** `fluvtree.solvers.advection` + `StreamPower` (phenomenological sibling) and
  `FixedBedProcess` (the non-erodible limit / anchor).
- Shares the two-layer sediment-mantle machinery with any future
  transport/deposition tracking.

## References

- Sklar & Dietrich (2004), *WRR* -- saltation-abrasion.
- Whipple, Hancock & Anderson (2000), *GSA Bull.* -- plucking vs abrasion vs
  cavitation regimes.
- **Chatanantavet & Parker (2009), *JGR Earth Surf.*** -- mechanistic abrasion +
  plucking + macroabrasion combined (closest to this enhancement's scope).
- Shobe, Tucker & Barnhart (2017), *GMD* -- SPACE: coupled detachment-transport,
  two-layer architecture precedent (Landlab).
- Lague (2010), *JGR* / Turowski et al. -- cover intermittency / probabilistic cover.
