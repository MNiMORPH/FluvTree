# Combined river + lake network: the input spec

Status: **proposal (2026-09-12).** Not built. This specifies what a *producer*
package must hand FluvTree so that rivers and lakes live on one network.

## The scope decision

River + lake **topology and geometry** are produced **outside** FluvTree and
delivered as an input. FluvTree supplies water and sediment and evolves `z`.

The producer may *use* `r.richdem`, but is not necessarily part of it; it could be
its own GRASS package. Section 4 gives the evidence for that separation.

The contract in one line:

| | owns |
|---|---|
| **producer** | what the valleys and basins *are* – the container, as of t=0 |
| **FluvTree** | what water and sediment *do* in them – the evolving state |

## 1. The unifying idea: one geometry object, not two

A lake is not a different kind of object from a valley reach. Both are: an
along-flow coordinate `x`, a floor elevation, and a **width that varies with
elevation**, `B(x, z)`.

FluvTree already has exactly this primitive. `ValleyGeometry` (`src/fluvtree/valley.py`)
stores a ragged, piecewise-linear `B(x, z)` and gives `V = ∫B dz` with `dV/dz = B`
holding exactly. Read it as a lake and it *is* the hypsometry that Fill–Spill–Merge
needs; read it as a valley and it is the sediment-storage geometry the diffusion
solver already uses.

So the producer has **one** geometry deliverable, not a river one and a lake one.

**Consequence for the spec, and it is load-bearing:** do not clip the `(z, B)`
levels at the present water surface or the present valley floor. A segment that is
a river today may be drowned tomorrow, and a lake that fills becomes a river. Each
node must carry enough elevation range that the *other* classification is
representable without going back to the DEM.

### 1.1 Depression links: hypsometry is authoritative, `B(x, z)` is derived

Laying `B(x, z)` along a depression's deepest axis is only faithful if that axis is
straight and the basin is prismatic about it. Real depressions are lobate, bent, or
dendritic (a drowned valley network), so a perpendicular-width sweep along the
thalweg will **not** integrate to the true storage.

Resolve it by inverting which quantity is primary. The depression's area–elevation
curve `A(z)` is a genuinely 2-D measurement and comes from the DEM; the link widths
are then whatever reproduces it. For a depression spanning links *l* with nodes *i*:

> **Σ B(x_i, z) · Δx_i  =  A(z)   at every z**

`B` on a depression link is therefore an **effective** width, not a measurable
perpendicular one. Volume is then exact by construction and independent of how bent
the axis is:

> V(z) = ∫ Σ B Δx dz' = ∫ A(z') dz' = the true depression volume.

This costs nothing in the solver, which already reads `B` per node and whose
`storage_volume` is documented as *"per unit valley length"* (`ValleyGeometry.storage_volume`), so
the `B · Δx` product is the quantity it already forms.

What the along-axis *distribution* of `A(z)` controls is only **where inside the
depression** material accumulates – the delta-progradation question – not **how
much** the depression holds. That is the right place for a modelling choice to
live, and it is FluvTree's, not the producer's.

**Anchoring the axis.** The producer does not have to invent it: `dephier` already
emits `pit_cell` (the deepest point) and `out_cell` (the sill). A depression link
runs from where the network enters the depression, through `pit_cell`, to
`out_cell`.

### 1.2 Meta-depressions are groupings of links, not links

`A(z)` differs through a depression, and a meta-depression's extent is its two
children plus the ground between them. That tempts one link per hierarchy node –
but a meta-depression link would be spatially disjoint (a piece on each side), and
a disjoint link is no longer the contiguous `B(x, z)` object a river link is. The
unification would break at exactly the point it is meant to hold.

**Resolve it by not putting meta-depressions on the graph at all.**

- **Links are cut at leaf granularity**: one contiguous depression link per leaf
  depression, plus the ordinary river links that connect them.
- **The hierarchy is a tree whose leaves reference links**, carried as data
  alongside the graph. A meta-depression is then a *set* of links – "above
  `out_elev`, links {L3, L5, L7} share one water surface" – and a set has no
  obligation to be a single link.

The gap between two merge partners is not a hole: it is the reach over the pour
point. `dephier` records that connection explicitly – `geolink` is *"the depression
on the other side of the pour point"* and `odep` *"the depression this one
overflows into"*. So a meta-depression's links form a **connected subgraph**, not
disconnected pieces.

**This requires one property of the network, and it is load-bearing: the network
must be COMPLETE.** Continuous from every channel head, through every depression,
out over its sill, to the outlet – no interior termini. Filling, breaching, and
least-cost-path routing all deliver that; **the spec requires the property, not the
method**. The depression hierarchy is computed on the *original* DEM either way.
Without completeness, flow terminates at every pit, the connecting reaches do not
exist as links, and the depression links really are disconnected – this whole
construction is then unavailable.

**The storage constraint becomes hierarchical, and is free.** With §1.1 applied per
link, for any node *D* of the tree:

> Σ over links in *D*'s subtree of B(x_i, z) · Δx_i  =  A_D(z)

Satisfying it at the leaves and connecting reaches satisfies it at every level
above, **provided the extents partition** – every cell of a meta-depression belongs
to exactly one leaf or one connecting region. That partition requirement is on the
producer, and is the one thing that must be checked rather than assumed.

### 1.3 The maximum water level is computed, not guessed

§1 asks for `(z, B)` levels over a "generous" elevation range. The hierarchy makes
that precise. Taking a **depression** to mean the top-level non-ocean-linked
meta-depression, that node's `out_elev` is the highest water level the system can
reach before spilling to the outlet. It is the exact upper bound for every `B(z)`
table beneath it, and it bounds how much of the tree needs shipping: the sub-tree
under each top-level non-ocean-linked node, and no further.

### 1.4 Junction degree is a function of water level

Many tributaries entering one lake is not an exception to the binary graph; it is
the hierarchy read backwards. At high water the tributaries meet one water surface
– an *n*-way junction. As the level drops the water body splits into sub-basins and
the tributaries re-separate into distinct binary confluences. The depression
hierarchy is precisely the record of which merges happen at which elevation.

So: **store the finest (low-water) decomposition, and merge at run time.** The
stored graph stays binary; FluvTree groups links into shared-surface bodies each
step, and the hierarchy is the precomputed answer for which groups form at which
level – a lookup, not a discovery. Tracing the tributary network to partition the
hierarchy in advance is then unnecessary.

The solver already permits the run-time result: confluence assembly sums over
`for j in up_ids` with no cap (`DiffusionSolver._assemble`), so an *n*-way shared surface
costs nothing.

### 1.5 What each link carries, and what FluvTree writes

Junction degree being state needs one stored field, not a list.

**One label per LINK: the smallest hierarchy node whose OWN region contains it.**
Every degree of freedom then inherits it unambiguously, because a graph node
belongs to its **downstream** link – `('jcn', s)` is named for the segment *leaving*
it (`network.py:from_segment_lists`; *"a confluence node is named by the single
segment leaving it"*) – and in the solver's layout each segment carries all its own
nodes. So no DOF is ever shared between two labels. See §1.6 for the one exception.

A *hierarchy* node's own region is its extent minus its children's – so a leaf's
own region is the basin itself, and a meta-depression's own region is just the
saddle ground between its two children. Those own-regions tile the landscape exactly (this *is*
the §1.2 partition), so every link falls in exactly one:

| link | `dep_label` |
|---|---|
| inside a leaf depression | that leaf |
| a connecting reach over a pour point | the **meta**-depression that drowns it |
| outside every depression | null / ocean-linked |

Downstream ownership also gives the right answer at both ends of a depression with
no special rule: at the **entry** margin the node's downstream link is the
depression, so it joins the depression; at the **sill** its downstream link is the
river, so it does not – which is correct, since nothing above `out_elev` is stored
in that depression anyway.

So tributary links do **not** carry duplicate labels of the same depression. They
carry *different* labels, and what they share at high water is a **common ancestor
in the tree**. "Which links are in meta-depression *M*?" is then a subtree query,
and membership is derived, never stored twice.

**The run-time state.** The filled bodies at any instant are an **antichain** of
the tree – a set of hierarchy nodes, none an ancestor of another. Each step:

1. take the active antichain and each active hierarchy node's water level (FSM's
   output);
2. for each active hierarchy node *D* at level `z_D`: every DOF whose (inherited) label lies
   in *D*'s subtree, and whose bed is below `z_D`, joins one body with surface `z_D`;
3. everything else is river.

A graph junction whose adjacent DOFs are all in the same body is not a hydraulic
confluence at all – it is an interior point of a flat water surface.
That is where the *n*-way behaviour comes from: the body groups the links, and the
graph never changes.

**Worked case.** Leaves *A*, *B*; parent *M*; saddle reach `L_S` between them;
tributaries `L1, L2 → A` and `L3 → B`.

- *Low water*, antichain `{A, B}`: body *A* = `{L_A}`. `L1` and `L2` meet it as two
  ordinary binary confluences. `L_S` is a river.
- *High water*, `A` spills over its `geolink` into `B` and they merge; antichain
  `{M}`: body *M* = subtree(*M*) = `{L_A, L_S, L_B}`. Now `L1, L2, L3` all end on
  one surface – three tributaries, one water body – while the stored graph is still
  the same binary graph, and the junctions inside the lake are interior points.

**The division of labour this implies**, and it is the guiding principle applied
unchanged: the producer ships what *can* happen – the hierarchy is every possible
lake configuration – and FluvTree writes what *is* happening. Stored and immutable:
`dep_label`, the tree, merge elevations, `A(z)` / `B(x, z)`. Written each step by
FluvTree: `seg_type`, `z_w`, and the body each DOF currently belongs to.

### 1.6 Short segments: cut by default, with a per-node fallback

Depression margins and tributary junctions are both fixed by terrain. When they
land close together, cutting the graph at the margin leaves a stub segment.

**Measured behaviour first** (steady state, gravel closure, uniform `Q`, analytic
fixed point `S → S0`; scripts in the session scratchpad):

| squeezed feature | length | `max|S-S0|/S0` | `cond(LHS)` |
|---|---|---|---|
| 1:1 segment in series | 2000 m | 3.7e-15 | 3.2e+01 |
| | 20 m | 9.1e-13 | 2.9e+05 |
| | 0.2 m | 4.0e-11 | 2.9e+09 |
| 2:1 tributary-to-confluence face | 500 m | 1.1e-14 | 2.9e+01 |
| | 1 m | 6.3e-12 | 1.5e+04 |
| | 0.01 m | 1.0e-10 | 1.5e+06 |

Two things follow. **(1)** There is no silent-wrong-answer regime here: the steady
state stays correct to ~1e-10 across five orders of magnitude of shortening, and
the confluence land area tends to a finite limit rather than vanishing.
**(2)** Conditioning degrades as `1/L` and nothing else does, leaving large headroom
in double precision. Node-count violations remain a clean failure
(`ValueError: confluence segment 2 needs >= 3 nodes (has 2)`, verified).

So the hazard is real but mild, and does not justify inventing a minimum-length
threshold. **Cut at margins by default and keep the label on links**; where a cut
would leave a stub next to a junction, decline it and carry a per-node label array
on the one segment that spans the margin. The network holds arbitrary per-segment
arrays, so the fallback costs nothing structurally.

**Assign by cell overlap, not by point-in-polygon.** A link overlies a depression
over a measurable *length*, which makes link-to-depression a principled
extent-against-extent comparison. A bare point has no such measure, so labelling
"is this node inside the polygon?" would be a knife-edge and would be a genuine
step down from link assignment. (This is why §1.5 keeps the label on links; the
paragraph below records why per-DOF assignment is nonetheless well posed where the
fallback needs it.)

A DOF is not a bare point here. It owns an explicit control volume: half the
upstream spacing plus half the downstream spacing, which the solver computes
outright – `land = B[1] * 0.5 * ((x[1] - x[0]) + (x[2] - x[1]))`
(`DiffusionSolver._assemble`, and the same form in `_compute_land_areas`). So the assignment rule is

> intersect each node's **cell** with the depression own-regions, exactly as one
> would intersect a link, and assign the cell by its overlap.

That restores the extent-against-extent comparison, and at **finer** resolution
than links: a link is a union of cells, so link-level assignment is precisely the
*coarsening* of cell-level assignment. Cell-level therefore dominates link-level
rather than trading against it.

The stronger reason to prefer it: the partition should be over **the same control
volumes the solver conserves over**. Then `Σ B(x_i, z) · Δx_i = A(z)` is a
statement about the very volumes that appear in the mass balance, not about a
different tiling that happens to cover the same ground.

**Cutting at the margin makes the partition exact – bracket the margin, never land
on it.** Cutting is topologically cleaner, and it does more than that: it removes
the straddling cell entirely, because the two segments' adjacent nodes are
separated by a **face**, and a face is exactly where control volumes meet. The
solver duplicates the junction value rather than sharing it – *"the confluence
value is duplicated across the mainstem's first node and each tributary's last
node, separated by a finite face; junction nodes on the graph are topological
markers ... not shared `z` degrees of freedom"* (`solvers/diffusion` module docstring) – and the
face sits at the midpoint of the gap (`land = B[-1] * 0.5 * ((x[-1] - x[-2]) +
(x[dj][0] - x[-1]))`, in `DiffusionSolver._assemble`).

So a margin at `x_m` is represented by **two nodes bracketing it**, at
`x_m - δ/2` and `x_m + δ/2`: the upstream cell is then wholly outside the
depression and the downstream cell wholly inside. No majority rule, no fraction,
no quantisation.

**The trap.** Placing a *single* node at `x_m` – the natural reading of "cut here"
– gives the two duplicated DOFs the same coordinate, hence zero face length. That
is a division by zero in `_face_conductance` and an exactly singular matrix, not a
degraded answer. It is easy to do by accident; it arose twice while building the
§1.6 experiments. **The producer must guarantee strictly increasing `x` across
every cut.**

**Where a fraction is still needed.** Only where the producer *declines* to cut –
near a junction, to avoid a stub. There the margin falls mid-cell and the cell is
split between two labels, 50/50 in the symmetric case. That is straightforward to
carry as a `(label, fraction)` pair per side, and it is the right thing to do
there. Note it is never *required* for volume: because `B` is derived from `A`
(§1.1), the producer absorbs any quantisation into the widths and volume stays
exact either way. The fraction buys sub-cell basin membership, which matters at
the outer boundary of the top-level meta-depression, where "inside vs draining to
the outlet" is the distinction being resolved.

**A correction to an earlier draft of this section.** It argued that per-link
labelling "misplaces the shoreline by up to a link length". That was wrong. The
shoreline is set by the `z_bed < z_w` test in §1.5, evaluated per node, and is
resolved at node spacing whatever carries the label. Label granularity governs
the *storage partition* and *basin membership near margins* – not the shoreline.

What this does **not** fix: two tributary junctions falling close together. That is
a property of the drainage network itself, present with or without depressions, and
is left open.

**What the measurements do not cover, stated so it is not mistaken for settled.**
They are **steady-state only**, on a gravel closure with uniform `Q` and a
rectangular valley. Untested: the transient approach (a tiny cell has a short
response time, and the Picard iteration under a dynamic source may behave
differently); a *lake* link, whose `B` is orders of magnitude larger than its
neighbours' and which was never in the test network at all; and the advection
solver. The conclusion "short segments degrade conditioning but not accuracy" is
therefore warranted for the case measured and is a hypothesis elsewhere.

### 1.7 Making the geometry updatable without a callback

The structure must support geometry that changes as the landscape evolves. Four
cases, and three of them stay inside FluvTree:

1. **Sediment fills a basin.** The floor rises; `A(z)` is the original curve with
   its lower part truncated. Raising the datum is internal.
2. **A sill erodes.** `out_elev` drops, truncating the top of the usable curve.
   Internal.
3. **Two depressions merge.** The combined curve is the *meta*-depression's `A(z)`
   – already a node in the hierarchy. **This is free if and only if the producer
   ships hypsometry for meta-depressions, not just leaves.** See §2.5.
4. **A new depression forms** where the network was a river. Internal *provided*
   the initial delivery carried `(z, B)` levels over a generous elevation range
   (§1). Ponding is then just integrating a curve that is already there.

The limit, stated honestly: if the planform changes enough that the **axis itself
is wrong** – a lake growing laterally into a sub-basin the graph does not
represent – no amount of stored `B(z)` recovers it, and the producer must be re-run.
That is a real boundary of the one-shot-file design, not something to paper over.

## 2. Deliverables

**Two tiers, and the distinction is the organising one for this whole file.**

- **Core – the topology.** What makes this a river-and-lake network and what nothing
  else can supply: the graph, where its nodes sit, the depression hierarchy, and the
  storage geometry that makes a depression fillable rather than merely labelled.
  §2.1, §2.2, §2.3, §2.5.
- **Convenience – whatever else is cheap to pull from the DEM.** Not part of the
  network; attached because the producer is holding the raster anyway and the
  alternative is recomputing it later from a DEM the modeller may no longer have.
  §2.4.

§2.6 is empty by the same test (boundary conditions are modelling choices), and
§2.7–2.8 stand outside both tiers: provenance is not payload, it is what makes the
payload trustworthy, and is required regardless of what else is present.

Constraints below marked *(verified)* are read off the current code, with the line
that enforces them.

### 2.1 Topology, per segment — *core*

1. `seg_id` – integer, stable.
2. `upstream_segment_IDs`, `downstream_segment_IDs` – the form
   `RiverNetwork.from_segment_lists` takes.
3. **At most one downstream segment.** *(verified: `RiverNetwork.downstream_segment` raises
   "non-convergent topology".)* The joined graph must still be a convergent tree.
   This is the locked FluvTree invariant (`DESIGN-structure-and-naming.md`, Scope),
   and lakes must not break it.
4. **One outlet.** *(verified: `RiverNetwork.from_arrays`, `DiffusionSolver.__init__`, and
   `advection.sweep_order` / `incise_n1_step` all take `mouth_segments()[0]`.)* Multiple outlets are not
   supported today; see §7.
5. `dep_label` – **per link**, the hierarchy node whose **own region** contains it,
   per §1.5; null outside every depression. Per node only on an uncut margin
   (§1.6), which the network holds as an ordinary per-segment array. This is the immutable field the
   producer ships. `seg_type` (`"river"` / `"lake"`) is **not** shipped as
   structure – FluvTree derives it each step from the active antichain – though an
   initial classification is useful as a t=0 state and as a check on the producer.

### 2.2 Along-flow geometry, per segment, per node — *core*

6. `x` – along-flow distance [m], monotonic.
7. `z_bed` – bed or basin-floor elevation [m].
8. **Node-count minima.** *(verified: `DiffusionSolver._check_confluence_sizes`.)* A segment
   downstream of a confluence needs ≥ 3 nodes; a tributary entering a confluence
   needs ≥ 2. The producer must honour this when it chooses breakpoints, or
   FluvTree raises at solver construction.

**Which nodes go in the segment array is per-process, and the producer must
declare which convention it emits.** This is settled FluvTree design, not an open
question – `c697e468`, "Correct reach-assembly convention: GRLP has no shared
junction dof":

> The shared-junction-node interpretation is **per-process** (natural for
> one-dof-per-node processes like TerraPIN), **not universal**. Verified against
> `grlp@366fb3e` that GRLP owns all nodes per segment and couples confluences via a
> matrix row, not a shared `z`. **The substrate stays agnostic and stores whatever
> fields each process needs.**

So the two readings are both legitimate and neither is *the* canonical shape:

- **all-nodes-per-segment** – what the alluvial diffusion solver consumes, where
  the junction is coupled by a matrix row across a finite face and there is no
  shared junction `z`;
- **interior-plus-shared-endpoints** – natural for a one-DOF-per-node process,
  which is what `reach_profile` serves.

`RiverNetwork` defines the topology and holds whatever fields a process needs; a
solver *consumes* that and builds its own state layout. Neither solver's internal
arrangement defines the interchange.

**The consequence for this file is a mandatory declaration, not a fixed shape.**
Emitting one convention while the consuming process expects the other is a silent
failure – short arrays read as complete, a subtly wrong profile, no error and no
warning. So the file carries `node_convention` (§2.8) and the reader **validates it
against the process being run** rather than assuming. Which convention a given
producer run should emit follows from the process it is feeding, and is therefore a
modelling choice, not a property of the terrain.

### 2.3 Cross-section geometry, per node, ragged — *core*

9. `z_levels`, `B_levels` – ascending `z`, paired widths. **≥ 2 levels per node**
   *(verified: `ValleyGeometry.__init__`.)*
10. `channel_width` – the backup width below the lowest tabulated level.

Store these as flat `z`, flat `B`, plus an `offsets` index. That is exactly
`ValleyGeometry`'s own in-memory layout, so the reader is close to a copy.

**The same table serves water and sediment, at two different porosities.** A
depression fills with water *or* with sediment through the identical `B(x, z)`
curve; only the solid fraction differs. `ValleyGeometry` bakes `(1 - lambda_p)`
into both `storage_volume` and `storage_jacobian` (`ValleyGeometry.storage_jacobian` / `.storage_volume`), so a
water-storage query needs the `lambda_p = 0` path. The producer therefore ships
**one geometry and no porosity** – `lambda_p` is FluvTree's sediment parameter,
not a property of the container, and must not be folded into the delivered
widths.

### 2.4 Convenience attributes — *not core*

11. `A` – contributing drainage area [m²]. **Optional, and recommended.**

Contributing area is a property *of* the network, not part of it: it accumulates
through the topology, and it is an input a particular process needs rather than
something that makes this a river-and-lake network. By the same test that empties
§2.6, it is not a required deliverable.

It differs from the §2.6 items in one way worth stating. `S0`, `z_bl` and the rest
have **no terrain basis** – the producer could not supply them if it tried.
Contributing area does, and the producer is uniquely placed to supply it: flow
accumulation is a **byproduct of the routing this spec already requires**, so
shipping it costs nothing, whereas recovering it later means going back to the DEM.
A file without it cannot feed the standard solvers unaided.

Hence *optional and recommended*, with one rule that generalises to every
attribute in this section: **the file declares which optional attributes it
carries**, as an explicit `optional_attributes` list in the provenance block. A
reader branches on that list rather than probing for variables, and an empty list
says "none, deliberately" rather than leaving absence ambiguous – the same reason
`min_depression_filter` records the literal `"none"`.

The criterion for this section is membership, not cost: **an attribute is optional
because it is not part of the network graph**, whatever it happened to cost to
compute. Cheapness is why these particular ones are worth attaching at creation
time; it is not what makes them optional.

**If supplied, supply area, not discharge.** FluvTree currently treats `Q` as fixed
for the run *(verified: `DiffusionSolver.__init__`, "Constant-in-time reach state")*, and
lakes break that – an endorheic basin passes no water downstream until it spills.
Area is time-invariant and routing-independent; discharge is neither. Give the
invariant quantity and let `Q` be derived as routing changes.

**What else is cheap from the same rasters.** These are proposals, not decisions –
each is nearly free given what the producer already computes, and each is dropped
without argument if unwanted:

12. **Sinuosity per link** – the planform length of the polyline over the
    straight-line distance between its ends. `DiffusionSolver` already takes
    `sinuosity` as a physics parameter and currently **defaults it to 1**
    (`DiffusionSolver.__init__`), where it enters `C0` as `sinuosity**-p`. The producer holds
    the actual polyline and can measure it. This is the strongest candidate: a real
    parameter FluvTree is presently guessing.
13. **Shreve / Strahler order per link** – pure topology, free once the graph exists,
    and the vocabulary `fluvtree.generate` already uses.
14. **Depression statistics** – `cell_count`, `total_elevation`. `dephier` emits them
    alongside `dep_vol` (already core, §2.5 item 16) at no extra cost.
15. **Local relief or hillslope gradient beside each link** – from
    `r.richdem.terrainattribute`. Not needed by any current process; the cheapest
    thing to have already measured if hillslope coupling ever lands.

Deliberately **not** on this list: valley width `B(x, z)` for *river* links. It is
core (§2.3), not convenience, and it is not cheap — it needs cross-sections cut
perpendicular to the channel, or TerraPIN. That is Gap 4 (§4.2), and calling it
convenient would hide the one piece of core geometry nothing currently produces.

### 2.5 The depression hierarchy — *core*

Per depression, including ones that are dry today:

12. `dep_label`
13. `parent`, `lchild`, `rchild` – the merge tree
14. `pit_elev` – floor; `out_elev` – **sill (spill) elevation**
15. `out_seg`, `out_node` – *where on the segment graph the sill sits*
16. `dep_vol` – capacity to the sill
17. endorheic / ocean-linked flag
18. `A(z)` – the area–elevation curve, for **every node of the tree**, leaf and
    meta alike. It is the authoritative geometry (§1.1); the link widths are
    derived from it. Shipping the meta curves makes a merge a lookup, and gives
    FluvTree an independent check that the §1.2 partition holds.
19. `odep` / `geolink` – the overflow connections, and the merge elevation per
    tree node. These are what turn the tree into a run-time lookup (§1.4).

Item 15 is the one that does not exist in any current output and is the crux of
the join (§5): the hierarchy knows the sill as a raster cell; FluvTree needs it as
a position on a reach.

### 2.6 Boundaries – deliberately empty

**No boundary condition is the producer's.** They sit on the topology; they are not
part of it, and none of them is a terrain product:

- `S0` / `Q_s_0` – upstream sediment supply has no terrain basis whatever. A DEM
  cannot tell you what it is.
- `z_bl` – whether base level is held, dropped or cycled is a modelling choice. The
  producer could offer the outlet's terrain elevation as a starting value, but that
  is a convenience, not a deliverable.
- `x_bl` – where the outlet ghost node sits. This looks geometric, but it is part of
  *how base level is imposed*, not of the terrain: it says where the condition is
  applied. The solver already defaults it by linear extrapolation
  (`DiffusionSolver._assemble`: `2 * x[-1] - x[-2]` when absent).
- **Initial water state** – how much water each depression holds at `t = 0`. It is
  "what *is* happening", which §1.5 assigns to FluvTree, and deriving it needs a
  water input that is itself a modelling choice. Recorded in §7.

The producer's whole obligation here is to make the boundaries **addressable**: say
which node is the outlet and which are channel heads. Both are topology, and both
are already covered by §2.1.

### 2.7 Georeferencing and provenance

20. `x`, `y` per node, and the CRS.
21. **The provenance block** – specified in full in §2.8.

### 2.8 The provenance block

Both routing engines are deterministic (§5.6), so the network is **reproducible**:
the same DEM through the same engine with the same settings yields the same graph.
That makes a real reproducibility contract possible rather than aspirational – but
only if the file carries everything that determines the answer. The organising
rule, as elsewhere: **the file states what made it, and the producer refuses to
write a file with any of it missing** – the same discipline as `trust/provenance.py`
raising when a step is skipped.

The test to apply to every candidate field: *could changing this change the graph?*
If yes it is mandatory. Timestamps and hostnames fail that test and are recorded
as context, never as identity.

**(a) Identity of this file**

    spec_version, created (ISO-8601 UTC), output_sha256

`output_sha256` covers the data payload excluding this block, so a re-run can be
proved identical rather than assumed so.

**(b) The input grid** – this alone fixes the discretisation

    dem_source, dem_sha256, dem_crs, dem_nodata
    dem_type           # CELL / FCELL / DCELL
    dem_units, dem_vertical_datum
    region: n, s, e, w, nsres, ewres, rows, cols

`dem_type` earns its place: **integer elevations produce vastly more equal-elevation
ties** than float ones, so it predicts how much of the network rests on a tie rule
rather than on terrain. A hash, not just a path – a filename is not an identity.

**(c) What generated it**

    producer_name, producer_version, producer_git_commit
    routing_engine     # "richdem" | "r.watershed"
    richdem_version, grass_version, r_richdem_version

`routing_engine` is load-bearing, not bookkeeping: the two engines break ties in
**opposite** directions (§5.6), so this field is what tells a later reader whether
two networks are comparable at all.

**(d) Every choice that can change the graph**

    optional_attributes      # explicit list of §2.4 attributes present,
                             #   e.g. ["A", "sinuosity"]; [] means none, and
                             #   means it deliberately, not by omission
    node_convention          # "all_nodes_per_segment" | "interior_plus_endpoints"
                             #   — §2.2; validated against the consuming process
    conditioning_method      # fill | breach | least_cost  (§1.2 completeness)
    flow_direction           # must be SFD/D8
    beautify_flats           # must be false  (r.watershed -b)
    flat_resolution          # "richdem:Barnes2014" | "r.watershed:astar"
    tie_break                # "LIFO" (richdem) | "FIFO" (r.watershed)
    stream_threshold, accumulation_units
    min_depression_filter    # "none" — recorded explicitly, §5.3
    # (no snap_tolerance: none is used — §5.5)
    stub_policy, stub_threshold                       # §1.6
    commands                 # the literal invocations, in order, verbatim

`commands` is what turns the block from a description into a reproduction recipe.

**(e) Counts – where the conventions actually bit**

Parameters say what was *configured*; counts say what it *did*. Both are required.

    n_links, n_nodes, n_heads, n_outlets
    n_depressions_leaf, n_depressions_meta
    n_depressions_dropped_by_filter
    n_cells_in_flats
    n_ties_outlet, n_ties_pit, n_ties_flat_direction
    n_cuts_declined_stub, n_margins_fractional

The three `n_ties_*` counts are the ones §5.6 exists for: they measure how much of
the network rests on a tie-break convention instead of on the landscape. A run
where they are large is not *wrong* – it is reproducible – but it is a run whose
flat-country topology would change under the other engine, and the reader deserves
to know that before comparing it with anything.

`n_depressions_dropped_by_filter` follows the standing rule that a filter reports
what it kept **and** what it removed.

**Both formerly-open thresholds are now settled at "none"** (§5.3, §5.5), so the
producer applies no filter and no snapping. `min_depression_filter` is still
recorded, with the literal value `"none"` – an absent field would read as an
omission, and the point of the block is that a reader can tell a decision from a
gap.

## 3. What FluvTree must do with it

`RiverNetwork.from_arrays` takes `(upstream_segment_IDs, downstream_segment_IDs,
x, z, Q, B, S0, x_bl, z_bl, Q_s_0=None)`. The new reader lands in the planned
`io/` layer, and must **import without GRASS installed** – the same discipline
that keeps `import fluvtree` matplotlib-free.

## 4. Can `r.richdem` produce this?

Partly, and the shape of what is missing is the argument for a separate package.

### 4.1 What exists today

| module | gives | covers |
|---|---|---|
| `r.richdem.dephier` | depression-labels raster, flow-directions raster, hierarchy vector | most of §2.5 |
| `r.richdem.fsm` | equilibrium water surface from a water input; updates `water_vol` | the water solve |
| `flowaccumulation`, `filldepressions`, `breachdepressions`, `resolveflats`, `terrainattribute` | conditioning and routing primitives | inputs to the above |

The hierarchy vector's attribute table (`librichdem/vector.py`, `_SCHEMA_COLS`)
already carries: `type` (leaf/meta), `pit_cell`, `out_cell`, `parent`, `odep`,
`lchild`, `rchild`, `geolink`, `pit_elev`, `out_elev`, `ocean_parent`,
`cell_count`, `dep_vol`, `water_vol`, `total_elevation`.

So **`out_elev` is the sill elevation, `pit_elev` the floor, `parent`/`lchild`/
`rchild` the merge tree, `dep_vol` the capacity** – items 12, 13, 14, 16 and 17 of
§2.5 are already produced and already round-trip
(`depressions_to_grass` / `depressions_from_grass`).

### 4.2 Four gaps

**Gap 1 – hypsometry is computed but never exposed.** `dephier` gives `dep_vol`,
a single integral up to the sill, plus `cell_count`. It does not give `A(z)` or
`V(z)`. FSM *does* map volume to water-surface elevation "from the depression's
hypsometry" (`r.richdem.fsm.md`), but that curve stays inside the C++ and is not a
GRASS output. Two routes: **(a)** derive it outside RichDEM from the labels raster
plus the DEM – for cells of label *L*, sort elevations, cumulative cell area gives
`A(z)`, its cumulative sum gives `V(z)`; or **(b)** expose it from RichDEM.
Route (a) needs no C++ change and is the one to try first.

**Gap 2 – the river side is a prototype, not a module.** The DEM → reach path is
`r.watershed` + `r.stream.extract` + `v.stream.network`, prototyped in this repo's
own `scripts/drainageNetworkGRASS.py`. That file is a working note, not a
producer: it hard-codes catchment cats (985, 1258), a stream threshold of 3000 and
an accumulation cut of 1e7, and its `v.stream.profiler` step is a bare shell invocation
(`v.stream.profiler cat=1258 …`) that is not valid Python – flagged in the file
itself by the preceding comment, "Run in shell or rewrite". It needs rewriting.

**Gap 3 – nothing joins streams to depressions.** No module today merges a
depression hierarchy and a stream network into one segment graph. This is the
actual new work, and §5 is where its decisions sit.

**Gap 4 – valley width `B(x, z)` for *river* segments is not a depression product
at all.** It needs cross-sections extracted perpendicular to the channel, or
TerraPIN. Nothing produces it today.

### 4.3 Why that argues for a separate package

`r.richdem` describes itself as a GRASS add-on collection *wrapping RichDEM*.
Gaps 2, 3 and 4 are not RichDEM wrappers – they are stream-network extraction, a
graph join, and cross-section geometry. Putting them in `r.richdem` would make
that collection something other than what its README says it is. A separate
package that **depends on** `r.richdem` for the depression side keeps both honest.

## 5. The join: where the real decisions are

These are decisions, not defaults. None should be picked silently.

**Settled during design** – kept here so the list is not read as still open:

1. ~~**Breakpoints.**~~ **Settled: cut at margins and sills by default**, bracketing
   the margin with two nodes rather than landing one on it; decline the cut where
   it would leave a stub beside a junction and carry a per-node label array on the
   spanning segment instead (§1.6).
2. ~~**Leaf depression → one lake segment.**~~ **Settled**: links at leaf
   granularity, meta-depressions are groupings of links rather than links of their
   own (§1.2).

**Still open:**

3. ~~**Minimum depression size.**~~ **Settled: no filter. Every depression the
   hierarchy finds is kept**, single-cell pits included.

   The filter is not, after all, unavoidable – it is avoidable by moving the
   decision from the producer to the dynamics. A one-cell pit fills in the first
   step and stops being a lake, so the criterion becomes *physical* (does it
   actually fill?) rather than *geometric* (is it bigger than some number?). That
   is strictly better than a preprocessor guess, and it removes the last invented
   threshold from the producer.

   **The producer's goal is exactness; pruning is a downstream optimisation.** Ship
   the exact network the terrain implies, and prune later if computation demands
   it – a prune applied then is a recorded, reversible choice made against a known
   cost, whereas a filter applied here is an irreversible guess baked into the
   data. The two consequences – a hierarchy node per pit (a data-volume and lookup
   cost, not a correctness one) and noise pits briefly labelling channel links as
   depressions at t=0 – are both visible through the §2.8 counts, which is part of
   why those counts are mandatory.
4. ~~**Endorheic basins and the ocean.**~~ **Settled.**

   **Endorheic-ness is dynamics, not topology.** Because §1.2 requires a *complete*
   network, conditioning already supplies an outflow path over every rim, so an
   endorheic basin is not a second terminus – it is a depression whose level never
   reaches its sill, leaving its outflow link unused. That falls out of FSM and asks
   nothing of the producer. The ocean maps onto the outlet, and `ocean_parent` /
   ocean-linkage already ship (§2.5 item 17).

   **A lake cannot spill across a catchment divide.** The divide is *defined by* the
   routing, so a lake whose lowest escape leads into the neighbouring basin simply
   has the divide drawn the other way; the hierarchy already takes `out_elev` as the
   lowest rim point, so the two are consistent by construction. Either the sill lies
   below the divide and the lake drains within its own catchment, or the lake is
   endorheic with respect to it. Inter-catchment spilling is therefore not a
   configuration a consistently-computed network can hold, and **multi-outlet
   support is not forced by it**.

   Proglacial outlet switching is real, but it is driven by ice-margin and
   land-surface elevation change – *terrain* change, not water-level change. That
   means recomputing the catchment, i.e. re-running the producer.

   **Sea level: cut, do not recompute – when it is inundation.** A falling sea level
   exposes shelf and extends the network seaward; a rising one drowns the lower
   network. Both are handled by the machinery already built: compute the network
   once over the **maximal subaerial extent** (the lowest stand of interest, so the
   DEM must carry bathymetry to that depth) and treat sea level as **just another
   water level**. The ocean is then the outermost water body, `ocean_parent` becomes
   *derived* rather than static exactly as `seg_type` is, and a depression below the
   current stand is a bay while the same depression above it is a lake. No new
   mechanism.

   The honest limit, and it is the same one as above: this handles **inundation**,
   not **deformation**. Isostatic rebound, shelf deposition and the erosion of
   exposed shelf change the terrain itself, and terrain change means a re-run.

   **The general rule this yields:** *level changes are handled by the hierarchy;
   terrain changes require re-running the producer.*

5. ~~**Snapping tolerance.**~~ **Settled: none needed. Place the node at the
   sill.** The §1.6 measurements are what license this: a node may be inserted
   arbitrarily close to an existing one without corrupting the answer – conditioning
   degrades as `1/L` and nothing else does – so there is never a reason to snap a
   sill to a nearby breakpoint rather than put a node where the sill actually is.

   **What remains is a degeneracy rule, not a tolerance** – and it is worth being
   exact about how one would ever hit it, because it is not a freak case.

   **The mechanism is mistaking a face for a node.** In the solver's layout a
   confluence is *not* a node: it is the **face** between the tributary's last node
   and the mainstem's first node. The natural GIS reading is the opposite – a
   junction is a point on the map, so put a node there – and taking both sides to
   that point makes the face length zero. Measured:

   | mainstem's first node | result |
   |---|---|
   | at the confluence point (8000 m) | `Matrix is exactly singular`, **z all NaN** |
   | one spacing downstream (10000 m) | finite, correct |

   Note the failure mode: `evolve` **does not raise**. It warns and returns NaN. So
   unlike the short-segment case of §1.6, this one *is* a silent-wrong-answer path,
   and it is the reason the validator exists. The test suite's own convention keeps
   the separation – tributaries end at 8000 while the confluence segment starts at
   10000, a full spacing apart.

   **To be clear about what is topology and what is discretisation:** a confluence
   *is* a node of the graph, and so are a sill and a depression margin – the
   producer emits all three as nodes, and `from_segment_lists` creates `('jcn', s)`
   for every one. What is *not* a node is the solver's **sample point**: the two
   segments meeting at that graph node each carry their own terminal sample, and it
   is those two samples that must not share an `x`. The graph node sits between
   them. §1.6's "bracket the margin" is a statement about sample placement around a
   node, not a denial that the node exists.

   A second, duller mechanism: sampling node positions from vector-line vertices
   without de-duplicating repeated ones. The prototype does something adjacent –
   `drainageNetworkGRASS.py` concatenates every vertex of every line into one flat
   array and takes a cumulative distance *across* line boundaries, so its junction
   spacings are not meaningful to begin with.

   **What is explicitly fine: coincidence between siblings.** Two tributaries
   entering one confluence may share identical `x`, and the tests do exactly that
   (`[_D * np.arange(1, 5.0)] * 4`) – each tributary's face is to the mainstem, not
   to the other tributary. So the constraint is not "all coordinates distinct"; it
   is **strictly increasing across every face**.
6. ~~**Convergence guarantee.**~~ **Settled: the Depression Hierarchy cannot
   produce a divergent depression.** `out_cell`, `out_elev`, `odep` and `parent`
   are each single-valued, `lchild`/`rchild` make it a binary tree, and `flowdirs`
   is *"a value [0,7]"* – one direction per cell
   (`struct Depression` and `GetDepressionHierarchy`,
   `richdem/depressions/depression_hierarchy.hpp`).

   **The guarantee is delivered by tie-break, not by impossibility**, and the
   header says so: *"If a depression has more than one outlet at the same level one
   of them is arbitrarily chosen; hopefully this happens only rarely in natural
   environments."* The same wording covers `out_cell` at equal elevation and
   `pit_cell` at equal depth.

   **"Arbitrary" here means physically undetermined, not nondeterministic.** The
   hierarchy is deterministic by explicit design, on two counts
   (`GetDepressionHierarchy`):

   - parallel nondeterminism is removed outright – *"Since the above runs in
     parallel, the ordering of the seed cells is nondeterministic. Let's fix
     that."* followed by `std::sort` on both seed lists;
   - equal-elevation ties follow a **documented LIFO rule** in the radix heap –
     *"If two or more cells are of equal elevation then the one added last (most
     recently) is returned from the queue first. This ensures that a single
     depression gets all the cells within a flat area."* The tie rule is not
     incidental; it exists so a flat lands wholly in one depression.

   So the same DEM yields the same hierarchy, and reproducibility is not at issue.
   What survives is **reporting**, and it belongs with provenance (§2.7): the
   producer should count the outlets, pits, and flat directions decided by a tie
   rather than by terrain, and report the count. Not because the answer might
   change between runs – it will not – but because those are the places where the
   network rests on a convention instead of on the landscape, and a quantised DEM
   (integer-metre lidar, low-relief glaciated terrain) produces many more of them
   than *"hopefully ... only rarely"* suggests.

   **`r.watershed` is also deterministic** (checked against OSGeo/grass `main`).
   Its A\* heap breaks equal-elevation ties on insertion order, in both
   implementations:

   ```c
   /* raster/r.watershed/ram/do_astar.c */
   static int cmp_pnt(CELL elea, CELL eleb, size_t addeda, size_t addedb)
   {
       if (elea == eleb) { return (addeda < addedb); }
       return (elea < eleb);
   }
   ```

   The `seg/` variant carries the same rule, and neither uses OpenMP – the A\*
   search is single-threaded.

   **But note the conventions are opposite.** RichDEM returns the cell added
   *last* among equal elevations; `r.watershed` returns the one added *first*. Both
   are deterministic; they are not the *same* determinism, so the two will resolve
   a given flat differently. Networks intended for comparison must be produced by
   one routing engine throughout – not one here and one there.

   **Two `r.watershed` flags are load-bearing:**

   - **`-s` is required.** The manual is explicit – *"SFD (D8) flow (default is
     MFD)"* – so the **default is multiple-flow-direction** and would break the
     convergence requirement above. (its `r.watershed` call already
     passes `flags='s'`, so the prototype is correct on this point.)
   - **`-b` must not be used.** *"Beautify flat areas: Flow direction in flat areas
     is modified to look prettier"* – a cosmetic alteration of exactly the flow
     directions this section depends on.

   One residual, and it is not the hierarchy's: the **stream network** is extracted
   separately, and `r.richdem.flowaccumulation` offers MFD algorithms (D-infinity,
   Quinn, Holmgren) alongside D8. The hierarchy's convergence guarantee does not
   extend to an MFD extraction. **Stream-network extraction must be
   single-flow-direction.**
7. **Tributary junctions falling close together.** A property of the drainage
   network itself, present with or without depressions – so the §1.6 measurements
   bound its severity but do not remove it. Untouched by any decision above.

## 6. File format — GeoPackage, with NetCDF as an optional companion

**One GeoPackage is the canonical file.** It is SQLite underneath, which turns out
to serve both audiences from one artifact:

- **Spatial layers, for people and GIS.** Links as `LineString`, nodes as `Point`,
  leaf depressions as `Polygon`, meta-depression pour points as `Point` — very
  nearly what `depressions_to_grass` already builds. The CRS is carried natively and
  the whole network opens in QGIS. That matters more than it sounds: the join
  (§4.2 Gap 3) is the piece most likely to be subtly wrong, and being able to *look*
  at a labelled network over the DEM is the fastest way to catch it.
- **Aspatial tables, for FluvTree.** GPKG permits ordinary attribute tables
  alongside the spatial ones, so the ragged payload lives in long form —
  `(node_id, level, z, B)` and `(dep_label, level, z, A)` — with the hierarchy and
  the provenance block as further tables.

**The dependency argument decides it.** FluvTree depends on numpy, scipy, networkx
and matplotlib (`pyproject.toml`). Requiring GDAL, `fiona` or `geopandas` for the
reader would be a heavy addition for a library whose §3 rule is that it must import
without GRASS. But if the numeric payload sits in aspatial tables with plain
columns, **FluvTree reads a GeoPackage with `sqlite3` from the standard library and
no new dependency at all** — the spatial layers exist for QGIS and never need to be
parsed by FluvTree. Two audiences, one file, neither needing the other's toolchain.

**NetCDF as an optional companion, not a rival.** Its advantage is real but narrow:
the ragged `(values, offsets)` layout maps straight onto `ValleyGeometry` with no
marshalling, and `scipy.io.netcdf_file` reads NetCDF-3 with — again — no new
dependency, scipy already being required. The case for emitting it arises if the
long tables prove slow: `B(x, z)` is per node, so a 10⁵-node network with tens of
levels per node is a table of 10⁶–10⁷ rows, and pulling that through `sqlite3`
row-wise is the one place this design could bite. **Measure before adding it.**

Ruled out: any GRASS-native format as the interchange, which would drag GRASS into
FluvTree's dependencies.

## 7. Open on the FluvTree side

Out of scope for this spec, recorded so they are not lost:

1. Fill–Spill–Merge on the graph, as a per-step diagnostic process.
2. Reclassification as a discrete event inside BDF2 – the solver caches index
   maps, land areas, boundaries and the two-level history at construction.
3. Dynamic `Q`. This lands on an already-recorded gap: "`update_Q` – discharge
   accumulation from the network – NOT clearly ported" (`GRLP-parity-and-gaps.md`).
4. Multi-outlet support (§2.1 item 4).
5. The lake sediment budget – trapping efficiency, and uniform infill versus delta
   progradation. A modelling choice, not a numerical one.
6. **Every boundary condition** (§2.6), including the ones currently read off
   `RiverNetwork.from_arrays`: `S0` / `Q_s_0`, `z_bl`, `x_bl`.
7. **A link with no flow.** The interior stencil divides by `Q[i]`, so discharge
   must be strictly positive – measured on a three-segment chain, `Q` of 1e-9 is
   fine but `Q = 0` gives `Matrix is exactly singular` and `z` all NaN, with
   warnings only. The producer never emits this: it ships contributing *area*,
   which is terrain and always positive. But whatever couples FSM to the solver
   must decide between a positive floor and excluding inactive links from the
   assembly – the latter being the same machinery as item 2.
8. **Initial water state**, and the ongoing water input. FSM needs to be told how
   much water is in the system – `r.richdem.fsm` takes a `water_depth` – and
   nothing currently supplies it: `z_w` appears in this spec only as a field
   FluvTree *writes*. Both the `t = 0` volume per depression and the per-step
   runoff forcing are modelling inputs. Evaporation, which a closed basin needs,
   belongs further out still – it is a property of a model that *uses* this graph,
   not of FluvTree and certainly not of the graph.

## 8. Naming

**`landgraph`.** Free on PyPI (checked: the project JSON returns 404).

The name is deliberately broader than rivers-and-lakes, and correctly so: what this
produces is a graph of the land surface, of which the river-and-lake network is the
first — not the only — thing that can be read off it. Nothing in §2's core tier is
specific to fluvial process.
