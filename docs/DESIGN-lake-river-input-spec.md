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
`storage_volume` is documented as *"per unit valley length"* (`valley.py:151`), so
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
`for j in up_ids` with no cap (`diffusion.py:254`), so an *n*-way shared surface
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
(`diffusion.py:284`, and the same form at `188-197`). So the assignment rule is

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
markers ... not shared `z` degrees of freedom"* (`diffusion.py:18-22`) – and the
face sits at the midpoint of the gap (`land = B[-1] * 0.5 * ((x[-1] - x[-2]) +
(x[dj][0] - x[-1]))`, `diffusion.py:269`).

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

## 2. Required deliverables

Constraints below marked *(verified)* are read off the current code, with the line
that enforces them.

### 2.1 Topology, per segment

1. `seg_id` – integer, stable.
2. `upstream_segment_IDs`, `downstream_segment_IDs` – the form
   `RiverNetwork.from_segment_lists` takes.
3. **At most one downstream segment.** *(verified: `network.py:161` raises
   "non-convergent topology".)* The joined graph must still be a convergent tree.
   This is the locked FluvTree invariant (`DESIGN-structure-and-naming.md`, Scope),
   and lakes must not break it.
4. **One outlet.** *(verified: `network.py:111`, `diffusion.py:122`,
   `advection.py:34,45` all take `mouth_segments()[0]`.)* Multiple outlets are not
   supported today; see §7.
5. `dep_label` – **per link**, the hierarchy node whose **own region** contains it,
   per §1.5; null outside every depression. Per node only on an uncut margin
   (§1.6), which the network holds as an ordinary per-segment array. This is the immutable field the
   producer ships. `seg_type` (`"river"` / `"lake"`) is **not** shipped as
   structure – FluvTree derives it each step from the active antichain – though an
   initial classification is useful as a t=0 state and as a check on the producer.

### 2.2 Along-flow geometry, per segment, per node

6. `x` – along-flow distance [m], monotonic.
7. `z_bed` – bed or basin-floor elevation [m].
8. **Node-count minima.** *(verified: `diffusion.py:174,179`.)* A segment
   downstream of a confluence needs ≥ 3 nodes; a tributary entering a confluence
   needs ≥ 2. The producer must honour this when it chooses breakpoints, or
   FluvTree raises at solver construction.

### 2.3 Cross-section geometry, per node, ragged

9. `z_levels`, `B_levels` – ascending `z`, paired widths. **≥ 2 levels per node**
   *(verified: `valley.py:51`.)*
10. `channel_width` – the backup width below the lowest tabulated level.

Store these as flat `z`, flat `B`, plus an `offsets` index. That is exactly
`ValleyGeometry`'s own in-memory layout, so the reader is close to a copy.

**The same table serves water and sediment, at two different porosities.** A
depression fills with water *or* with sediment through the identical `B(x, z)`
curve; only the solid fraction differs. `ValleyGeometry` bakes `(1 - lambda_p)`
into both `storage_volume` and `storage_jacobian` (`valley.py:149,181`), so a
water-storage query needs the `lambda_p = 0` path. The producer therefore ships
**one geometry and no porosity** – `lambda_p` is FluvTree's sediment parameter,
not a property of the container, and must not be folded into the delivered
widths.

### 2.4 Hydrology, per node

11. `A` – contributing drainage area [m²].

**Deliver area, not only discharge.** FluvTree currently treats `Q` as fixed for
the run *(verified: `diffusion.py:92`, "Constant-in-time reach state")*. Lakes
break that: an endorheic basin passes no water downstream until it spills, and
evaporation removes water. Area is time-invariant and routing-independent;
discharge is neither. Give FluvTree the invariant quantity and let it recompute
`Q` as routing changes.

### 2.5 The depression hierarchy

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

### 2.6 Boundaries

20. outlet segment and node, `x_bl`, `z_bl`.
21. per-head `S0` or a network-level `Q_s_0` (the existing
    `RiverNetwork.from_arrays` boundary set).

### 2.7 Georeferencing and provenance

22. `x`, `y` per node, and the CRS.
23. **Provenance block**: DEM source and resolution, computational region,
    RichDEM / `r.richdem` versions, and **every threshold applied, with the count
    of features it kept and dropped**. Include the **arbitrary-tie-break counts**
    per §5.6 – outlets, pit cells, and flat directions resolved by an arbitrary
    pick – and the flow-routing algorithm used, which must be single-direction.

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
an accumulation cut of 1e7, and line 50 is a bare shell invocation
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

3. **Minimum depression size.** `dephier` finds every pit, including single-cell
   DEM noise. A filter is unavoidable, and the criterion – volume, area, depth, or
   a combination – **is yours to set.** I am not choosing one.
4. **Endorheic basins and the ocean.** How `ocean_parent` maps onto the outlet
   boundary, and what a basin that never spills means for downstream `Q`.
5. **Snapping tolerance.** Sill cells will rarely fall exactly on a stream-reach
   breakpoint.
6. ~~**Convergence guarantee.**~~ **Settled: the Depression Hierarchy cannot
   produce a divergent depression.** `out_cell`, `out_elev`, `odep` and `parent`
   are each single-valued, `lchild`/`rchild` make it a binary tree, and `flowdirs`
   is *"a value [0,7]"* – one direction per cell
   (`richdem/depressions/depression_hierarchy.hpp:38-80, 269`).

   **The guarantee is delivered by tie-break, not by impossibility**, and the
   header says so: *"If a depression has more than one outlet at the same level one
   of them is arbitrarily chosen; hopefully this happens only rarely in natural
   environments."* The same wording covers `out_cell` at equal elevation and
   `pit_cell` at equal depth.

   **"Arbitrary" here means physically undetermined, not nondeterministic.** The
   hierarchy is deterministic by explicit design, on two counts
   (`depression_hierarchy.hpp:395-407`):

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
     convergence requirement above. (`scripts/drainageNetworkGRASS.py:41` already
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

## 6. File format

**One self-describing NetCDF file.** The ragged cross-sections are the awkward
part, and NetCDF carries them cleanly as flat value arrays plus an `offsets`
index – mirroring `ValleyGeometry`'s layout exactly. The provenance block of §2.7
becomes file attributes.

Alternative: a directory of Parquet + a JSON manifest. Not recommended: any
GRASS-native format as the interchange, which would drag GRASS into FluvTree's
dependencies.

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

## 8. Naming

The producer package is unnamed. That is Andy's call.
