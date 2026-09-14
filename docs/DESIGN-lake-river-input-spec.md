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
5. `seg_type` – `"river"` or `"lake"`, the initial classification.

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

Item 15 is the one that does not exist in any current output and is the crux of
the join (§5): the hierarchy knows the sill as a raster cell; FluvTree needs it as
a position on a reach.

### 2.6 Boundaries

18. outlet segment and node, `x_bl`, `z_bl`.
19. per-head `S0` or a network-level `Q_s_0` (the existing
    `RiverNetwork.from_arrays` boundary set).

### 2.7 Georeferencing and provenance

20. `x`, `y` per node, and the CRS.
21. **Provenance block**: DEM source and resolution, computational region,
    RichDEM / `r.richdem` versions, and **every threshold applied, with the count
    of features it kept and dropped**.

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

1. **Breakpoints.** Segments break at confluences; do they also break at lake
   margins and at sills? (Almost certainly yes, but it interacts with the
   node-count minima in §2.2.)
2. **Leaf depression → one lake segment.** A *meta*-depression is then not a
   segment at all: it is a potential future merge, carried in the tree.
3. **Minimum depression size.** `dephier` finds every pit, including single-cell
   DEM noise. A filter is unavoidable, and the criterion – volume, area, depth, or
   a combination – **is yours to set.** I am not choosing one.
4. **Endorheic basins and the ocean.** How `ocean_parent` maps onto the outlet
   boundary, and what a basin that never spills means for downstream `Q`.
5. **Snapping tolerance.** Sill cells will rarely fall exactly on a stream-reach
   breakpoint.
6. **Convergence guarantee.** A depression spilling toward two neighbours, or an
   ambiguously routed flat, breaks the convergent-tree invariant. The join needs
   an explicit check and a documented resolution rule, not a silent pick.

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
