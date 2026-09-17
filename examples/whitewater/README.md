# Sand-bed lower Whitewater — a sediment pulse, then a dam

**1860–1940** post-settlement land clearance delivers sediment along the valley,
aggrading it by ~3.6 m. **1933–1938** a lock and dam closes on the Mississippi,
raising base level at the Whitewater's mouth by an effective **2 m** — less than the
pool itself, because the Whitewater progrades a delta into it. The river then has until
**2360** to clear what was dumped on it, and largely cannot.

| year | mean stored | peak | base level | Q |
|---|---|---|---|---|
| 1932 (pre-dam) | +3.13 m | +3.24 m | — | ×1.00 |
| 1940 | +3.50 m | +3.60 m | +2.00 m | ×1.00 |
| 1990 | +3.49 m | +3.65 m | +2.00 m | ×1.10 |
| 2030 | +3.37 m | +3.64 m | +2.00 m | ×1.20 |
| 2360 | **+2.47 m** | +3.43 m | +2.00 m | ×1.20 |

**The valley keeps most of it** — 1.03 m removed in 420 years, from a peak of 3.50 m. That is the
result of taking the valley width seriously (below): the bed diffusivity is low
enough that a metres-thick alluvial fill is effectively permanent on human
timescales. It is also what the Whitewater shows: the post-settlement alluvium is
still there.

The Whitewater is a classic post-settlement-alluvium valley and enters the Mississippi
near Winona, so both forcings are in the spirit of its history — but **the supply rate
is illustrative, not a measured value**, and the dam is a generic 2.5 m rise rather
than a specific structure.

**Discharge rises** to +10% by 1990 and +20% by 2030, then holds. Those two anchors
define a single straight line, which reaches zero in 1950 — so the ramp is clipped
there rather than given an invented start date.

**The base-level rise is ramped over five years, not stepped.** FluvTree cannot
simulate ponding, so an instantaneous change in `z_bl` would be a numerical shock
standing in for a physical process the model does not represent.

Run in order:

```
python3 ww_scenario.py wwsand2m.json     # spin up, run the scenario, store frames
python3 ww_animate.py  out.mp4           # render
```

`wwsand2m.json` comes from `v.stream.networkx` on the `streams_sand` vector in the
`Whitewater` GRASS location, sampling `accumulation` and `Whitewater_5m_nonulls`.

## The network

Nine links, five heads, one outlet. Reaches identified from `streams_all` by gauge
position, not by guessing from direction:

| cat | reach | extent | valley width |
|---|---|---|---|
| 6 | Middle Fork — the mainstem above Elba | lowest 1 km | 500 m (ramp start) |
| 14 | **North Fork** — USGS 05376000 sits 77 m from it; its 270 km² matches the gauge's stated 101 mi² = 262 km² | lowest 5 km | 200 m |
| 8 | South Fork | lowest 5 km | 200 m |
| 16 | Elba → South Fork confluence | full | ramp |
| 17, 19, 20 | mainstem to the mouth | full | ramp to 1000 m |
| 18 | **Beaver Creek** — joins 0.67 km north of gauge 05376800, described as "0.5 mi above mouth of Beaver Creek" | lowest 4 km (it is only ~3.6 km, so uncut) | 200 m |
| 15 | Trout Creek | lowest 4 km | 50 → 200 m |

The sand reach begins ~1 km above Elba on the Middle Fork, so the forks join *inside*
the modelled domain.

## Head sediment supply is pinned as a FLUX, not a slope

FluvTree's head boundary is a slope `S0`, and the flux it implies is
`Q_s = k_Qs · I · Q · S0^p` — **proportional to Q**. So raising discharge 20% with
`S0` fixed silently raises the upstream sediment supply by 20% too, which is not what
"supply reverts to pre-disturbance" means: upland sediment yield is set by land use,
not by how much water the river carries.

The scenario therefore pins `Q_s` at its pre-disturbance value and recomputes `S0`
from it whenever `Q` changes. This is not cosmetic — it roughly triples what the
river removes, because the extra water now does net erosive work instead of arriving
with proportionally more sediment.

## Choices that had to be made, and why

**Sample the raster at the vector's own resolution.** `streams_sand` has ~2 m vertex
spacing. Sampling the *5 m* accumulation against it put **34% of nodes off-channel**,
with contributing area collapsing to 25 m² — one cell — in the middle of the main
stem. Sampling the 2 m accumulation instead gives **0%**. `Q = 0` is fatal: the
interior stencil divides by `Q[i]`, so it yields an exactly singular matrix and `z`
all NaN, with only a warning.

**The junction is a face, not a shared node.** Adjacent links initially shared an `x`
coordinate, which is a zero-length face and the same singular-matrix failure. The
junction vertex goes to the *downstream* link and each upstream link ends one spacing
earlier.

**The raw DEM trace cannot be the initial condition.** 25% of its intervals rise.
Isotonic regression is the minimal fix for that, but it removes violations by pooling
them into exactly-flat runs — and the sand closure is **singular at `S = 0`**
(conductance goes as `S^(-1/6)`). So the profile is built by integrating each link's
own reach-average slope upstream from base level: strictly decreasing everywhere,
real relief structure preserved, every slope taken from the data. It differs from the
DEM trace by 0.68 m RMS.

**Width comes from the closure, at the reach-average slope.** `SandClosure` computes
channel width itself. Evaluated at node-to-node DEM slopes it gave widths from 0 to
2282 m, because the law goes as `S^(7/6)`. One slope per link — from the data, no
smoothing window chosen — gives 18–133 m.

**Discharge, from the gauge.** USGS **05376800**, *Whitewater River near Beaver, MN*,
271 mi² = 702 km², 4.7 mi north of Elba.

- 14 water years of annual means (1976–1998): **168.4 cfs = 4.77 m³/s**, i.e.
  **214 mm/yr** of runoff
- 38 water years of annual peaks (1938–1999), Weibull plotting positions:
  **Q₂ = 3355 cfs = 95.0 m³/s**; Q₁.₅ = 2210 cfs = 62.6 m³/s

The formative discharge is taken as Q₂, and intermittency set so `I·Q` reproduces the
mean annual flow: **I = 4.77 / 95.0 = 0.0502**. Scaled to other reaches by
contributing area (constant specific discharge), which gives 4.1–111 m³/s across the
network.

**Valley width is `B`, not channel width.** `B` in the solver is the width sediment is
spread over. An earlier version of this script used `SandClosure.channel_width`, which
is a *diagnostic* — 18–133 m instead of 500–1000 m. That made the bed roughly an
order of magnitude too responsive and is the single largest correction this example
has had. Width now ramps linearly from 500 m at the upstream end to 1000 m at the
mouth.

**Response time.** With the valley width, spin-up takes **50,000 years** to
equilibrate, against ~5,000 with a channel width. That is why the post-settlement fill
does not clear.

**Why raising discharge does so little, and a trap.** In this closure
`B = k n Q S^(7/6)`, so **B is proportional to Q**, while the bed diffusivity goes as
`Q/B`. Scaling discharge *with the closure's own width relation* therefore leaves the
response **exactly unchanged** — the widening cancels the extra transport. The run
holds `B` fixed, which is the defensible short-term case for cohesive banks, so the
15% does raise `kappa` by 15%. It still barely moves the outcome, for the separate
reason above: base level, not capacity, is the limit.

## Why sand only

`DiffusionSolver` takes **one closure for the whole network** — `p` and `k_Qs` are
single scalars fixed at construction ("Uniform closure -> one scalar for the whole
network"). Gravel and sand differ precisely in `p` (7/6 versus 5/6), so a mixed
gravel–sand network is **not** a matter of configuration: it needs per-segment
closure properties threaded through the assembly. That is a real solver change, not
a small one, so this example stays entirely in the sand-bed reach.
