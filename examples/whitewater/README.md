# Sand-bed lower Whitewater — a sediment pulse, then a dam

**1860–1940** post-settlement land clearance delivers sediment along the valley,
aggrading it by ~3.6 m. **1933–1938** a lock and dam closes on the Mississippi,
raising base level at the Whitewater's mouth by 2.5 m. The river then has until
**2360** to clear what was dumped on it, and largely cannot.

| year | mean stored | peak |
|---|---|---|
| 1932 (pre-dam) | +3.08 m | +3.24 m |
| 1940 | +3.46 m | +3.60 m |
| 2100 | +3.45 m | +3.64 m |
| 2360 | **+3.35 m** | +3.64 m |

**The valley keeps almost all of it** — 0.11 m removed in 420 years. That is the
result of taking the valley width seriously (below): the bed diffusivity is low
enough that a metres-thick alluvial fill is effectively permanent on human
timescales. It is also what the Whitewater shows: the post-settlement alluvium is
still there.

The Whitewater is a classic post-settlement-alluvium valley and enters the Mississippi
near Winona, so both forcings are in the spirit of its history — but **the supply rate
is illustrative, not a measured value**, and the dam is a generic 2.5 m rise rather
than a specific structure.

**1985–1995** water discharge rises 15%, also ramped. Its effect is small — +2.74 m
stored in 2360 against +2.77 m without it. Once base level is raised the system is
**base-level-limited, not transport-limited**: the river's capacity to work is no
longer what decides how much sediment stays.

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

Five links, 36 km, three heads and two binary junctions, one outlet. Reach slopes
0.00053–0.00279, contributing areas 30–830 km². Trout Creek enters as one of the
tributaries.

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
