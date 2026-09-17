# Sand-bed lower Whitewater — a sediment pulse and its aftermath

An 80-year sediment supply pulse aggrades the valley by ~3.6 m; 420 years of
recovery clear about two thirds of it. The Whitewater is a classic
post-settlement-alluvium valley, so the scenario is in the spirit of that history —
but **the supply rate is illustrative, not a measured value.**

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

**Parameters.** `D = 0.3 mm`, Manning `n = 0.03`, `tau_crit_bank = 5 Pa` — the values
in FluvTree's own sand tests. Discharge is contributing area × `1e-7 m/s` with
intermittency `0.05`; those two are illustrative.

## Why sand only

`DiffusionSolver` takes **one closure for the whole network** — `p` and `k_Qs` are
single scalars fixed at construction ("Uniform closure -> one scalar for the whole
network"). Gravel and sand differ precisely in `p` (7/6 versus 5/6), so a mixed
gravel–sand network is **not** a matter of configuration: it needs per-segment
closure properties threaded through the assembly. That is a real solver change, not
a small one, so this example stays entirely in the sand-bed reach.
