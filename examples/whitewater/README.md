# Sand-bed lower Whitewater — a sediment pulse, then a dam

**1860–1940** post-settlement land clearance delivers sediment along the valley,
aggrading it by ~3.6 m. **1933–1938** a lock and dam closes on the Mississippi,
raising base level at the Whitewater's mouth by 2.5 m. The river then has until
**2360** to clear what was dumped on it, and largely cannot.

| year | mean stored | peak |
|---|---|---|
| 1932 (pre-dam) | +2.71 m | +3.23 m |
| 1940 | +3.08 m | +3.59 m |
| 2100 | +2.99 m | +3.30 m |
| 2360 | **+2.77 m** | +2.95 m |

Without the dam the same pulse relaxes to **+1.05 m** over the same interval. Raising
base level flattens the gradient near the mouth, transport capacity falls, and the
post-settlement alluvium is **locked in** rather than exported. That contrast is the
point of the run.

The Whitewater is a classic post-settlement-alluvium valley and enters the Mississippi
near Winona, so both forcings are in the spirit of its history — but **the supply rate
is illustrative, not a measured value**, and the dam is a generic 2.5 m rise rather
than a specific structure.

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
