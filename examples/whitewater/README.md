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

**Parameters.** `D = 0.3 mm`, Manning `n = 0.03`, `tau_crit_bank = 5 Pa` — the values
in FluvTree's own sand tests. Discharge is contributing area × `1e-7 m/s` with
intermittency `0.05`, giving flood discharges of 4.4–82 m³/s and time-averaged
0.22–4.1 m³/s. That is ~158 mm/yr of runoff, on the low side for southeastern
Minnesota. Both are illustrative.

**Response time.** Bed diffusivity `kappa = p k_Qs I Q S^(p-1) / ((1-lambda_p) B)` is
1.5–2.6e-2 m²/s on the main stem, so `L^2/kappa` is ~770 yr over 25 km. Spin-up
e-folds at ~450 yr, consistent. If the response looks fast, the levers are the reach
slopes (`kappa ~ S^(-1/6)`, so flatter reaches are *more* diffusive) and `k_Qs`
through `tau_crit_bank`.

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
