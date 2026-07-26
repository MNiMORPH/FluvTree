# TerraPIN ↔ GRLP coupling — design note

Status: **design settled; TerraPIN side ready; coupling to be built here in FluvTree.**
Author: Andy Wickert + Claude Code. Date: 2026-07-26.

## Context

FluvTree is the integrating framework: a directed convergent graph of a river
network. **GRLP** (the gravel long-profile model) is being ported into FluvTree
and evolves each reach's **long profile** (bed elevation `z` down-valley).
**TerraPIN** resolves the **cross-section** at each node — the valley evolving by
incision, aggradation, lateral migration, avulsion, and hillslope (talus)
processes.

The coupling makes GRLP's *prescribed* valley width `B` **emergent** from
TerraPIN, and lets TerraPIN feed back the **lateral** sediment that a 1-D
long-profile model cannot compute (valley widening eroding the walls).

Repos:
- TerraPIN: `~/models/TerraPIN` (PyPI `terrapin-valley`, import `terrapin`).
- GRLP: `~/models/GRLP` (`grlp.grlp.LongProfile`) — being ported into FluvTree;
  its history may be rewritten, so treat the standalone version as a moving target.

## The division of labour

- **GRLP owns the vertical** — the long-profile `z` evolution (sediment transport
  down-valley → `dz`).
- **TerraPIN owns the lateral** — the cross-section: width, wall retreat, talus,
  migration, avulsion.

**One `StandardTerrapin` cross-section per GRLP node.** Each evolves independently
(they do not talk to each other directly; the network coupling is GRLP + the
lateral driver).

## Data exchange (per node, per step)

| direction | quantity | GRLP hook | TerraPIN hook |
|---|---|---|---|
| GRLP → TerraPIN | `dz` (bed change) | `self.z` change | `sweep(x_ch, z_ch + dz)` (auto incise/aggrade) |
| TerraPIN → GRLP | valley width `B` | `set_B(B=array)` (feeds diffusivity) | `compute_valley_width()` |
| TerraPIN → GRLP | **lateral** solid sediment | `set_source_sink_distributed(ssd)` | `sediment_out` after a lateral op |
| external → TerraPIN | channel position (migration) | — | `migrate(x)` / `avulse(x)` |

## The sediment split (the crux — avoids double-counting)

GRLP's `dz` already carries the **vertical** sediment budget. When TerraPIN
applies that `dz` (a vertical `sweep`), the `sediment_out` it reports is *GRLP's,
already counted* — **discard it.** Only the **lateral** `sediment_out` (from
`migrate` / `avulse` / `retreat` — valley widening, wall shedding) is the *new*
contribution that feeds GRLP's `ssd`.

TerraPIN separates these naturally: each op sets `sediment_out` for that op, so
read it right after a *lateral* op and ignore it after the *vertical* sweep.
(Equivalently, drive the BMI with lateral and vertical in **separate `update()`
calls** and read `sediment_out` after the lateral one.)

## The coupling loop (order matters)

```python
# per timestep dt, over the reach's N nodes:
for i, t in enumerate(terrapins):
    lateral_move(i)                         # external driver -> t.migrate/avulse (may be a no-op)
    lateral_sed[i] = t.sediment_out         # SOLID; the lateral contribution
    widths[i] = t.compute_valley_width()
grlp.set_B(B=widths)
grlp.set_source_sink_distributed(to_ssd(lateral_sed, dt, dx))
grlp.evolve_threshold_width_river(1, dt)    # updates grlp.z  (vertical)
for i, t in enumerate(terrapins):
    t.sweep(t.x_ch, grlp.z[i])              # apply dz; sediment_out here is GRLP's -> ignore
```

## Alignment / conventions

- **Porosity `λ_p`.** GRLP has its own `self.lambda_p`; TerraPIN stamps a per-body
  `lambda_p` (bedrock 0, alluvium/colluvium 0.35). Set them consistently.
  TerraPIN's `sediment_out` is already a **solid** volume, which matches GRLP's
  convention (GRLP fluffs by `1/(1 - λ_p)` internally).
- **`z`.** Initialize each `z_ch` = `grlp.z[i]`.
- **Units.** TerraPIN `sediment_out` / `area_out` are **areas** (m² = volume per
  unit down-valley length). Convert to GRLP's `ssd` (a source/sink rate) using the
  node spacing `dx` and step `dt`. **[open — reconcile precisely.]**
- **Emergent `valley_width`** is measured wall-to-wall against the *confining* walls
  (bedrock + original alluvium), so floor deposits (belt, floodplain, talus) do not
  narrow it. Measured at `z_ch` (width-at-bed).

## The BMI (drive TerraPIN from FluvTree)

`terrapin.bmi.BmiStandardTerrapin` wraps **one** cross-section (scalar grid;
FluvTree composes N). gFlex-style: subclasses `bmipy.Bmi` when installed, else a
plain object (usable without the optional dep).

- **Inputs** (`set_value` then `update()`):
  `channel_bottom__elevation` → vertical `sweep`; `channel__x_position` → `migrate`.
- **Outputs** (`get_value`): `valley__width`,
  `channel_solid_sediment__volume` (solid, the river's load),
  `channel_bulk_sediment__area`.
- **Time** is nominal (event-driven; each `update()` applies the queued move).
- **Config** is a dict or JSON path (no yaml dependency): domain, `bedrock_top`,
  `surface`, `channel {x,width,depth,z}`, `repose`, `lambda_p`, `porosities`.

FluvTree can drive TerraPIN either through the BMI (one per node) or by calling
`StandardTerrapin` directly — the BMI just standardizes it.

## Open items

- **Lateral driver.** An *input* to the coupling (stochastic migration, or a rule).
  A simple default is worth having for testing — candidate FluvTree component.
- **Units reconciliation** between TerraPIN areas and GRLP `ssd` (see above).
- **`valley_width` level** for aggraded valleys (currently at `z_ch`).
- **Porting GRLP into FluvTree** — the coupling targets FluvTree's GRLP.
- **Optional `[bmi]` extra** in TerraPIN's `pyproject.toml` (bmipy) — not yet added.
- **TerraPIN talus/hillslope refinements** — TerraPIN issue #10 (non-flat strath,
  overflow, colluvium age, undercut→river displacement, local wall-top).

## What's ready on the TerraPIN side (already built + tested)

- `sweep(x1, z1)` — the unified motion (incise / migrate / aggrade as cases; auto
  erode/deposit); vertical `sweep` is the `dz` hook.
- `migrate` / `avulse` / `retreat` — lateral + hillslope, producing lateral sediment.
- `compute_valley_width()` — emergent, asymmetric, wall-to-wall.
- Per-body `lambda_p`; `sediment_out` = net **solid**, `area_out` = net bulk.
- `BmiStandardTerrapin` — the BMI above.
