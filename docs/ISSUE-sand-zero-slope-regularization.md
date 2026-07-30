# Issue (draft): regularize the sand-bed (cohesive) engine as slope → 0

Status: **parked** — idea path, to be filed as a GitHub issue once alluvial-core has a
remote. Set aside deliberately to stay on the FluvTree integration; the interim
(threshold taper) concept is kept, the hydraulic-radius multiplier is the preferred
refinement for later (Andy has further thoughts to add).

## Problem

The sand-bed closure has flux exponent `p = 5/6`, so the implicit sediment-flux
**conductance goes as `S^(p-1) = S^(-1/6)`, which diverges as `S → 0`**. A flat
reach (`S = 0`, e.g. a fresh flat bed or a reach at base level) makes the matrix
exactly singular (`spsolve` → NaN). Gravel (`p = 7/6`, `S^(+1/6) → 0`) has no such
problem; this is sand-specific. Sand runs today only from a **sloped initial
condition**, which is a workaround, not a fix.

## Root cause (physics, not just numerics)

The threshold-width closure pins the bed shear stress at `(1+ε)·τ_crit_bank`
**using depth** (`τ_b = ρg·h·S`, with `h ∝ 1/S`). That assumes a wide channel
(`R_h ≈ h`). As `S → 0` the closed form drives `b → 0` and `h → ∞` to hold that
stress — physically absurd (a zero-width, infinitely-deep slot). Verified:

```
 S       b [m]      h [m]        τ_b(depth)   τ_b(R_h = bh/(b+2h))
 3e-3    27         0.08         2.40 Pa      2.39 Pa
 1e-4    0.51       2.45         2.40 Pa      0.23 Pa   <- already sub-threshold
 1e-6    0.0024     245          2.40 Pa      1.2e-5 Pa
 0       0          ∞            2.40 Pa      0
```

The pinning is a **wide-channel artifact**: once `b ≪ h`, `R_h = bh/(b+2h) → b/2`,
not `h`, so the true bed stress collapses. Critical Shields stress for 1 mm sand:
`τ_c = τ*_c·(ρ_s-ρ)·g·D ≈ 0.5–1.0 Pa` (0.73 Pa at τ*_c = 0.045).

## Approach A — threshold taper (interim / "min-Q_s" machinery)

Multiply transport (and conductance) by a smooth factor of the **hydraulic-radius**
Shields stress vs the grain threshold: 1 where wide, → 0 as `τ_b(R_h) → τ_c`.

```python
def taper(closure, Q, S, tau_star_c=0.045, r_full=2.0, rho=1000., g=9.805):
    b, h = closure.channel_width(Q, S), closure.channel_depth(S)
    Rh = b*h/(b+2*h)
    tau_c = tau_star_c*(closure.rho_s-closure.rho)*g*closure.D
    r = rho*g*Rh*S / tau_c
    t = np.clip((r-1)/(r_full-1), 0, 1)
    return t*t*(3-2*t)          # smoothstep; discharge-aware; = 1 in the normal regime
```

Pros: exactly 1 in the normal regime (no change to the calibrated model), physically
motivated shutoff, discharge-aware (big rivers transport to lower slope). Cons: two
science knobs (`τ*_c`, `r_full`); hard clamp needs a `0·∞` guard at exact `S=0`.

## Approach B — hydraulic-radius multiplier (preferred, for later)

The depth→`R_h` correction is a **pure geometric ratio** of the width and depth the
model already solves for, applied as a Picard-frozen multiplier — no reconstructed
slope-stress, no new nonlinearity:

```python
def m_Rh(closure, Q, S, S_eps=1e-12):
    Ss = np.where(S > 0, S, S_eps)
    b, h = closure.channel_width(Q, Ss), closure.channel_depth(Ss)
    return b/(b+2*h)           # R_h/h ; ->1 wide, ->0 narrow-deep
# conductance:  C0*Q * max(S,S_eps)**(p-1) * m_Rh(...)   (S_eps floors S^(p-1) only)
```

Verified behaviour (sand closure, D=1mm):

```
 S       b/h        m=R_h/h    conductance uncorr.   corrected
 3e-3    1.5e4      0.9999     3.69e13               3.69e13
 3e-4    99         0.980      5.41e13               5.30e13
 1e-4    9.2        0.82       6.50e13               5.34e13
 1e-6    4e-4       0.0002     1.4e14 (blowing up)   2.99e10
 0       —          0.000      ∞                     finite
```

Pros: **parameter-free** (pure geometry); `m ≈ 1` for any real wide channel
(≤0.2% for b/h > 50) so the calibrated model is essentially untouched, and the small
residual *is* the true hydraulic radius the depth form dropped; regularizes the
singularity; no new nonlinearity (frozen coefficient like `C1`). Nice consequence:
`m < 1` steepens the equilibrium slope needed to carry the flux through a
low-conveyance reach, so the model **maintains a minimum transporting slope** instead
of flattening to zero — arguably the right physics near base level.

## Open questions / decisions

- Power on `m`: `m^1` (direct stress correction) vs `m^{3/2}` (if propagated through
  the transport–stress relation). Only sharpens the transition; ≈1 in the normal
  regime either way. **Andy has further thoughts here.**
- Wiring: `TransportClosure.flux_multiplier(Q, S)` returning 1 by default (gravel
  stays **bit-exact** with GRLP), sand returns `b/(b+2h)`; multiplied into
  `compute_Q_s` (flux) *and* `_face_conductance` / `C1` (conductance), consistently,
  with the `S_eps` floor.
- Keep A or replace with B? B is cleaner and parameter-free; A is the fallback.

## Acceptance

- A sand reach relaxes to a true zero-slope base level with no NaN / no singular
  matrix.
- Gravel remains bit-identical to GRLP (the multiplier is 1 for the gravel closure).
- Steady state in the normal (wide-channel) regime unchanged to <~1% vs current sand.
