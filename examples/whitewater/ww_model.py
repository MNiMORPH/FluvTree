"""Set up the sand-bed lower Whitewater as a FluvTree network."""
import sys, numpy as np
sys.path.insert(0, "/home/awickert/models/FluvTree/src")
sys.path.insert(0, "/tmp/claude-1000/-home-awickert-models-FluvTree/4fd69dbf-d3b7-42b4-aeb8-7a9098b98b09/scratchpad")
from build_ww import build
from fluvtree.processes import build_network
from fluvtree.closures.sand import SandClosure

# --- discharge, from USGS 05376800 "Whitewater River near Beaver, MN" -------
# Drainage area 271 mi2 = 702 km2, 4.7 mi north of Elba.
#   annual means, 14 water years 1976-1998: mean 168.4 cfs = 4.77 m3/s
#     -> 214 mm/yr of runoff
#   annual peaks, 38 water years 1938-1999, Weibull plotting positions:
#     Q_2 = 3355 cfs = 95.0 m3/s   (Q_1.5 = 2210 cfs = 62.6 m3/s)
# Take the formative discharge as Q_2 and set intermittency so that I*Q
# reproduces the mean annual flow: I = 4.77 / 95.0 = 0.0502.
# Scaled to other reaches by contributing area, i.e. constant specific discharge.
Q2_BEAVER, A_BEAVER = 95.0, 702.0e6
RUNOFF   = Q2_BEAVER / A_BEAVER      # 1.353e-7 m/s -> formative discharge
INTERMIT = 4.77 / Q2_BEAVER          # 0.0502

# --- valley width ----------------------------------------------------------
# B in the solver is the width sediment is spread over — the VALLEY, not the
# channel. (An earlier version of this script used the closure's channel_width,
# which is a diagnostic, not this. That made the bed roughly an order of
# magnitude too responsive.) Ramped linearly along the network.
B_UPSTREAM, B_MOUTH = 500.0, 1000.0
D, MANN, TAU = 0.3e-3, 0.03, 5.0     # the values in FluvTree's own sand tests
CLOSURE = SandClosure(D=D, n=MANN, tau_crit_bank=TAU)

def pava_decreasing(y, w=None):
    """Closest non-increasing sequence in least squares (pool adjacent violators).

    A river long profile falls downstream. The raw DEM does not — 25% of 5 m
    intervals rise, which is real microtopography plus survey noise, but it makes
    local slopes meaningless and the sand closure's width goes as S**(7/6).
    Isotonic regression is the MINIMAL adjustment that restores the physical
    constraint: no smoothing window to choose, no information discarded beyond the
    monotonicity violation itself.
    """
    y = np.asarray(y, float).copy()
    w = np.ones_like(y) if w is None else np.asarray(w, float).copy()
    vals, wts, cnt = [], [], []
    for yi, wi in zip(y, w):
        vals.append(yi); wts.append(wi); cnt.append(1)
        while len(vals) > 1 and vals[-2] < vals[-1]:      # violation: rising
            v2, w2, c2 = vals.pop(), wts.pop(), cnt.pop()
            v1, w1, c1 = vals.pop(), wts.pop(), cnt.pop()
            vals.append((v1 * w1 + v2 * w2) / (w1 + w2))
            wts.append(w1 + w2); cnt.append(c1 + c2)
    out = np.empty_like(y); i = 0
    for v, c in zip(vals, cnt):
        out[i:i + c] = v; i += c
    return out


def make(json_path, dx=100.0):
    order, idx, to, up, segs = build(json_path, dx=dx)
    up_ids   = [[idx[u] for u in up.get(c, [])] for c in order]
    down_ids = [[idx[to[c]]] if to[c] != 0 else [] for c in order]

    x  = [segs[c]["x"] for c in order]
    z_dem = [segs[c]["z"].copy() for c in order]          # raw DEM trace, kept for reference
    z_iso = [pava_decreasing(zz) for zz in z_dem]
    Q  = [np.maximum.accumulate(np.abs(segs[c]["A"])) * RUNOFF for c in order]

    # width from the closure's own hydraulic geometry, at the local bed slope
    # Reach-average slope, from the data (end-to-end over the isotonic profile).
    S_reach = [float((z_iso[i][0] - z_iso[i][-1]) / (x[i][-1] - x[i][0]))
               for i in range(len(order))]

    # Valley width, ramped linearly from the upstream end of the network to the mouth.
    x_up = min(float(xx[0]) for xx in x)
    B = [B_UPSTREAM + (B_MOUTH - B_UPSTREAM) * (xx - x_up) / (0.0 - x_up) for xx in x]

    # head slope from the local bed, outlet base level from the DEM
    # INITIAL CONDITION. The raw DEM trace cannot be used directly: 25% of its
    # intervals rise, and isotonic regression removes that only by pooling violators
    # into exactly-flat runs — and the sand closure is SINGULAR at S = 0
    # (conductance ~ S**(-1/6)). So build the profile by integrating each link's own
    # reach-average slope upstream from base level. That keeps the real relief
    # structure — steep tributaries, gentle main stem — while being strictly
    # decreasing everywhere, and every slope in it came from the data.
    mouth = [i for i, d in enumerate(down_ids) if not d][0]
    z_bl_dem = float(z_dem[mouth][-1])
    z = [None] * len(order)

    def fill(i, z_end):
        """z_end is the elevation at this link's DOWNSTREAM end."""
        xi = x[i]
        z[i] = z_end + S_reach[i] * (xi[-1] - xi)
        for u in up_ids[i]:
            z[u] = None                       # placeholder; filled below
            fill(u, z[i][0] + S_reach[i] * (x[i][0] - x[u][-1]))

    fill(mouth, z_bl_dem)

    heads = [i for i, u in enumerate(up_ids) if not u]
    S0 = [float(max(S_reach[i], 1e-5)) for i in sorted(heads)]
    x_bl = float(x[mouth][-1] + (x[mouth][-1] - x[mouth][-2]))
    z_bl = float(z[mouth][-1] - S_reach[mouth] * (x_bl - x[mouth][-1]))

    net = build_network(up_ids, down_ids, x, z, Q, B,
                        S0=S0, x_bl=x_bl, z_bl=z_bl)
    return net, order, idx, dict(x=x, z=z, z_dem=z_dem, z_iso=z_iso, Q=Q, B=B,
                                 S_reach=S_reach, S0=S0, x_bl=x_bl, z_bl=z_bl)

if __name__ == "__main__":
    net, order, idx, m = make(sys.argv[1])
    print(f"segments {len(order)}  (cats {order})")
    dz = np.concatenate([a - b for a, b in zip(m["z_dem"], m["z"])])
    print(f"IC vs raw DEM trace: max |dz| {np.abs(dz).max():.2f} m, "
          f"RMS {np.sqrt((dz**2).mean()):.2f} m")
    mn = min(float(np.min(np.diff(-zz))) for zz in m["z"])
    print(f"strictly decreasing everywhere: {mn > 0}  (min downstream drop {mn:.5f} m)")
    for i, c in enumerate(order):
        print(f"  seg {i} (cat {c:>3}): {len(m['x'][i]):>4} nodes  "
              f"Q {m['Q'][i][0]:6.2f}..{m['Q'][i][-1]:6.2f} m3/s   "
              f"S {m['S_reach'][i]:.5f}  B_valley {m['B'][i].min():6.0f}..{m['B'][i].max():6.0f} m")
    print(f"  S0 (per head) {['%.5f' % s for s in m['S0']]}")
    print(f"  x_bl {m['x_bl']:.1f}   z_bl {m['z_bl']:.3f}")
