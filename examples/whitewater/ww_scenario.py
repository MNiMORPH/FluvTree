"""Sediment-supply pulse and recovery on the sand-bed lower Whitewater.

Phase 1  spin-up to steady state (discarded; the pre-disturbance river)
Phase 2  PULSE: uniform sediment delivery to the valley floor, 80 yr.
         The Whitewater is a classic post-settlement-alluvium valley — upland
         agricultural erosion delivered sediment along the whole valley, burying
         the floor by metres. The rate here is chosen to produce metre-scale
         aggradation over decades; it is illustrative, not a measured value.
Phase 3  RECOVERY: delivery stops, the river reworks what was dumped on it.
"""
import sys, numpy as np, warnings, pickle
sys.path.insert(0, "/home/awickert/models/FluvTree/src")
sys.path.insert(0, "/tmp/claude-1000/-home-awickert-models-FluvTree/4fd69dbf-d3b7-42b4-aeb8-7a9098b98b09/scratchpad")
from ww_model import make, CLOSURE, INTERMIT, RUNOFF
from fluvtree.solvers.diffusion import DiffusionSolver

YR = 3.15576e7
START_YEAR = 1860               # post-settlement land clearance under way
PULSE_END  = 1940               # supply falls away as conservation takes hold
END_YEAR   = 2360
PULSE_M_PER_YR = 0.045          # illustrative supply rate during the pulse
# Lock and dam closure. The Whitewater enters the Mississippi near Winona, so a
# navigation pool raises its base level. Ramped over five years rather than
# imposed instantly: FluvTree has no ponding, so a step change in z_bl would be a
# numerical shock standing in for a physical process it cannot represent.
DAM_START, DAM_END, DAM_RISE = 1933, 1938, 2.5
# A 15% rise in water discharge, ramped across a decade.
# NOTE: this closure sets width as B = k*n*Q*S^(7/6), so B is PROPORTIONAL to Q,
# while the bed diffusivity goes as Q/B. Scaling both together leaves the response
# EXACTLY unchanged. The discharge increase therefore only does anything if the
# channel does not widen in proportion — so B is held fixed here, which is the
# defensible short-term case for cohesive banks. Over long enough for the banks to
# adjust, the effect would cancel.
Q_START, Q_END, Q_FACTOR = 1985, 1995, 1.15
DT_YR = 2

net, order, idx, m = make(sys.argv[1])
sol = DiffusionSolver(net, CLOSURE, intermittency=INTERMIT)
warnings.simplefilter("ignore")

# --- spin-up ---------------------------------------------------------------
for _ in range(5):
    sol.evolve(nt=40, dt=25 * YR, niter=3)
z_steady = [net.get_segment_field(s, "z").copy() for s in net.segment_ids]
sol._zold2 = None; sol._dt_prev = None          # restart the BDF2 history cleanly

z_bl0 = sol.z_bl
Q0 = [q.copy() for q in sol.Q]

def set_discharge(year):
    f = np.clip((year - Q_START) / float(Q_END - Q_START), 0.0, 1.0)
    scale = 1.0 + (Q_FACTOR - 1.0) * f
    sol.Q = [q * scale for q in Q0]          # B deliberately NOT rescaled
    return scale

def base_level(year):
    """2.5 m of rise, ramped linearly across the closure years."""
    f = np.clip((year - DAM_START) / float(DAM_END - DAM_START), 0.0, 1.0)
    return z_bl0 + DAM_RISE * f

def phase_of(year):
    if year < PULSE_END and year < DAM_START:      return "pulse"
    if DAM_START <= year <= DAM_END:               return "dam"
    if year < PULSE_END:                           return "pulse"
    if Q_START <= year <= Q_END:                   return "wetter"
    return "recovery"

frames = []
def snap(year):
    frames.append(dict(year=year, phase=phase_of(year), z_bl=sol.z_bl, qscale=qs,
                       z=[net.get_segment_field(s, "z").copy() for s in net.segment_ids]))

year = float(START_YEAR)
qs = set_discharge(year)
snap(year)
while year < END_YEAR:
    year += DT_YR
    sol.z_bl = base_level(year)
    qs = set_discharge(year)
    src = PULSE_M_PER_YR / YR if year <= PULSE_END else None
    sol.evolve(nt=1, dt=DT_YR * YR, niter=3, source_rate=src)
    snap(year)

out = dict(frames=frames, z_steady=z_steady, x=m["x"], Q=m["Q"], B=m["B"],
           z_dem=m["z_dem"], order=order, z_bl0=z_bl0,
           start_year=START_YEAR, pulse_end=PULSE_END, end_year=END_YEAR,
           dam=(DAM_START, DAM_END, DAM_RISE), qchange=(Q_START, Q_END, Q_FACTOR),
           dt_yr=DT_YR, pulse_rate=PULSE_M_PER_YR)
pickle.dump(out, open("/tmp/claude-1000/-home-awickert-models-FluvTree/4fd69dbf-d3b7-42b4-aeb8-7a9098b98b09/scratchpad/ww_frames.pkl", "wb"))

dz = [np.concatenate([a - b for a, b in zip(f["z"], z_steady)]) for f in frames]
yrs = np.array([f["year"] for f in frames])
def at(y):
    k = int(np.argmin(np.abs(yrs - y))); return dz[k]
print(f"{len(frames)} frames, {START_YEAR}-{END_YEAR}")
print(f"peak aggradation anywhere {max(float(np.max(d)) for d in dz):+.2f} m")
for y in (1932, 1938, 1940, 1984, 1996, 2100, END_YEAR):
    d = at(y)
    print(f"  {y}: mean {d.mean():+.2f} m   max {d.max():+.2f} m   "
          f"z_bl {frames[int(np.argmin(np.abs(yrs-y)))]['z_bl'] - z_bl0:+.2f} m   "
          f"Q x{frames[int(np.argmin(np.abs(yrs-y)))]['qscale']:.2f}")
