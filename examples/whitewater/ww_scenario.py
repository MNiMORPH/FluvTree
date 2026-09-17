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
PULSE_M_PER_YR = 0.045          # illustrative supply rate during the pulse
PULSE_YR, RECOVER_YR = 80, 420
DT_YR = 2

net, order, idx, m = make(sys.argv[1])
sol = DiffusionSolver(net, CLOSURE, intermittency=INTERMIT)
warnings.simplefilter("ignore")

# --- spin-up ---------------------------------------------------------------
for _ in range(5):
    sol.evolve(nt=40, dt=25 * YR, niter=3)
z_steady = [net.get_segment_field(s, "z").copy() for s in net.segment_ids]
sol._zold2 = None; sol._dt_prev = None          # restart the BDF2 history cleanly

frames = []
def snap(t, phase):
    frames.append(dict(t=t, phase=phase,
                       z=[net.get_segment_field(s, "z").copy() for s in net.segment_ids]))

snap(0.0, "pulse")
t = 0.0
for _ in range(PULSE_YR // DT_YR):
    sol.evolve(nt=1, dt=DT_YR * YR, niter=3,
               source_rate=PULSE_M_PER_YR / YR)
    t += DT_YR; snap(t, "pulse")
for _ in range(RECOVER_YR // DT_YR):
    sol.evolve(nt=1, dt=DT_YR * YR, niter=3)
    t += DT_YR; snap(t, "recovery")

out = dict(frames=frames, z_steady=z_steady, x=m["x"], Q=m["Q"], B=m["B"],
           z_dem=m["z_dem"], order=order,
           pulse_yr=PULSE_YR, dt_yr=DT_YR, pulse_rate=PULSE_M_PER_YR)
pickle.dump(out, open("/tmp/claude-1000/-home-awickert-models-FluvTree/4fd69dbf-d3b7-42b4-aeb8-7a9098b98b09/scratchpad/ww_frames.pkl", "wb"))

dz = [np.concatenate([a - b for a, b in zip(f["z"], z_steady)]) for f in frames]
peak = max(float(np.max(d)) for d in dz)
low  = min(float(np.min(d)) for d in dz)
iend = np.concatenate([a - b for a, b in zip(frames[PULSE_YR // DT_YR]["z"], z_steady)])
print(f"{len(frames)} frames, {t:.0f} yr")
print(f"peak aggradation {peak:+.2f} m   maximum incision below steady state {low:+.2f} m")
print(f"at the end of the pulse: mean {iend.mean():+.2f} m, max {iend.max():+.2f} m")
print(f"at the end:               mean {dz[-1].mean():+.2f} m, max {dz[-1].max():+.2f} m")
