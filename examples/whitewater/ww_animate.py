"""Animate the sand-bed lower Whitewater through a sediment pulse and recovery."""
import json, pickle, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.collections import LineCollection

SP = "/tmp/claude-1000/-home-awickert-models-FluvTree/4fd69dbf-d3b7-42b4-aeb8-7a9098b98b09/scratchpad"
D = pickle.load(open(f"{SP}/ww_frames.pkl", "rb"))
NET = sys.argv[2] if len(sys.argv) > 2 else f"{SP}/wwsand3.json"
J = json.load(open(NET))
frames, zs, X, order = D["frames"], D["z_steady"], D["x"], D["order"]
DAM_A, DAM_B, DAM_RISE = D["dam"]
QA, QB, QF = D["qchange"]
PULSE_END, Z_BL0 = D["pulse_end"], D["z_bl0"]

# planform, resampled to the model nodes
link = {l["source"]: l for l in J["links"]}
PX, PY = [], []
for i, c in enumerate(order):
    L = link[c]
    s = np.asarray(L["s_downstream"], float); s -= s[0]
    sn = np.linspace(s[0], s[-1], len(X[i]))
    PX.append(np.interp(sn, s, np.asarray(L["x"], float)))
    PY.append(np.interp(sn, s, np.asarray(L["y"], float)))

INK, MUTED = "#1a1a1a", "#7a7a7a"
AGG, POOL, FLOW = "#b07d3a", "#2d6e8e", "#4f7f52"   # sediment / pool / discharge
dzmax = max(float(np.max(np.concatenate([a - b for a, b in zip(f["z"], zs)])))
            for f in frames)
mean_series = [float(np.mean(np.concatenate([a - b for a, b in zip(f["z"], zs)])))
               for f in frames]
tt = np.array([f["year"] for f in frames])

fig = plt.figure(figsize=(12.8, 7.2), facecolor="white")
gs = fig.add_gridspec(3, 2, width_ratios=[2.25, 1], height_ratios=[1.45, 1.25, 0.85],
                      left=0.058, right=0.995, top=0.865, bottom=0.075,
                      wspace=0.13, hspace=0.50)
axP = fig.add_subplot(gs[0, 0])       # long profile
axD = fig.add_subplot(gs[1, 0])       # deviation from steady state
axT = fig.add_subplot(gs[2, 0])       # time series
axM = fig.add_subplot(gs[:, 1])       # map, full height

# --- static: long profile reference ---------------------------------------
for i in range(len(order)):
    axP.plot(X[i] / 1000, zs[i], lw=1.1, ls=(0, (4, 3)), color=MUTED, zorder=2)
allz = np.concatenate(zs)
axP.set_xlim(min(x.min() for x in X) / 1000 - 0.4, 0.6)
axP.set_ylim(allz.min() - 2.0, allz.max() + 5.0)
axP.tick_params(labelbottom=False)
pool_line, = axP.plot([], [], lw=1.6, color=POOL, zorder=5)
pool_fill = [None]
dam_note = axP.text(0.42, 0.055, "", transform=axP.transAxes, fontsize=9.5,
                    color=POOL, ha="left")
axP.set_ylabel("bed elevation (m)", fontsize=10, color=MUTED)
axP.tick_params(labelsize=9, colors=MUTED)
axP.grid(axis="y", color="#ededed", lw=0.8); axP.set_axisbelow(True)
for sp in ("top", "right"): axP.spines[sp].set_visible(False)
for sp in ("left", "bottom"): axP.spines[sp].set_color("#d0d0d0")
axP.plot([], [], lw=1.1, ls=(0, (4, 3)), color=MUTED, label="pre-disturbance steady state")
axP.plot([], [], lw=2.2, color=INK, label="bed now")
axP.fill_between([], [], color=AGG, alpha=0.55, label="sediment stored in the valley")
axP.legend(loc="upper right", frameon=False, fontsize=9.5, labelcolor=INK)

lines = [axP.plot([], [], lw=2.2, color=INK, solid_capstyle="round", zorder=4)[0]
         for _ in order]
fills = [None] * len(order)

# --- static: deviation from steady state ----------------------------------
axD.axhline(0, lw=1.0, color=MUTED, zorder=1)
dlines = [axD.plot([], [], lw=2.0, color=AGG, solid_capstyle="round", zorder=3)[0]
          for _ in order]
dfills = [None] * len(order)
axD.set_xlim(axP.get_xlim())
axD.set_ylim(-0.3, dzmax * 1.12)
axD.set_ylabel("sediment stored (m)", fontsize=9.5, color=MUTED)
axD.set_xlabel("distance from the mouth (km)", fontsize=9.5, color=MUTED)
axD.tick_params(labelsize=8.5, colors=MUTED)
axD.grid(axis="y", color="#ededed", lw=0.8); axD.set_axisbelow(True)
for sp in ("top", "right"): axD.spines[sp].set_visible(False)
for sp in ("left", "bottom"): axD.spines[sp].set_color("#d0d0d0")
axD.set_title("thickness of fill, reach by reach", fontsize=10.5, color=INK,
              loc="left", pad=5)

# --- static: map -----------------------------------------------------------
axM.set_aspect("equal")
axM.set_title("sediment stored, by reach", fontsize=10.5, color=INK, loc="left", pad=6)
axM.tick_params(labelsize=8, colors=MUTED)
axM.set_xticklabels([])
axM.set_ylabel("northing (km)", fontsize=8.5, color=MUTED)
for sp in axM.spines.values(): sp.set_color("#d0d0d0")
segsM = []
for i in range(len(order)):
    pts = np.column_stack([PX[i] / 1000, PY[i] / 1000]).reshape(-1, 1, 2)
    lc = LineCollection(np.concatenate([pts[:-1], pts[1:]], axis=1),
                        cmap="YlOrBr", norm=plt.Normalize(0, dzmax),
                        linewidths=4.0, capstyle="round", joinstyle="round", zorder=3)
    lc.set_array(np.zeros(len(PX[i]) - 1))
    axM.add_collection(lc)
    segsM.append(lc)
axM.set_xlim(min(p.min() for p in PX) / 1000 - 0.4, max(p.max() for p in PX) / 1000 + 0.4)
axM.set_ylim(min(p.min() for p in PY) / 1000 - 0.4, max(p.max() for p in PY) / 1000 + 0.4)
cb = fig.colorbar(segsM[0], ax=axM, fraction=0.055, pad=0.03,
                  orientation="horizontal", location="bottom")
cb.set_label("metres of sediment", fontsize=9, color=MUTED)
cb.ax.tick_params(labelsize=8, colors=MUTED); cb.outline.set_visible(False)

# --- static: time series ---------------------------------------------------
axT.axvspan(tt.min(), PULSE_END, color=AGG, alpha=0.15, lw=0)
axT.axvspan(DAM_A, DAM_B, color=POOL, alpha=0.30, lw=0)
axT.axvspan(QA, QB, color=FLOW, alpha=0.26, lw=0)          # the ramp
axT.axvspan(QB, tt.max(), color=FLOW, alpha=0.11, lw=0)    # held, not ended
axT.plot(tt, mean_series, lw=2.0, color=AGG)
axT.set_xlim(tt.min(), tt.max()); axT.set_ylim(-0.12, max(mean_series) * 1.12)
axT.set_xlabel("year", fontsize=9, color=MUTED)
axT.set_ylabel("mean sediment stored (m)", fontsize=9, color=MUTED)
axT.tick_params(labelsize=8, colors=MUTED)
axT.grid(axis="y", color="#ededed", lw=0.8); axT.set_axisbelow(True)
for sp in ("top", "right"): axT.spines[sp].set_visible(False)
for sp in ("left", "bottom"): axT.spines[sp].set_color("#d0d0d0")
axT.text(tt.min() + 8, max(mean_series) * 0.93, "supply pulse", ha="left",
         fontsize=9.5, color=AGG)
axT.annotate("lock and dam", xy=(DAM_B, max(mean_series) * 0.42),
             xytext=(DAM_B + 62, max(mean_series) * 0.26), fontsize=9.5, color=POOL,
             arrowprops=dict(arrowstyle="->", color=POOL, lw=1.1))
axT.annotate(f"discharge to +{100*(QF-1):.0f}%, then held",
             xy=(QB + 30, max(mean_series) * 0.28),
             xytext=(QB + 60, max(mean_series) * 0.08), fontsize=9.5, color=FLOW,
             arrowprops=dict(arrowstyle="->", color=FLOW, lw=1.1))
marker, = axT.plot([], [], "o", ms=7, color=INK, zorder=5)

title = fig.text(0.058, 0.945, "", ha="left", fontsize=15.5, color=INK)
sub = fig.text(0.058, 0.902, "", ha="left", fontsize=10.5, color=MUTED)

def update(k):
    f = frames[k]
    dz = [a - b for a, b in zip(f["z"], zs)]
    for i in range(len(order)):
        lines[i].set_data(X[i] / 1000, f["z"][i])
        if fills[i] is not None:
            fills[i].remove()
        fills[i] = axP.fill_between(X[i] / 1000, zs[i], f["z"][i],
                                    color=AGG, alpha=0.55, lw=0, zorder=3)
        segsM[i].set_array(0.5 * (dz[i][:-1] + dz[i][1:]))
        dlines[i].set_data(X[i] / 1000, dz[i])
        if dfills[i] is not None:
            dfills[i].remove()
        dfills[i] = axD.fill_between(X[i] / 1000, 0, dz[i], color=AGG,
                                     alpha=0.35, lw=0, zorder=2)
    marker.set_data([f["year"]], [mean_series[k]])

    # the navigation pool at the mouth
    zb = f["z_bl"]
    xr = axP.get_xlim()[1]
    pool_line.set_data([0.0, xr], [zb, zb])
    if pool_fill[0] is not None:
        pool_fill[0].remove()
    pool_fill[0] = axP.fill_between([0.0, xr], [zb, zb], axP.get_ylim()[0],
                                    color=POOL, alpha=0.22, lw=0, zorder=1)
    rise = zb - Z_BL0
    dam_note.set_text(f"navigation pool  ·  base level {rise:+.2f} m" if rise > 0.01 else "")

    ph = {"pulse": "land clearance upstream",
          "dam": "lock and dam closing",
          "wetter": f"discharge rising {100*(QF-1):.0f}%",
          "recovery": "the river cannot clear it"}[f["phase"]]
    title.set_text("Whitewater River watershed")
    q = f.get("qscale", 1.0)
    sub.set_text(f"{f['year']:.0f}   ·   {ph}   ·   mean {mean_series[k]:+.2f} m stored"
                 f"   ·   peak {max(float(np.max(d)) for d in dz):+.2f} m"
                 + (f"   ·   Q x{q:.2f}" if q > 1.001 else ""))
    return lines + dlines + segsM + [marker]

ani = FuncAnimation(fig, update, frames=len(frames), interval=70, blit=False)
out = sys.argv[1] if len(sys.argv) > 1 else f"{SP}/whitewater_sand.mp4"
ani.save(out, writer=FFMpegWriter(fps=14, bitrate=2800), dpi=110)
print("wrote", out, f"({len(frames)} frames)")
