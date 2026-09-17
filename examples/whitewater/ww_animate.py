"""Animate the sand-bed lower Whitewater through a sediment pulse and recovery."""
import json, pickle, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

SP = "/tmp/claude-1000/-home-awickert-models-FluvTree/4fd69dbf-d3b7-42b4-aeb8-7a9098b98b09/scratchpad"
D = pickle.load(open(f"{SP}/ww_frames.pkl", "rb"))
J = json.load(open(f"{SP}/wwsand2m.json"))
frames, zs, X, order = D["frames"], D["z_steady"], D["x"], D["order"]
DAM_A, DAM_B, DAM_RISE = D["dam"]
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
AGG, POOL = "#b07d3a", "#2d6e8e"         # sediment stored / impounded water
dzmax = max(float(np.max(np.concatenate([a - b for a, b in zip(f["z"], zs)])))
            for f in frames)
mean_series = [float(np.mean(np.concatenate([a - b for a, b in zip(f["z"], zs)])))
               for f in frames]
tt = np.array([f["year"] for f in frames])

fig = plt.figure(figsize=(12.8, 7.2), facecolor="white")
gs = fig.add_gridspec(2, 2, width_ratios=[2.25, 1], height_ratios=[2.5, 1],
                      left=0.058, right=0.995, top=0.865, bottom=0.085,
                      wspace=0.13, hspace=0.42)
axP = fig.add_subplot(gs[0, 0])       # long profile
axT = fig.add_subplot(gs[1, 0])       # time series
axM = fig.add_subplot(gs[:, 1])       # map, full height

# --- static: long profile reference ---------------------------------------
for i in range(len(order)):
    axP.plot(X[i] / 1000, zs[i], lw=1.1, ls=(0, (4, 3)), color=MUTED, zorder=2)
allz = np.concatenate(zs)
axP.set_xlim(min(x.min() for x in X) / 1000 - 0.4, 0.6)
axP.set_ylim(allz.min() - 2.0, allz.max() + 5.0)
axP.set_xlabel("distance from the mouth (km)", fontsize=10, color=MUTED)
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

# --- static: map -----------------------------------------------------------
for i in range(len(order)):
    axM.plot(PX[i] / 1000, PY[i] / 1000, lw=5.0, color="#ececec", zorder=1)
axM.set_aspect("equal")
axM.set_title("sediment stored, by reach", fontsize=10.5, color=INK, loc="left", pad=6)
axM.tick_params(labelsize=8, colors=MUTED)
axM.set_xticklabels([])
axM.set_ylabel("northing (km)", fontsize=8.5, color=MUTED)
for sp in axM.spines.values(): sp.set_color("#d0d0d0")
segsM = [axM.scatter(PX[i] / 1000, PY[i] / 1000, s=16, c=np.zeros_like(PX[i]),
                     cmap="YlOrBr", vmin=0, vmax=dzmax, zorder=3, linewidths=0)
         for i in range(len(order))]
cb = fig.colorbar(segsM[0], ax=axM, fraction=0.055, pad=0.03,
                  orientation="horizontal", location="bottom")
cb.set_label("metres of sediment", fontsize=9, color=MUTED)
cb.ax.tick_params(labelsize=8, colors=MUTED); cb.outline.set_visible(False)

# --- static: time series ---------------------------------------------------
axT.axvspan(tt.min(), PULSE_END, color=AGG, alpha=0.15, lw=0)
axT.axvspan(DAM_A, DAM_B, color=POOL, alpha=0.30, lw=0)
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
axT.annotate("lock and dam", xy=(DAM_B, max(mean_series) * 0.40),
             xytext=(DAM_B + 70, max(mean_series) * 0.30), fontsize=9.5, color=POOL,
             arrowprops=dict(arrowstyle="->", color=POOL, lw=1.1))
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
        segsM[i].set_array(dz[i])
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

    ph = {"pulse": "sediment pouring in",
          "dam": "lock and dam closing",
          "recovery": "the river cannot clear it"}[f["phase"]]
    title.set_text(f"Lower Whitewater, sand bed — {ph}")
    sub.set_text(f"{f['year']:.0f}    ·    mean {mean_series[k]:+.2f} m stored"
                 f"    ·    peak {max(float(np.max(d)) for d in dz):+.2f} m")
    return lines + segsM + [marker]

ani = FuncAnimation(fig, update, frames=len(frames), interval=70, blit=False)
out = sys.argv[1] if len(sys.argv) > 1 else f"{SP}/whitewater_sand.mp4"
ani.save(out, writer=FFMpegWriter(fps=14, bitrate=2800), dpi=110)
print("wrote", out, f"({len(frames)} frames)")
