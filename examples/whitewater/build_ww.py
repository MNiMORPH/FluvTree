"""Build a FluvTree RiverNetwork from the sand-bed lower Whitewater."""
import json, numpy as np

DX = 100.0          # node spacing [m]

def build(path, dx=DX, clip=None):
    """clip: {cat: km} — keep only the lowest `km` of that link."""
    clip = clip or {}
    d = json.load(open(path))
    links = {l["source"]: l for l in d["links"]}
    to = {l["source"]: l["target"] for l in d["links"]}
    up = {}
    for s, t in to.items():
        up.setdefault(t, []).append(s)
    order = sorted(links)
    idx = {c: i for i, c in enumerate(order)}

    # distance from the outlet, walking downstream, so x increases downstream
    def dist_to_outlet(c):
        d0 = 0.0
        while to[c] != 0:
            L = links[to[c]]
            s = np.asarray(L["s_downstream"], float)
            d0 += s[-1] - s[0]
            c = to[c]
        return d0

    segs = {}
    for c in order:
        L = links[c]
        s = np.asarray(L["s_downstream"], float); s = s - s[0]
        keep = clip.get(c)
        if keep is not None:                       # retain only the lowest `keep` m
            m = s >= (s[-1] - keep)
            L = dict(L)
            for k in ("s_downstream", "z", "A", "x", "y"):
                L[k] = list(np.asarray(L[k], float)[m])
            s = np.asarray(L["s_downstream"], float); s = s - s[0]
        z = np.asarray(L["z"], float)
        A = np.asarray(L["A"], float)
        n = max(3, int(round((s[-1] - s[0]) / dx)) + 1)
        sn = np.linspace(s[0], s[-1], n)
        zn = np.interp(sn, s, z)
        An = np.interp(sn, s, np.abs(A))
        d_out_end = dist_to_outlet(c)            # distance from this link's end to outlet
        x = -(d_out_end + (s[-1] - sn))          # negative upstream, 0 at outlet
        # planform of the CLIPPED link, on the model's nodes — the map must not
        # draw the full-length line stretched onto fewer nodes
        pxn = np.interp(sn, s, np.asarray(L["x"], float))
        pyn = np.interp(sn, s, np.asarray(L["y"], float))
        segs[c] = dict(x=x, z=zn, A=An, px=pxn, py=pyn)

    # The junction is a FACE, not a shared node: give the junction vertex to the
    # DOWNSTREAM link and end each upstream link one spacing earlier. Two nodes at
    # one coordinate give a zero-length face -> exactly singular matrix, z all NaN,
    # warning only. (See LandGraph docs/SPEC-network-file.md §1.6.)
    for c in order:
        if to[c] != 0:
            for k in ("x", "z", "A", "px", "py"):
                segs[c][k] = segs[c][k][:-1]
    return order, idx, to, up, segs

if __name__ == "__main__":
    import sys
    order, idx, to, up, segs = build(sys.argv[1])
    for c in order:
        S = segs[c]
        print(f"link {c:>3}: {len(S['x']):>4} nodes  x {S['x'][0]:9.0f}..{S['x'][-1]:9.0f}  "
              f"z {S['z'][0]:6.2f}..{S['z'][-1]:6.2f}  A {S['A'][0]/1e6:6.1f}..{S['A'][-1]/1e6:6.1f} km2")
    print("\ntopology:", {c: to[c] for c in order})
    print("upstream:", {c: up.get(c, []) for c in order})
