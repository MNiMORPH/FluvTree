"""
Benchmark the FluvTree GRLP process against stable GRLP (ground truth).

The adapter reuses GRLP's solver numerics, so driving the same network from the
same initial conditions through the FluvTree graph must reproduce a standalone
``grlp.Network`` -- to machine precision. This regression protects the graph
<-> solver state round-trip (segment ordering, head-slope ordering, boundary
placement): a plumbing bug there would show up as drift from ground truth.

grlp is a hard dependency of this process; skipped if unavailable.
"""

import numpy as np
import pytest

grlp = pytest.importorskip("grlp")

from fluvtree.processes import GRLPProcess, build_grlp_network


_D = 2000.0

# topology + Q (x arrays built per segment length below)
SPECS = {
    "symmetric_confluence": dict(
        x=[_D * np.arange(1, 5.0), _D * np.arange(1, 5.0), _D * np.arange(5, 9.0)],
        Q=[5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=_D * 9,
    ),
    "multi_level": dict(
        x=[_D * np.arange(1, 5.0), _D * np.arange(1, 5.0),
           _D * np.arange(1, 5.0), _D * np.arange(1, 5.0),
           _D * np.arange(5, 9.0), _D * np.arange(9, 13.0),
           _D * np.arange(13, 17.0)],
        Q=[5 * np.ones(4)] * 4 + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
        up=[[], [], [], [], [0, 1], [4, 2], [5, 3]],
        down=[[4], [4], [5], [6], [5], [6], []], x_bl=_D * 17,
    ),
}

S0 = 0.015
B = 100.0
Z_BL = 0.0


def _initial_z(spec):
    """A fixed, non-flat initial bed for each segment (descends downstream)."""
    return [np.linspace(60.0, 20.0, len(xi)) for xi in spec["x"]]


def _apply_euler(net):
    net.set_time_integration(1)
    net.set_iteration_tolerance(None)
    net.set_niter(3)


def _apply_bdf2(net):
    net.set_time_integration(2)
    net.set_iteration_tolerance(1.0e-4)


SCHEMES = {"euler": _apply_euler, "bdf2": _apply_bdf2}


def _ground_truth(spec, z0, configure):
    n = len(spec["x"])
    n_heads = sum(1 for i in range(n) if len(spec["up"][i]) == 0)
    net = grlp.Network()
    net.initialize(
        x_bl=spec["x_bl"], z_bl=Z_BL, S0=[S0] * n_heads, Q_s_0=None,
        upstream_segment_IDs=spec["up"], downstream_segment_IDs=spec["down"],
        x=[xi.copy() for xi in spec["x"]], z=[zi.copy() for zi in z0],
        Q=[qi.copy() for qi in spec["Q"]],
        B=[B * np.ones(len(xi)) for xi in spec["x"]],
    )
    configure(net)
    net.get_z_lengths()
    return net


@pytest.mark.parametrize("scheme", list(SCHEMES))
@pytest.mark.parametrize("name", list(SPECS))
def test_process_matches_standalone_grlp(name, scheme):
    spec = SPECS[name]
    z0 = _initial_z(spec)
    configure = SCHEMES[scheme]

    gt = _ground_truth(spec, z0, configure)
    rn = build_grlp_network(
        spec["up"], spec["down"], spec["x"], z0, spec["Q"],
        [B * np.ones(len(xi)) for xi in spec["x"]], S0, spec["x_bl"], Z_BL)
    proc = GRLPProcess(rn, configure=configure)

    dt = 3.15e10
    for _ in range(6):
        gt.evolve_threshold_width_river_network(nt=1, dt=dt)
        proc.step(dt, nt=1)

    for i, s in enumerate(proc._segs):
        z_gt = gt.segments[i].z
        z_fluv = rn.get_segment_field(s, "z")
        assert np.allclose(z_gt, z_fluv, rtol=0, atol=1e-9), \
            f"{name}/{scheme} seg {s}: max|dz|={np.max(np.abs(z_gt - z_fluv)):.3e}"
