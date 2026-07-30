"""
The plotting layer: the GRLP-style long profile, and the dependency invariant.

``fluvtree.plot`` is the only module that touches matplotlib. These tests pin that
``long_profile`` draws the network (a reach line + a junction/base-level connector
per segment) and that ``FluvTree.plot`` delegates to it -- and, crucially, that a
bare ``import fluvtree`` pulls in no matplotlib (the core stays dependency-light).
"""

import os
import subprocess
import sys

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")           # headless: no display needed

import fluvtree as ft
from fluvtree.plot import long_profile, slope_area, planform
from fluvtree.processes import build_network


_D = 2000.0


def _confluence():
    up, down = [[], [], [0, 1]], [[2], [2], []]
    x = [_D * np.arange(1, 5.0), _D * np.arange(1, 5.0), _D * np.arange(5, 9.0)]
    Q = [5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)]
    B = [100.0 * np.ones(4)] * 3
    z = [0.015 * (_D * 9 - xi) for xi in x]
    return build_network(up, down, x, z, Q, B, S0=0.015, x_bl=_D * 9, z_bl=0.0)


def test_long_profile_draws_reaches_and_connectors():
    net = _confluence()
    ax = long_profile(net)
    # each segment contributes its reach line + one connector (to the downstream
    # reach, or to base level at the mouth): two lines per segment.
    assert len(ax.lines) == 2 * len(net.segment_ids)
    assert ax.get_xlabel() == "Downstream distance [m]"
    assert ax.get_ylabel() == "Elevation [m]"


def test_long_profile_reaches_base_level_at_mouth():
    net = _confluence()
    ax = long_profile(net)
    # some drawn point sits exactly at base level (the mouth->z_bl connector)
    outlet = net.edge_of(net.mouth_segments()[0])[1]
    z_bl = net.get_node_field(outlet, "z_bl")
    ys = np.concatenate([ln.get_ydata() for ln in ax.lines])
    assert np.any(np.isclose(ys, z_bl))


def test_model_plot_delegates_to_long_profile():
    m = ft.FluvTree(_confluence()).add(ft.GravelBed(D=0.05))
    m.run(nt=50, dt=3.15e10)
    ax = m.plot(lw=2)
    assert ax.get_ylabel() == "Elevation [m]"
    assert len(ax.lines) == 2 * len(m.network.segment_ids)


def test_slope_area_draws_and_fits():
    net = _confluence()
    ax = slope_area(net)                          # defaults to slope vs Q
    assert ax.get_xlabel() == "Discharge [m$^3$/s]"
    assert ax.get_ylabel() == "Slope"
    assert len(ax.lines) == 2                      # scatter + the fit line
    assert len(slope_area(_confluence(), fit=False).lines) == 1
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"


def test_planform_lanes_reaches_and_connectors():
    net = _confluence()
    ax = planform(net)
    # a line per reach + a connector per non-mouth reach
    n = len(net.segment_ids)
    n_mouths = len(net.mouth_segments())
    assert len(ax.lines) == 2 * n - n_mouths
    assert ax.get_xlabel() == "Downstream distance [m]"
    # the two channel heads sit on distinct lanes; the confluence between them
    head_lanes = sorted({ax.lines[i].get_ydata()[0] for i in range(len(ax.lines))})
    assert len(head_lanes) >= 2


def test_model_plot_conveniences_delegate():
    m = ft.FluvTree(_confluence()).add(ft.GravelBed(D=0.05))
    m.run(nt=30, dt=3.15e10)
    assert m.plot_slope_area().get_ylabel() == "Slope"
    assert m.plot_planform().get_xlabel() == "Downstream distance [m]"


def test_import_fluvtree_is_matplotlib_free():
    # the core must not pull in matplotlib; only fluvtree.plot may.
    env = dict(os.environ)
    src = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
    env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
    code = ("import fluvtree, sys; "
            "mpl = [m for m in sys.modules if m.split('.')[0] == 'matplotlib']; "
            "assert not mpl, mpl")
    subprocess.run([sys.executable, "-c", code], check=True, env=env)
