"""
Validate the transport closures against the models they are extracted from.

Gravel is checked against a live ``grlp`` (pip-installed, authoritative); sand
against a live ``srlp`` loaded from the sibling checkout if present. The closure
must reproduce each model's ``k_Qs``, flux exponent, and width/depth relations
exactly -- this is what makes the extraction a re-homing, not a re-derivation.
"""

import importlib.util
import os

import numpy as np
import pytest

from fluvtree.closures.base import TransportClosure
from fluvtree.closures.gravel import GravelClosure
from fluvtree.closures.sand import SandClosure


def _load_srlp():
    """Load the srlp module from the sibling checkout, or skip."""
    path = os.path.expanduser("~/models/SRLP/srlp/srlp.py")
    if not os.path.exists(path):
        pytest.skip("SRLP checkout not present")
    spec = importlib.util.spec_from_file_location("srlp", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# Shared flux form
# --------------------------------------------------------------------------- #

def test_shared_flux_form_and_exponents():
    gravel, sand = GravelClosure(), SandClosure(D=0.3e-3, n=0.03, tau_crit_bank=5.0)
    assert gravel.p == pytest.approx(7 / 6.0)
    assert sand.p == pytest.approx(5 / 6.0)
    # the implicit conductance goes as S**(p-1)
    assert gravel.conductance_exponent == pytest.approx(1 / 6.0)
    assert sand.conductance_exponent == pytest.approx(-1 / 6.0)


def test_sediment_discharge_signed_form():
    c = GravelClosure()
    Q, S = 10.0, 0.01
    # downhill (positive S) -> positive Q_s = k_Qs * I * Q * S**p
    assert c.sediment_discharge(Q, S) == pytest.approx(c.k_Qs * Q * S ** c.p)
    # sign follows the slope
    assert c.sediment_discharge(Q, -S) == pytest.approx(-c.k_Qs * Q * S ** c.p)
    # intermittency scales linearly
    assert c.sediment_discharge(Q, S, 0.5) == pytest.approx(0.5 * c.sediment_discharge(Q, S))


def test_base_class_requires_implementation():
    with pytest.raises(NotImplementedError):
        TransportClosure().k_Qs


# --------------------------------------------------------------------------- #
# Gravel vs GRLP (authoritative)
# --------------------------------------------------------------------------- #

def test_gravel_k_Qs_matches_grlp():
    grlp = pytest.importorskip("grlp")
    seg = grlp.LongProfile()
    seg.basic_constants()
    seg.bedload_lumped_constants()
    assert GravelClosure().k_Qs == pytest.approx(seg.k_Qs, rel=0, abs=0)


def test_gravel_width_and_depth_match_grlp():
    grlp = pytest.importorskip("grlp")
    D = 0.05
    Q = np.array([10.0, 20.0])
    S = np.array([0.01, 0.008])
    seg = grlp.LongProfile()
    seg.basic_constants()
    seg.bedload_lumped_constants()
    seg.D, seg.Q, seg.S = D, Q, S
    seg.compute_channel_width()
    seg.compute_flow_depth()
    c = GravelClosure(D=D)
    assert np.allclose(c.channel_width(Q, S), seg.b, rtol=0, atol=1e-12)
    assert np.allclose(c.channel_depth(S), seg.h, rtol=0, atol=1e-12)


# --------------------------------------------------------------------------- #
# Sand vs SRLP
# --------------------------------------------------------------------------- #

def test_sand_k_Qs_matches_srlp():
    srlp = _load_srlp()
    D, n, tau = 0.3e-3, 0.03, 5.0
    seg = srlp.LongProfile()
    seg.basic_constants()
    seg.set_D(D)
    seg.set_Mannings_roughness(n)
    seg.set_tau_crit_bank(tau)
    seg.sediment_lumped_constants()
    assert SandClosure(D=D, n=n, tau_crit_bank=tau).k_Qs == pytest.approx(seg.k_Qs, rel=0, abs=0)


def test_sand_width_and_depth_match_srlp():
    srlp = _load_srlp()
    D, n, tau = 0.3e-3, 0.03, 5.0
    Q = np.array([50.0, 80.0])
    S = np.array([2e-4, 1.5e-4])
    seg = srlp.LongProfile()
    seg.basic_constants()
    seg.set_D(D)
    seg.set_Mannings_roughness(n)
    seg.set_tau_crit_bank(tau)
    seg.Q, seg.S = Q, S
    seg.compute_channel_width()
    seg.compute_flow_depth()
    c = SandClosure(D=D, n=n, tau_crit_bank=tau)
    assert np.allclose(c.channel_width(Q, S), seg.b, rtol=0, atol=1e-12)
    assert np.allclose(c.channel_depth(S), seg.h, rtol=0, atol=1e-12)
