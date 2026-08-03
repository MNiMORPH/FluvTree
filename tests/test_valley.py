"""
Valley-storage geometry: piecewise-linear, ragged ``B(x, z)`` cross-sections.

Checks the width interpolation, the stored-volume integral, and the one property
the volume-first solve depends on -- ``dV/dz == storage_jacobian`` -- plus the
out-of-range backups (below -> channel width, above -> held top width) and the
construction guards.
"""

import warnings

import numpy as np
import pytest

from fluvtree.valley import ValleyGeometry

_LAM = 0.35


def test_rectangular_reproduces_constant_B():
    # constant width -> V = (1 - lambda_p) * B * z, the current constant-B model
    vg = ValleyGeometry([np.array([0.0, 50.0])], [np.array([100.0, 100.0])], _LAM)
    z = np.array([30.0])
    assert np.allclose(vg.valley_width(z), 100.0)
    assert np.allclose(vg.storage_volume(z), (1 - _LAM) * 100.0 * 30.0)


def test_trapezoid_width_and_volume_analytic():
    # B(z) = 100 + 2 m z, so V = (1 - lambda_p) * (100 z + m z^2)
    m = 0.5
    vg = ValleyGeometry([np.array([0.0, 50.0])],
                        [np.array([100.0, 100.0 + 2 * m * 50])], _LAM)
    for zq in (10.0, 25.0, 50.0):
        q = np.array([zq])
        assert np.allclose(vg.valley_width(q), 100.0 + 2 * m * zq)
        assert np.allclose(vg.storage_volume(q),
                           (1 - _LAM) * (100.0 * zq + m * zq ** 2))


def test_jacobian_is_dVdz():
    # the consistency the volume-first solve needs: storage_jacobian == d/dz storage_volume
    m = 0.5
    vg = ValleyGeometry([np.array([0.0, 50.0])],
                        [np.array([100.0, 100.0 + 2 * m * 50])], _LAM)
    z = np.array([27.3])
    h = 1e-4
    fd = (vg.storage_volume(z + h) - vg.storage_volume(z - h)) / (2 * h)
    assert np.allclose(fd, vg.storage_jacobian(z), rtol=1e-6)


def test_ragged_multiple_nodes():
    # node 0 has 2 levels, node 1 has 3 -- different counts, one call
    vg = ValleyGeometry(
        [np.array([0.0, 40.0]), np.array([0.0, 20.0, 60.0])],
        [np.array([80.0, 120.0]), np.array([50.0, 90.0, 200.0])], _LAM)
    zq = np.array([20.0, 20.0])
    assert np.allclose(vg.valley_width(zq), [100.0, 90.0])
    assert np.allclose(vg.storage_volume(zq),
                       [(1 - _LAM) * (80.0 * 20 + 0.5 * (120 - 80) / 40 * 20 ** 2),
                        (1 - _LAM) * 0.5 * (50.0 + 90.0) * 20])


def test_below_range_uses_channel_width_and_warns():
    vg = ValleyGeometry([np.array([0.0, 50.0])], [np.array([100.0, 150.0])],
                        _LAM, channel_width=8.0)
    with pytest.warns(RuntimeWarning, match="below the tabulated"):
        assert np.allclose(vg.valley_width(np.array([-5.0])), 8.0)


def test_above_range_holds_top_and_warns():
    vg = ValleyGeometry([np.array([0.0, 50.0])], [np.array([100.0, 150.0])], _LAM)
    with pytest.warns(RuntimeWarning, match="above the tabulated"):
        assert np.allclose(vg.valley_width(np.array([70.0])), 150.0)


def test_boundary_levels_do_not_warn():
    # exactly at the lowest or highest tabulated level is in range, not out of it
    vg = ValleyGeometry([np.array([0.0, 50.0])], [np.array([100.0, 150.0])], _LAM)
    with warnings.catch_warnings():
        warnings.simplefilter("error")   # any warning fails the test
        assert np.allclose(vg.valley_width(np.array([0.0])), 100.0)
        assert np.allclose(vg.valley_width(np.array([50.0])), 150.0)


def test_below_range_without_backup_raises():
    vg = ValleyGeometry([np.array([0.0, 50.0])], [np.array([100.0, 150.0])], _LAM)
    with pytest.raises(ValueError, match="channel_width"):
        vg.valley_width(np.array([-1.0]))


def test_construction_guards():
    with pytest.raises(ValueError, match=">= 2"):
        ValleyGeometry([np.array([0.0])], [np.array([100.0])], _LAM)
    with pytest.raises(ValueError, match="ascending"):
        ValleyGeometry([np.array([50.0, 0.0])], [np.array([100.0, 150.0])], _LAM)
