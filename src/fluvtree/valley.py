"""
Valley-storage geometry: the width primitive ``B(x, z)`` and the stored sediment
volume it implies. General to FluvTree -- a valley stores sediment regardless of
grain size, so this is *not* transport-specific and lives alongside the network
(referenced onto it), not on a closure.

The solver conserves stored sediment *volume* ``V``, not bed elevation ``z``; ``z``
is recovered through this geometry. The cross-section at each ``x`` is given as a set
of ``(z, B)`` pairs and read as a **piecewise-linear** width ``B(z)``. With ``B``
linear between pairs the volume ``V = (1 - lambda_p) * integral_0^z B dz'`` is
piecewise-*quadratic* and ``dV/dz = (1 - lambda_p) * B`` holds exactly -- the
consistency the volume-first solve needs.

The cross-sections are **ragged** (each ``x`` may carry a different number of
levels), stored as two flat value arrays ``z``, ``B`` sharing one ``offsets`` index
(``z`` and ``B`` are paired). Width and volume evaluate vectorized over all ``x``.

Outside the tabulated range the valley function should be *extended* (e.g. TerraPIN
updating it). As a backup when it has not been: below the lowest level the width
drops to one **channel width** (and warns); above the highest it holds the top width
constant (and warns).
"""

import warnings

import numpy as np


class ValleyGeometry(object):
    """
    Piecewise-linear, ragged valley cross-section ``B(x, z)`` over ``n`` nodes.

    Parameters
    ----------
    z_levels, B_levels : sequence of 1-D arrays, one per node
        The ``(z, B)`` cross-section at each node: elevations (ascending) and the
        valley width at each. Paired, so ``B_levels[i]`` matches ``z_levels[i]``.
        At least two levels per node (needed to interpolate).
    lambda_p : float
        Sediment porosity; storage carries the ``(1 - lambda_p)`` solid fraction.
    channel_width : float or array of length ``n``, optional
        Backup width used *below* the lowest tabulated level (a bed that has incised
        out of the valley sits in a channel-width slot). Kept as a stored value here,
        not pulled from a closure, so the geometry stays transport-agnostic. Required
        only if a query can fall below range.
    """

    def __init__(self, z_levels, B_levels, lambda_p, channel_width=None):
        counts = np.array([len(z) for z in z_levels])
        if np.any(counts < 2):
            raise ValueError("each node needs >= 2 (z, B) levels to interpolate")
        self.n = len(z_levels)
        self.lambda_p = float(lambda_p)
        self.offsets = np.concatenate(([0], np.cumsum(counts)))   # length n + 1
        self.z = np.concatenate([np.asarray(z, float) for z in z_levels])
        self.B = np.concatenate([np.asarray(b, float) for b in B_levels])
        if self.z.shape != self.B.shape:
            raise ValueError("z and B level counts must match at every node")
        # each node's levels must be sorted ascending in z
        for i in range(self.n):
            zi = self.z[self.offsets[i]:self.offsets[i + 1]]
            if np.any(np.diff(zi) <= 0):
                raise ValueError("node %d: z levels must be strictly ascending" % i)
        self.channel_width = channel_width
        self._counts = counts
        self._cumA = self._cumulative_area()   # solid-fraction-free area to each level

    # ------------------------------------------------------------------ #
    # Setup
    # ------------------------------------------------------------------ #

    def _cumulative_area(self):
        """Per-node cumulative cross-sectional area (integral of ``B``) from each
        node's lowest level up to each of its levels -- so a volume query needs only
        one partial trapezoid on top. Solid fraction ``(1 - lambda_p)`` is applied at
        read time, not here."""
        z, B = self.z, self.B
        # interval area from level l-1 to l; zero at each node's first level and
        # across node boundaries (no interval spans two nodes).
        ia = np.zeros_like(z)
        ia[1:] = 0.5 * (B[1:] + B[:-1]) * (z[1:] - z[:-1])
        ia[self.offsets[:-1]] = 0.0                      # node starts carry no area
        G = np.cumsum(ia)                                # global (un-reset) cumulative
        G_start = np.repeat(G[self.offsets[:-1]], self._counts)
        return G - G_start                               # reset to 0 at each node start

    def _locate(self, z_query):
        """For each node, the number ``k`` of levels at or below ``z_query`` (so the
        query sits in local interval ``[k-1, k]``); ``k == 0`` is below range,
        ``k == count`` is above range."""
        z_query = np.asarray(z_query, float)
        if z_query.shape != (self.n,):
            raise ValueError("z_query must have one value per node (length %d)"
                             % self.n)
        q_exp = np.repeat(z_query, self._counts)
        k = np.add.reduceat((self.z <= q_exp).astype(np.intp), self.offsets[:-1])
        # A query exactly at the top level is in range, not above it: keep it in the
        # last interior interval (interpolates to the top width) rather than tripping
        # the above-range backup.
        z_top = self.z[self.offsets[1:] - 1]
        at_top = (k == self._counts) & (z_query <= z_top)
        k = np.where(at_top, self._counts - 1, k)
        return z_query, k

    # ------------------------------------------------------------------ #
    # Geometry primitives (mirror grlp.Segment.valley_width / storage_*)
    # ------------------------------------------------------------------ #

    def valley_width(self, z_query):
        """Valley width ``B`` at bed elevation ``z_query`` (one value per node),
        piecewise-linearly interpolated. Below range -> channel width (warns); above
        range -> the top width held constant (warns)."""
        z_query, k = self._locate(z_query)
        counts = self._counts
        starts = self.offsets[:-1]
        below = k == 0
        above = k == counts
        interior = ~(below | above)

        out = np.empty(self.n, float)

        # interior: linear interpolation on interval [k-1, k]
        i = np.nonzero(interior)[0]
        g_lo = starts[i] + k[i] - 1
        g_hi = starts[i] + k[i]
        z_lo, z_hi = self.z[g_lo], self.z[g_hi]
        B_lo, B_hi = self.B[g_lo], self.B[g_hi]
        out[i] = B_lo + (B_hi - B_lo) * (z_query[i] - z_lo) / (z_hi - z_lo)

        # above range: hold the top width
        ia = np.nonzero(above)[0]
        out[ia] = self.B[starts[ia] + counts[ia] - 1]

        # below range: one channel width (backup for an un-extended valley function)
        ib = np.nonzero(below)[0]
        if ib.size:
            if self.channel_width is None:
                raise ValueError(
                    "bed below the tabulated valley at node(s) %s and no "
                    "channel_width backup was given" % ib.tolist())
            cw = np.broadcast_to(np.asarray(self.channel_width, float), (self.n,))
            out[ib] = cw[ib]

        self._warn_out_of_range(below, above)
        return out

    def storage_jacobian(self, z_query):
        """``dV/dz = (1 - lambda_p) * B(x, z)`` at the query elevation."""
        return (1.0 - self.lambda_p) * self.valley_width(z_query)

    def storage_volume(self, z_query):
        """Stored sediment volume per unit valley length,
        ``V = (1 - lambda_p) * integral B dz'`` up to ``z_query``. The datum is each
        node's lowest level (only differences of ``V`` enter the solve)."""
        z_query, k = self._locate(z_query)
        counts = self._counts
        starts = self.offsets[:-1]
        widths = self.valley_width(z_query)              # also handles the warnings
        below = k == 0
        above = k == counts
        interior = ~(below | above)

        area = np.empty(self.n, float)

        # interior: cumulative area to level k-1 + partial trapezoid up to the query
        i = np.nonzero(interior)[0]
        g_km1 = starts[i] + k[i] - 1
        z_km1, B_km1 = self.z[g_km1], self.B[g_km1]
        area[i] = self._cumA[g_km1] + 0.5 * (B_km1 + widths[i]) * (z_query[i] - z_km1)

        # above range: all intervals + a rectangle at the (constant) top width
        ia = np.nonzero(above)[0]
        g_top = starts[ia] + counts[ia] - 1
        area[ia] = self._cumA[g_top] + widths[ia] * (z_query[ia] - self.z[g_top])

        # below range: rectangle at the channel width, measured down from level 0
        ib = np.nonzero(below)[0]
        g_bot = starts[ib]
        area[ib] = widths[ib] * (z_query[ib] - self.z[g_bot])   # negative below datum

        return (1.0 - self.lambda_p) * area

    # ------------------------------------------------------------------ #
    # Internal
    # ------------------------------------------------------------------ #

    def _warn_out_of_range(self, below, above):
        nb, na = int(below.sum()), int(above.sum())
        if nb:
            warnings.warn(
                "valley width below the tabulated range at %d node(s): using the "
                "channel-width backup. The valley function should be extended "
                "(e.g. by TerraPIN)." % nb, RuntimeWarning)
        if na:
            warnings.warn(
                "valley width above the tabulated range at %d node(s): holding the "
                "top width constant. The valley function should be extended "
                "(e.g. by TerraPIN)." % na, RuntimeWarning)
