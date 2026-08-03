"""
Random drainage-network topology by the Shreve (1974) algorithm.

A network of a given *magnitude* (number of source links) is drawn as a random
binary link string -- external (source) links and internal links -- and unfolded
into the upstream/downstream segment-ID lists FluvTree builds a network from. A
magnitude-``m`` network has ``m`` sources and ``m - 1`` confluences, i.e. ``2m - 1``
segments. Segment lengths and drainage areas are optional (constant, or drawn from
any object with a ``.rvs(size=...)`` method, e.g. a ``scipy.stats`` distribution).

Ported from GRLP's ``build_synthetic_network.Shreve_Random_Network`` (topology only;
the GRLP-hydrology field generators -- discharge, width, ssd, z -- are a separate
concern). Draws use a NumPy ``Generator`` so a ``seed`` makes a run reproducible.

Reference: Shreve, R. L. (1974), Water Resources Research.
"""

import numpy as np

from fluvtree.network import RiverNetwork


def _downstream_IDs(down_ls, i):
    """All segment IDs from ``i`` to the outlet (inclusive), from a raw
    downstream-ID list -- generation-internal, before a RiverNetwork exists."""
    IDs = []

    def _recurse(j):
        IDs.append(j)
        if down_ls[j]:
            _recurse(down_ls[j][0])

    _recurse(i)
    return IDs


def _upstream_IDs(up_ls, i):
    """All segment IDs at or above ``i``, from a raw upstream-ID list."""
    IDs = []

    def _recurse(j):
        IDs.append(j)
        for u in up_ls[j]:
            _recurse(u)

    _recurse(i)
    return IDs


class ShreveRandomNetwork(object):
    """
    A random drainage topology of a given magnitude (Shreve 1974).

    Parameters
    ----------
    magnitude : int
        Number of source (external) links -- the channel heads.
    segment_length : float or frozen distribution, optional
        A constant length for every segment, or an object with ``.rvs(size=1)``
        (e.g. a ``scipy.stats`` distribution) drawn per segment.
    max_length : float, optional
        If given, segment lengths are rescaled so the longest source-to-outlet path
        equals ``max_length``.
    segment_length_area_ratio : float or frozen distribution, optional
        Drainage area per unit length; with ``supply_area`` sets ``segment_areas``.
    supply_area : float or frozen distribution, optional
        Source-segment drainage area (head inflow area).
    topology : list of int, optional
        A prebuilt Shreve binary link string (external = 1, internal = 0) to unfold
        instead of drawing a new one.
    seed : int or numpy Generator, optional
        Seeds the topology draw for reproducibility.

    Attributes
    ----------
    upstream_segment_IDs, downstream_segment_IDs : list of list of int
        The per-segment topology (the input to ``RiverNetwork.from_segment_lists``).
    segment_lengths, segment_areas : list of float or None
        Set only if the corresponding parameters were given.
    """

    def __init__(self, magnitude=None, segment_length=None, max_length=None,
                 segment_length_area_ratio=None, supply_area=None,
                 segment_lengths=None, segment_length_area_ratios=None,
                 topology=None, seed=None):
        self.magnitude = magnitude
        self.links = topology
        self.segment_length = segment_length
        self.max_length = max_length
        self.segment_length_area_ratio = segment_length_area_ratio
        self.supply_area = supply_area
        self.segment_lengths = segment_lengths
        self.segment_length_area_ratios = segment_length_area_ratios
        self.source_areas = None
        self.segment_areas = None
        self._rng = (seed if isinstance(seed, np.random.Generator)
                     else np.random.default_rng(seed))

        if not self.links:
            self.build_network_topology()
        self.build_lists()
        if not segment_lengths and (segment_length or max_length):
            self.set_segment_lengths()
        if segment_length_area_ratio or segment_length_area_ratios:
            self.set_segment_areas()

    @property
    def N_external(self):
        """Number of external (source) links."""
        return sum(self.links)

    @property
    def N_internal(self):
        """Number of internal links."""
        return len(self.links) - sum(self.links)

    @property
    def probability_external(self):
        """Likelihood the next link is external (Shreve 1974)."""
        try:
            p = (self.N_internal - self.N_external)
            p *= (self.magnitude - self.N_external)
            p /= (self.N_internal - self.N_external + 1.0)
            p /= (2.0 * self.magnitude - 2.0 - self.N_internal - self.N_external)
        except ZeroDivisionError:
            p = np.inf
        return p

    def build_network_topology(self):
        """Draw the topology as a binary link string (external = 1, internal = 0),
        by the Shreve (1974) algorithm."""
        if self.magnitude is not None and self.magnitude < 2:
            self.links = []                 # a single segment (one source = outlet)
            return
        self.links = []
        k = 1
        last = 1
        while True:
            if self._rng.random() > self.probability_external:      # internal link
                self.links.append(0)
                if last:
                    k += 1
                else:
                    last = 0
            else:                                                   # external link
                self.links.append(1)
                if not last:
                    last = 1
                if last:
                    if k > 1:
                        k -= 1
                    else:
                        break
        self.links.remove(0)                                        # drop extra head

    def build_lists(self):
        """Unfold the link string into upstream/downstream segment-ID lists."""
        self.upstream_segment_IDs = [[]]
        self.downstream_segment_IDs = [[]]
        down_segs = []
        down_seg = 0
        seg = 1
        for i, l in enumerate(self.links):
            self.upstream_segment_IDs[down_seg].append(seg)
            self.upstream_segment_IDs.append([])
            self.downstream_segment_IDs.append([down_seg])
            if i < len(self.links) - 1:
                if not l:                       # internal link: continue upstream
                    down_segs.append(down_seg)
                    down_seg = seg
                else:                           # external: back down to a free link
                    while len(self.upstream_segment_IDs[down_seg]) > 1:
                        down_seg = down_segs.pop(-1)
                seg += 1

    def _draw(self, spec):
        """A constant ``spec``, or one draw from a ``.rvs``-bearing distribution."""
        try:
            return float(spec)
        except TypeError:
            return float(spec.rvs(size=1)[0])

    def set_segment_lengths(self):
        n = len(self.upstream_segment_IDs)
        if self.segment_length:
            self.segment_lengths = [self._draw(self.segment_length)
                                    for _ in range(n)]
        else:
            self.segment_lengths = list(np.ones(n))
        if self.max_length:
            lengths = [sum(self.segment_lengths[j]
                           for j in _downstream_IDs(self.downstream_segment_IDs, i))
                       for i in range(n)]
            scale = self.max_length / max(lengths)
            self.segment_lengths = [L * scale for L in self.segment_lengths]

    def set_segment_areas(self):
        n = len(self.upstream_segment_IDs)
        self.source_areas = [
            self._draw(self.supply_area) if not self.upstream_segment_IDs[i] else 0.0
            for i in range(n)]
        if not self.segment_length_area_ratios:
            self.segment_length_area_ratios = [
                self._draw(self.segment_length_area_ratio) for _ in range(n)]
        self.segment_areas = [
            self.segment_lengths[i] * self.segment_length_area_ratios[i]
            for i in range(n)]


def random_network(magnitude, seed=None, **kwargs):
    """
    A random drainage network of the given magnitude, as a :class:`RiverNetwork`.

    Draws a Shreve (1974) topology and returns it as a *topology-only* network (via
    :meth:`RiverNetwork.from_segment_lists`); attach ``x, z, Q, B`` with
    :meth:`RiverNetwork.from_arrays` or the planned field generators to make it
    runnable. ``kwargs`` are passed to :class:`ShreveRandomNetwork` (e.g.
    ``segment_length``, ``max_length``).
    """
    gen = ShreveRandomNetwork(magnitude=magnitude, seed=seed, **kwargs)
    return RiverNetwork.from_segment_lists(gen.upstream_segment_IDs,
                                           gen.downstream_segment_IDs)
