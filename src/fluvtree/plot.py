"""
Views of a river network -- the plotting layer.

This is the *only* module that imports matplotlib. The numerics core
(``network`` / ``solvers`` / ``closures``) stays dependency-light and never
imports this module; the arrow points one way (``plot`` reads the core), so bare
``import fluvtree`` pulls in no plotting. Plot functions take a
:class:`~fluvtree.network.RiverNetwork` and draw onto a matplotlib ``Axes``;
:meth:`fluvtree.model.FluvTree.plot` is the convenience that calls
:func:`long_profile` on its network.

Three views: :func:`long_profile` (elevation vs downstream distance, reaches joined
across confluences -- the canonical GRLP plot), :func:`slope_area` (the log-log
slope diagnostic, with an optional power-law fit), and :func:`planform` (a
schematic network map). :class:`fluvtree.model.FluvTree` has ``plot`` /
``plot_slope_area`` / ``plot_planform`` conveniences over them.
"""

import numpy as np
import matplotlib.pyplot as plt

from fluvtree.analysis import network as _analysis


def long_profile(network, ax=None, color="C0", connect_baselevel=True, **kwargs):
    """
    Plot the network long profile: elevation vs downstream distance.

    Each reach is drawn from its own ``x``/``z`` arrays and the reaches are joined
    across confluences (and, by default, down to base level at the outlet) with
    connector segments -- so a branching network reads as one continuous profile,
    as in GRLP's examples.

    Parameters
    ----------
    network : RiverNetwork
        The canonical network carrying per-segment ``x``, ``z`` and (for the
        base-level connector) the outlet node's ``x_bl``/``z_bl``.
    ax : matplotlib Axes, optional
        Axes to draw on; a new figure/axes is made if omitted.
    color : optional
        A single colour for the whole network, so it reads as one profile.
    connect_baselevel : bool, optional
        Draw a connector from each mouth reach down to base level (default True).
    **kwargs
        Passed through to ``ax.plot`` (e.g. ``lw``, ``ls``, ``alpha``).

    Returns
    -------
    matplotlib Axes
        The axes drawn on (not shown -- call ``plt.show()`` or use in a notebook).
    """
    if ax is None:
        _, ax = plt.subplots()
    for s in network.segment_ids:
        x = network.get_segment_field(s, "x")
        z = network.get_segment_field(s, "z")
        ax.plot(x, z, color=color, **kwargs)
        d = network.downstream_segment(s)
        if d is not None:                       # join to the downstream reach
            xd = network.get_segment_field(d, "x")
            zd = network.get_segment_field(d, "z")
            ax.plot([x[-1], xd[0]], [z[-1], zd[0]], color=color, **kwargs)
        elif connect_baselevel:                 # mouth reach -> base level
            outlet = network.edge_of(s)[1]
            try:
                x_bl = network.get_node_field(outlet, "x_bl")
                z_bl = network.get_node_field(outlet, "z_bl")
            except KeyError:
                pass
            else:
                ax.plot([x[-1], x_bl], [z[-1], z_bl], color=color, **kwargs)
    ax.set_xlabel("Downstream distance [m]")
    ax.set_ylabel("Elevation [m]")
    return ax


_ABSCISSA_LABEL = {"Q": "Discharge [m$^3$/s]", "A": "Drainage area [m$^2$]"}


def slope_area(network, against="Q", ax=None, fit=True, **kwargs):
    """
    The log-log slope diagnostic: bed slope vs an area-like abscissa.

    Slope ``|dz/dx|`` is computed at each interior interval of every reach and
    plotted against ``against`` at the interval midpoint. FluvTree networks carry
    discharge, so ``against`` defaults to ``"Q"`` (slope-discharge -- the
    transport-relevant form); pass ``against="A"`` for the classic slope-area if a
    drainage-area field is present. With ``fit`` (default), a straight line is fit
    in log-log space and its exponent annotated (for ``"A"`` its negative is the
    concavity index).

    Parameters
    ----------
    network : RiverNetwork
    against : str, optional
        Per-segment field for the abscissa (default ``"Q"``; ``"A"`` for area).
    ax : matplotlib Axes, optional
    fit : bool, optional
        Fit and draw a log-log power law (default True).
    **kwargs
        Passed to the scatter ``ax.loglog`` call.

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    A, S = _analysis.slope_area(network, against)   # analysis computes, plot draws
    kwargs.setdefault("marker", ".")
    kwargs.setdefault("linestyle", "none")
    ax.loglog(A, S, **kwargs)
    if fit and S.size >= 2:
        slope, intercept = np.polyfit(np.log10(A), np.log10(S), 1)
        af = np.array([A.min(), A.max()])
        ax.loglog(af, 10.0 ** intercept * af ** slope, "-", color="k", lw=1)
        ax.annotate("exponent = %.3f" % slope, xy=(0.05, 0.05),
                    xycoords="axes fraction")
    ax.set_xlabel(_ABSCISSA_LABEL.get(against, against))
    ax.set_ylabel("Slope")
    return ax


def _schematic_lanes(network, spacing=1.0):
    """Assign each reach a schematic y-lane: channel heads get evenly-spaced lanes
    in depth-first order, and every confluence sits at the mean of its tributaries'
    lanes (a dendrogram layout -- no overlaps on a convergent tree)."""
    lane = {}
    counter = [0]

    def assign(seg):
        ups = sorted(network.upstream_segments(seg))
        if not ups:
            lane[seg] = counter[0] * spacing
            counter[0] += 1
        else:
            for u in ups:
                assign(u)
            lane[seg] = float(np.mean([lane[u] for u in ups]))

    for mouth in sorted(network.mouth_segments()):
        assign(mouth)
    return lane


def planform(network, ax=None, spacing=1.0, **kwargs):
    """
    A schematic map of the network: real downstream distance across, branches apart.

    Each reach is drawn at its own y-lane over its true ``x`` range and joined to
    its downstream reach at the confluence, so the network's branching structure is
    legible (the FluvTree analogue of GRLP's ``Network.plot`` schematic). The
    ``y`` axis is schematic -- lanes, not coordinates -- since FluvTree's ``x`` is
    downstream distance, not a map coordinate; a true planform awaits map-coordinate
    fields (e.g. from DEM extraction).

    Returns
    -------
    matplotlib Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    lane = _schematic_lanes(network, spacing=spacing)
    for s in network.segment_ids:
        x = np.asarray(network.get_segment_field(s, "x"), float)
        line, = ax.plot(x, np.full(x.shape, lane[s]), **kwargs)
        d = network.downstream_segment(s)
        if d is not None:                    # join to the downstream reach's lane
            xd = network.get_segment_field(d, "x")
            ax.plot([x[-1], xd[0]], [lane[s], lane[d]],
                    color=line.get_color(), **kwargs)
    ax.set_xlabel("Downstream distance [m]")
    ax.set_yticks([])
    ax.set_ylabel("branches (schematic)")
    return ax
