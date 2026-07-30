"""
Views of a river network -- the plotting layer.

This is the *only* module that imports matplotlib. The numerics core
(``network`` / ``solvers`` / ``closures``) stays dependency-light and never
imports this module; the arrow points one way (``plot`` reads the core), so bare
``import fluvtree`` pulls in no plotting. Plot functions take a
:class:`~fluvtree.network.RiverNetwork` and draw onto a matplotlib ``Axes``;
:meth:`fluvtree.model.FluvTree.plot` is the convenience that calls
:func:`long_profile` on its network.

Currently: :func:`long_profile` -- the canonical GRLP-style long profile
(elevation vs downstream distance, reaches joined across confluences). Planform and
slope-area views are the natural next additions.
"""

import matplotlib.pyplot as plt


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
