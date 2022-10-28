"""
Functions for plotting LOFAR stations with matplotlib
"""
from .db import LofarAntennaDatabase
from .geo import (
    localnorth_to_etrs,
    geographic_array_from_xyz,
    xyz_from_geographic,
)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon


db = LofarAntennaDatabase()

__all__ = [
    "plot_hba",
    "plot_lba",
    "plot_station",
    "plot_superterp",
    "plot_core",
    "add_background",
]


def plot_hba(
    station_name, ax=None, centre=None, subfield="", labels=False, tilestyle="lines"
):
    """
    Plot LOFAR HBA tiles for one station

    Args:
        station_name: Station name, without suffix. E.g. "CS001"
        ax: existing matplotlib axes object to use
        centre: etrs coordinates of origin. Default: HBA phase centre of station.
        subfield: '0', '1' or ''to suffix after 'HBA' for e.g. CS002HBA1)
        labels: add labels

    Example:
        >>> from lofarantpos.plotutil import plot_hba
        >>> plot_hba("CS001")
    """
    if centre is None:
        centre = db.phase_centres[station_name + "HBA"]

    if ax is None:
        fig, ax = plt.subplots()
        ax.set_xlabel("Local East (m)")
        ax.set_ylabel("Local North (m)")
        ax.set_aspect("equal")

    etrs_to_xyz = localnorth_to_etrs(centre).T

    if station_name + "HBA" + subfield not in db.hba_rotations:
        plot_hba(
            station_name,
            ax=ax,
            centre=centre,
            subfield="0",
            labels=labels,
            tilestyle=tilestyle,
        )
        plot_hba(
            station_name,
            ax=ax,
            centre=centre,
            subfield="1",
            labels=labels,
            tilestyle=tilestyle,
        )
        return

    etrs_delta = db.antenna_etrs(station_name + "HBA" + subfield) - centre
    xys = (etrs_to_xyz @ etrs_delta.T)[:2, :].T

    theta = db.hba_rotations[station_name + "HBA" + subfield]
    c, s = np.cos(theta), np.sin(theta)
    rot_mat = db.pqr_to_localnorth(station_name + "HBA")[:2, :2] @ np.array(
        ((c, s), (-s, c))
    )
    x_dir = rot_mat @ [5.15, 0]
    y_dir = rot_mat @ [0, 5.15]

    tile = np.array(
        [
            -0.5 * x_dir - 0.5 * y_dir,
            +0.5 * x_dir - 0.5 * y_dir,
            +0.5 * x_dir + 0.5 * y_dir,
            -0.5 * x_dir + 0.5 * y_dir,
            -0.5 * x_dir - 0.5 * y_dir,
        ]
    )

    for num, xy in enumerate(xys):
        x, y = xy
        if tilestyle == "lines":
            ax.plot((x + tile)[:, 0], (y + tile)[:, 1], "k")
        else:
            # Plot transparent lines, to keep extent
            ax.plot((x + tile)[:, 0], (y + tile)[:, 1], "k", alpha=0)
        if tilestyle == "filled":
            ax.add_patch(Polygon((tile + [x, y]), fill=True, facecolor="k"))
        if labels:
            if subfield == "1":
                num += 24
            if tilestyle == "filled":
                ax.text(x, y, str(num), va="center", ha="center", color="w")
            else:
                ax.text(x, y, str(num), va="center", ha="center")


def plot_lba(station_name, ax=None, centre=None, labels=False):
    """
    Plot LOFAR LBA locations for one station

    Args:
        station_name: Station name, without suffix. E.g. "CS001"
        ax: existing matplotlib axes object to use
        centre: etrs coordinates of origin. Default: LBA phase centre of station.
        labels: add labels

    Example:
        >>> from lofarantpos.plotutil import plot_lba
        >>> plot_lba("IE613", labels=True)
    """
    if centre is None:
        centre = db.phase_centres[station_name + "LBA"]

    if ax is None:
        fig, ax = plt.subplots()
        ax.set_xlabel("Local East (m)")
        ax.set_ylabel("Local North (m)")
        ax.set_aspect("equal")

    etrs_to_xyz = localnorth_to_etrs(centre).T
    etrs_delta = db.antenna_etrs(station_name + "LBA") - centre
    xys = (etrs_to_xyz @ etrs_delta.T)[:2, :].T

    ax.plot(xys[:, 0], xys[:, 1], "k.")

    if labels:
        for num, xy in enumerate(xys):
            x, y = xy
            ax.text(x, y, str(num))


def plot_cabinet(station_name, ax=None, centre=None, labels=False):
    """
    Plot LOFAR cabinet location for one station

    Args:
        station_name: Station name, without suffix. E.g. "CS001"
        ax: existing matplotlib axes object to use
        centre: etrs coordinates of origin. Default: LBA phase centre of station.
        labels: add label

    Example:
        >>> from lofarantpos.plotutil import plot_hba, plot_cabinet
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots()
        >>> plot_hba("CS002", ax=ax)
        >>> plot_cabinet("CS002", ax=ax, labels=True)
    """
    if centre is None:
        centre = db.phase_centres[station_name + "LBA"]

    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")

    etrs_to_xyz = localnorth_to_etrs(centre).T
    etrs_delta = db.cabinet_etrs[station_name] - centre
    x, y, _ = etrs_to_xyz @ etrs_delta

    ax.plot(x, y, "ks")
    if labels:
        ax.text(x, y, station_name + " cabinet")


def plot_station(
    station_name, ax=None, centre=None, labels=False, tilestyle="lines", background=None
):
    """
    Plot a LOFAR station

    Args:
        station_name: Station name, without suffix. E.g. "CS001"
        ax: existing matplotlib axes object to use
        labels: add labels
        centre: centre of projection, as ETRS xyz coordinate. Default is station's LBA phase centre.
        tilestyle: style for HBA tiles ("filled", "lines" or None)
        background: name of background to draw, e.g. "openstreetmap" or "luchtfoto" (Dutch stations only)

    Example:
        >>> from lofarantpos.plotutil import plot_station
        >>> plot_station("CS002")
        >>> plot_station("CS011", background="openstreetmap")
        >>> plot_station("RS210", background="luchtfoto")
        >>> plot_station("IE613", background="Stamen_Toner")
    """
    if centre is None:
        centre = db.phase_centres[station_name + "LBA"]

    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.set_xlabel("Local East (m)")
        ax.set_ylabel("Local North (m)")
        ax.set_title("LOFAR station " + station_name)

    plot_lba(station_name, ax=ax, centre=centre, labels=labels)
    plot_hba(station_name, ax=ax, centre=centre, labels=labels, tilestyle=tilestyle)
    plot_cabinet(station_name, ax=ax, centre=centre, labels=labels)

    if background is not None:
        add_background(ax, centre, background)


def plot_superterp(
    ax=None, labels=False, circle=True, centre=None, tilestyle="lines", background=None
):
    """
    Plot the LOFAR superterp

    Args:
        ax: existing matplotlib axes object to use
        labels: add labels
        circle: plot a surrounding circle
        centre: ETRS xyz coordinates of centre (default is the centre of the superterp)
        tilestyle: style for HBA tiles ("lines", "filled" or None)
        background: name of background to draw, e.g. "openstreetmap" or "luchtfoto" (Dutch stations only)

    Example:
        >>> from lofarantpos.plotutil import plot_superterp
        >>> plot_superterp(background='openstreetmap', tilestyle=None)
    """
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.set_xlabel("Local East (m)")
        ax.set_ylabel("Local North (m)")
        ax.set_title("LOFAR superterp")

    if circle:
        ax.add_patch(Circle((0, 0), radius=185, fill=False, edgecolor="k"))

    if centre is None:
        centre = db.phase_centres["CS002LBA"]
    for station_name in "CS002", "CS003", "CS004", "CS005", "CS006", "CS007":
        plot_station(
            station_name,
            ax=ax,
            centre=centre,
            labels=labels,
            tilestyle=tilestyle,
            background=None,
        )

    if background is not None:
        add_background(ax, centre, background)


def plot_core(
    ax=None, labels=False, circle=True, centre=None, tilestyle="lines", background=None
):
    """
    Plot the LOFAR core

    Args:
        ax: existing matplotlib axes object to use
        labels: add labels for each dipole
        circle: plot a circle around the superterp
        centre: ETRS xyz coordinates of centre (default is centre of the superterp)
        tilestyle: style for HBA tiles ("lines", "filled", or None)
        background: name of background to draw, e.g. "openstreetmap" or "luchtfoto"

    Example:
        >>> from lofarantpos.plotutil import plot_core
        >>> plot_core()
    """
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.set_xlabel("Local East (m)")
        ax.set_ylabel("Local North (m)")
        ax.set_title("LOFAR core")

    if centre is None:
        centre = db.phase_centres["CS002LBA"]

    plot_superterp(
        ax=ax,
        labels=labels,
        circle=circle,
        centre=centre,
        tilestyle=tilestyle,
        background=None,
    )

    core_stations = set(
        [
            station_name[:5]
            for station_name in db.phase_centres.keys()
            if station_name[:2] == "CS"
        ]
    )
    superterp_stations = set(["CS002", "CS003", "CS004", "CS005", "CS006", "CS007"])
    lofar_centre = db.phase_centres["CS002LBA"]

    for station in core_stations - superterp_stations:
        plot_station(
            station,
            ax=ax,
            centre=centre,
            labels=labels,
            tilestyle=tilestyle,
            background=None,
        )

    add_background(ax, centre, background, zoom=15)


def _full_extent_to_xy(plotter, centre):
    """
    Convert plotter extents to local north coordinates w.r.t. centre

    Args:
        plotter: tilemapbase.plotter object
        centre: coordinates of centre, in ETRS xyz
    """
    try:
        import tilemapbase.mapping
    except ImportError:
        raise RuntimeError("To plot background maps, tilemapbase is required")

    scale = 2 ** plotter.zoom
    xmin_deg, ymin_deg = tilemapbase.mapping.to_lonlat(
        *plotter.extent.project(plotter.xtilemin / scale, plotter.ytilemin / scale)
    )
    xmax_deg, ymax_deg = tilemapbase.mapping.to_lonlat(
        *plotter.extent.project(
            (plotter.xtilemax + 1) / scale, (plotter.ytilemax + 1) / scale
        )
    )

    lower_left_etrs = xyz_from_geographic(np.deg2rad(xmin_deg), np.deg2rad(ymin_deg), 0)
    upper_right_etrs = xyz_from_geographic(
        np.deg2rad(xmax_deg), np.deg2rad(ymax_deg), 0
    )

    etrs_to_xyz = localnorth_to_etrs(centre).T
    lower_left_xyz = etrs_to_xyz @ (lower_left_etrs - centre)
    upper_right_xyz = etrs_to_xyz @ (upper_right_etrs - centre)

    return (
        lower_left_xyz[0],
        upper_right_xyz[0],
        upper_right_xyz[1],
        lower_left_xyz[1],
    )


def add_background(ax, centre, background, zoom=18):
    """
    Add openstreetmap background to axes. Assumes the axes extents are given in metres w.r.t. a local coordinates system
    with origin centre (given in ETRS XYZ)

    Args:
        ax: existing matplotlib Axes object
        centre: ETRS XYZ coordinates of the centre
        background: name of background to draw, e.g. "openstreetmap" or "luchtfoto" (Dutch stations only)
        zoom: zoom level for the background tiles

    Example:
        >>> from lofarantpos.plotutil import plot_station, add_background
        >>> from lofarantpos.db import LofarAntennaDatabase
        >>> db = LofarAntennaDatabase()
        >>> centre = db.phase_centres["CS002LBA"]
        >>> fig, ax = plt.subplots()
        >>> plot_station("CS002", ax=ax)
        >>> plot_station("CS001", ax=ax, centre=centre)
        >>> add_background(ax, centre, 'OSM')
    """
    try:
        import tilemapbase
    except ImportError:
        raise RuntimeError("To plot background maps, tilemapbase is required")

    if background == "osm" or background == "openstreetmap":
        t = tilemapbase.tiles.build_OSM()
    elif background in [
        "Carto_Dark",
        "Carto_Dark_Labels",
        "Carto_Dark_No_Labels",
        "Carto_Light",
        "Carto_Light_Labels",
        "Carto_Light_No_Labels",
        "Stamen_Terrain",
        "Stamen_Terrain_Background",
        "Stamen_Terrain_Labels",
        "Stamen_Terrain_Lines",
        "Stamen_Toner",
        "Stamen_Toner_Background",
        "Stamen_Toner_Hybrid",
        "Stamen_Toner_Labels",
        "Stamen_Toner_Lines",
        "Stamen_Toner_Lite",
        "Stamen_Watercolour",
    ]:
        t = getattr(tilemapbase.tiles, background)
    elif background == "luchtfoto":
        t = tilemapbase.tiles.Tiles(
            "https://service.pdok.nl/hwh/luchtfotorgb/wmts/v1_0/2022_ortho25/EPSG:3857/{zoom}/{x}/{y}.jpeg",
            "LUFO2022",
        )
    elif isinstance(background, tilemapbase.tiles.Tiles):
        t = background
    else:
        raise ValueError("Background not recognized: " + str(background))

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    tilemapbase.init(create=True)

    xmin_deg, ymin_deg, _ = np.rad2deg(
        geographic_array_from_xyz(
            centre + (localnorth_to_etrs(centre) @ [xmin, ymin, 0])
        )
    )[0]
    xmax_deg, ymax_deg, _ = np.rad2deg(
        geographic_array_from_xyz(
            centre + (localnorth_to_etrs(centre) @ [xmax, ymax, 0])
        )
    )[0]

    extent = tilemapbase.Extent.from_lonlat(xmin_deg, xmax_deg, ymin_deg, ymax_deg)

    plotter = tilemapbase.Plotter(extent, t, zoom=zoom)

    im = plotter.as_one_image()

    ax.imshow(im, extent=_full_extent_to_xy(plotter, centre), zorder=0)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
