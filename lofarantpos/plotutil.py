from .db import LofarAntennaDatabase
from .geo import (
    localnorth_to_etrs,
    geographic_from_xyz,
    geographic_array_from_xyz,
    xyz_from_geographic,
)

db = LofarAntennaDatabase()

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

import tilemapbase
import tilemapbase.mapping


def plot_hba(station_name, ax=None, centre=None, subfield="", labels=False, tiles=True):
    """
    Plot LOFAR HBA tiles for one station

    Args:
        station_name: Station name, without suffix. E.g. "CS001"
        ax: existing matplotlib axes object to use
        centre: etrs coordinates of origin. Default: LBA phase centre of station.
        subfield: '0', '1' or ''to suffix after 'HBA' for e.g. CS002HBA1)
        labels: add labels

    Example:
        >>> from lofarantpos.plotutil import plot_hba
        >>> plot_hba("CS001")
    """
    if centre is None:
        centre = db.phase_centres[station_name + "LBA"]

    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")

    etrs_to_xyz = localnorth_to_etrs(centre).T

    if station_name + "HBA" + subfield not in db.hba_rotations:
        plot_hba(
            station_name, ax=ax, centre=centre, subfield="0", labels=labels, tiles=tiles
        )
        plot_hba(
            station_name, ax=ax, centre=centre, subfield="1", labels=labels, tiles=tiles
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
        if tiles:
            ax.plot((x + tile)[:, 0], (y + tile)[:, 1], "k")
        else:
            # Plot transparent, to keep extent
            ax.plot((x + tile)[:, 0], (y + tile)[:, 1], "k", alpha=0)
        if labels:
            if subfield == "1":
                num += 24
            ax.text(x, y, str(num))


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
        ax.set_aspect("equal")

    etrs_to_xyz = localnorth_to_etrs(centre).T
    etrs_delta = db.antenna_etrs(station_name + "LBA") - centre
    xys = (etrs_to_xyz @ etrs_delta.T)[:2, :].T

    ax.plot(xys[:, 0], xys[:, 1], "ko")

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

    ax.plot(x, y, "ro")
    if labels:
        ax.text(x, y, station_name + " cabinet")


def plot_station(
    station_name, ax=None, centre=None, labels=False, tiles=True, osm_background=False
):
    """
    Plot a LOFAR station

    Args:
        station_name: Station name, without suffix. E.g. "CS001"
        ax: existing matplotlib axes object to use
        labels: add labels

    Example:
        >>> from lofarantpos.plotutil import plot_station
        >>> plot_station("CS002")
        >>> plot_station("CS002", osm_background=True)
    """
    if centre is None:
        centre = db.phase_centres[station_name + "LBA"]

    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.set_title("LOFAR station " + station_name)

    plot_lba(station_name, ax=ax, centre=centre, labels=labels)
    plot_hba(station_name, ax=ax, centre=centre, labels=labels, tiles=tiles)
    plot_cabinet(station_name, ax=ax, centre=centre, labels=labels)

    if osm_background:
        add_osm_background(ax, centre)


def plot_superterp(
    ax=None, labels=False, circle=True, tiles=True, osm_background=False
):
    """
    Plot the LOFAR superterp

    Args:
        ax: existing matplotlib axes object to use
        labels: add labels
        plot_circle: plot a surrounding circle

    Example:
        >>> from lofarantpos.plotutil import plot_superterp
        >>> plot_superterp(osm_background=True, tiles=False)
    """
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        ax.set_xlabel("Local East (m)")
        ax.set_ylabel("Local North (m)")
        ax.set_title("LOFAR superterp")

    if circle:
        circle1 = ax.add_patch(Circle((0, 0), radius=185, fill=False, edgecolor="k"))

    centre = db.phase_centres["CS002LBA"]
    for station_name in "CS002", "CS003", "CS004", "CS005", "CS006", "CS007":
        plot_station(
            station_name,
            ax=ax,
            centre=centre,
            labels=labels,
            tiles=tiles,
            osm_background=False,
        )

    if osm_background:
        add_osm_background(ax, centre)


def _full_extent_to_xy(plotter, centre):
    """
    Convert plotter extents to local north coordinates w.r.t. centre

    Args:
        plotter: tilemapbase.plotter object
        centre: coordinates of centre, in ETRS xyz
    """
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


def add_osm_background(ax, centre, zoom=17):
    """
    Add openstreetmap background to axes. Assumes the axes extents are given in metres w.r.t. a local coordinates system
    with origin centre (given in ETRS XYZ)

    Args:
        ax: existing matplotlib Axes object
        centre: ETRS XYZ coordinates of the centre
        zoom: zoom level for the used tiles

    Example:
        >>> from lofarantpos.plotutil import plot_station, add_osm_background
        >>> from lofarantpos.db import LofarAntennaDatabase
        >>> db = LofarAntennaDatabase()
        >>> centre = db.phase_centres["CS002LBA"]
        >>> fig, ax = plt.subplots()
        >>> plot_station("CS002", ax=ax)
        >>> plot_station("CS001", ax=ax, centre=centre)
        >>> add_osm_background(ax, centre=centre)
    """
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

    t = tilemapbase.tiles.build_OSM()

    extent = tilemapbase.Extent.from_lonlat(xmin_deg, xmax_deg, ymin_deg, ymax_deg)

    plotter = tilemapbase.Plotter(extent, t, zoom=zoom)

    im = plotter.as_one_image()

    ax.imshow(im, extent=_full_extent_to_xy(plotter, centre), zorder=0)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
