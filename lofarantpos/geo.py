"""Functions for geographic transformations commonly used for LOFAR"""

from numpy import sqrt, sin, cos, arctan2, array, cross, dot, float64
from numpy.linalg.linalg import norm


def normalized_earth_radius(latitude_rad):
    """Compute the normalized radius of the WGS84 ellipsoid at a given latitude"""
    wgs84_f = 1. / 298.257223563
    return 1.0 / sqrt(cos(latitude_rad) ** 2 + ((1.0 - wgs84_f) ** 2) * (sin(latitude_rad) ** 2))


def geographic_from_xyz(xyz_m):
    """Compute longitude, latitude and height (from the WGS84 ellipsoid) of a given point

    Args:
        xyz_m (Union[array, list]): xyz-coordinates (in m) of the given point.

    Returns:
        dict: Dictionary with 'lon_rad', 'lat_rad', 'height_m'
    """
    wgs84_a = 6378137.0
    wgs84_f = 1. / 298.257223563
    wgs84_e2 = wgs84_f * (2.0 - wgs84_f)

    x_m, y_m, z_m = xyz_m
    lon_rad = arctan2(y_m, x_m)
    r_m = sqrt(x_m ** 2 + y_m ** 2)
    # Iterate to latitude solution
    phi_previous = 1e4
    phi = arctan2(z_m, r_m)
    while abs(phi - phi_previous) > 1.6e-12:
        phi_previous = phi
        phi = arctan2(z_m + wgs84_e2 * wgs84_a * normalized_earth_radius(phi) * sin(phi),
                      r_m)
    lat_rad = phi
    height_m = r_m * cos(lat_rad) + z_m * sin(lat_rad) - wgs84_a * sqrt(1.0 - wgs84_e2 * sin(lat_rad) ** 2)
    return {'lon_rad': lon_rad, 'lat_rad': lat_rad, 'height_m': height_m}


def geographic_array_from_xyz(xyz_m):
    r'''
    xyz_m is a (N,3) array.
    Compute lon, lat, and height
    Output an (N, 3) array
    '''
    wgs84_a = 6378137.0
    wgs84_f = 1./298.257223563
    wgs84_e2 = wgs84_f*(2.0 - wgs84_f)
    
    x_m, y_m, z_m = xyz_m.T
    lon_rad = arctan2(y_m, x_m)
    r_m = sqrt(x_m**2 + y_m**2)
    # Iterate to latitude solution
    phi_previous = 1e4
    phi = arctan2(z_m, r_m)
    while (abs(phi -phi_previous) > 1.6e-12).any():
        phi_previous = phi
        phi = arctan2(z_m + wgs84_e2*wgs84_a*normalized_earth_radius(phi)*sin(phi),
                      r_m)
    lat_rad = phi
    height_m = r_m*cos(lat_rad) + z_m*sin(lat_rad) - wgs84_a*sqrt(1.0 - wgs84_e2*sin(lat_rad)**2)
    return vstack((lon_rad, lat_rad, height_m)).T



def xyz_from_geographic(lon_rad, lat_rad, height_m):
    """Compute cartesian xyz coordinates from a longitude, latitude and height (from
    the WGS84 ellipsoid)

    Args:
        lon_rad (Union[float, array]): longitude in radians
        lat_rad (Union[float, array]): latitude in radians
        height_m (Union[float, array]): height in meters

    Returns:
        array: xyz coordinates in meters

    Examples:
        >>> import numpy
        >>> if float(".".join(numpy.__version__.split('.')[:2]))>=1.14: numpy.set_printoptions(legacy=True)
        >>> xyz_from_geographic(-0.1382, 0.9266, 99.115)
        array([ 3802111.62491...,  -528822.82583...,  5076662.15079...])
        >>> coords = array([[-0.1382, 0.9266,  99.115],\
                            [ 0.2979, 0.9123, 114.708]])
        >>> xyz_from_geographic(coords[:,0], coords[:,1], coords[:,2]).T
        array([[ 3802111.62491...,  -528822.82583...,  5076662.15079...],
               [ 3738960.12012...,  1147998.32536...,  5021398.44437...]])
    """
    wgs84_a = 6378137.0
    wgs84_f = 1. / 298.257223563
    wgs84_e2 = wgs84_f * (2.0 - wgs84_f)
    c = normalized_earth_radius(lat_rad)
    f = wgs84_f
    a = wgs84_a
    s = c * (1 - f) ** 2
    return array([(a * c + height_m) * cos(lat_rad) * cos(lon_rad),
                  (a * c + height_m) * cos(lat_rad) * sin(lon_rad),
                  (a * s + height_m) * sin(lat_rad)], dtype=float64)


def normal_vector_ellipsoid(lon_rad, lat_rad):
    return array([cos(lat_rad) * cos(lon_rad),
                  cos(lat_rad) * sin(lon_rad),
                  sin(lat_rad)])


def normal_vector_meridian_plane(xyz_m):
    x_m, y_m, _ = xyz_m
    return array([y_m, -x_m, 0.0]) / sqrt(x_m ** 2 + y_m ** 2)


def projection_matrix(xyz0_m, normal_vector):
    r_unit = normal_vector
    meridian_normal = normal_vector_meridian_plane(xyz0_m)
    q_unit = cross(meridian_normal, r_unit)
    q_unit /= norm(q_unit)
    p_unit = cross(q_unit, r_unit)
    p_unit /= norm(p_unit)
    return array([p_unit, q_unit, r_unit]).T


def transform(xyz_m, xyz0_m, mat):
    """Perform a coordinate transformation on an array of points

    Args:
        xyz_m (array): Array of points
        xyz0_m (array): Origin of transformation
        mat (array): Transformation matrix

    Returns:
        array: Array of transformed points
    """
    offsets = xyz_m - xyz0_m
    return array([dot(mat, offset) for offset in offsets])


LOFAR_XYZ0_m = array([3826574.0, 461045.0, 5064894.5])
LOFAR_REF_MERIDIAN_NORMAL = normal_vector_meridian_plane(LOFAR_XYZ0_m)
