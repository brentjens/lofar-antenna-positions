"""Functions for geographic transformations commonly used for LOFAR"""

from numpy import sqrt, sin, cos, arctan2, array, cross, dot, float64, vstack, transpose, shape
from numpy.linalg import norm


def normalized_earth_radius(latitude_rad):
    """Compute the normalized radius of the WGS84 ellipsoid at a given latitude"""
    wgs84_f = 1. / 298.257223563
    return 1.0 / sqrt(cos(latitude_rad) ** 2 + ((1.0 - wgs84_f) ** 2) * (sin(latitude_rad) ** 2))


def geographic_from_xyz(xyz_m):
    """Compute longitude, latitude and height (from the WGS84 ellipsoid) of a given point or list
       of points

    Args:
        xyz_m (Union[array, list]): xyz-coordinates (in m) of the given point.

    Returns:
        Dict[Union[np.array, float]: Dictionary with 'lon_rad', 'lat_rad', 'height_m',
                                     values are float for a single input, arrays for multiple inputs

    Example:
        >>> from pprint import pprint
        >>> xyz_m = [3836811, 430299, 5059823]
        >>> pprint(geographic_from_xyz(xyz_m))
        {'height_m': -0.28265954554080963,
         'lat_rad': 0.9222359279580563,
         'lon_rad': 0.11168348969295486}
        >>> xyz2_m = array([3828615, 438754, 5065265])
        >>> pprint(geographic_from_xyz([xyz_m, xyz2_m]))
        {'height_m': array([-0.28265955, -0.74483879]),
         'lat_rad': array([0.92223593, 0.92365033]),
         'lon_rad': array([0.11168349, 0.11410087])}
    """
    lon_rad, lat_rad, height_m = geographic_array_from_xyz(xyz_m).T
    # For backward compatibility, return floats (rather than shape 1 arrays) for single input
    if shape(xyz_m) == (3,):
        lon_rad = lon_rad[0]
        lat_rad = lat_rad[0]
        height_m = height_m[0]
    return {'lon_rad': lon_rad, 'lat_rad': lat_rad, 'height_m': height_m}


def geographic_array_from_xyz(xyz_m):
    r'''
    xyz_m is a (N,3) array.
    Compute lon, lat, and height
    Output an (N, 3) array with latitude (rad), longitude (rad) and height (m)
    '''
    wgs84_a = 6378137.0
    wgs84_f = 1./298.257223563
    wgs84_e2 = wgs84_f*(2.0 - wgs84_f)
    
    x_m, y_m, z_m = transpose(xyz_m)
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


def localnorth_to_etrs(centerxyz_m):
    """
    Compute a matrix that transforms from a local coordinate system tangent
    to the WGS84 ellipsoid to ETRS89 ECEF XYZ coordinates.

    Args:
        centerxyz_m (array): xyz-coordinates of the center of the local coordinate system

    Returns:
        array: 3x3 rotation matrix

    Example:
        >>> import numpy
        >>> print(numpy.array_str(
        ...     localnorth_to_etrs(numpy.array([3801633.868, -529022.268, 5076996.892]))
        ...     , precision = 8, suppress_small = True))
        [[ 0.13782846 -0.79200355  0.59475516]
         [ 0.99045611  0.11021248 -0.08276408]
         [ 0.          0.60048613  0.79963517]]
    """
    center_lonlat = geographic_from_xyz(centerxyz_m)
    ellipsoid_normal = normal_vector_ellipsoid(center_lonlat['lon_rad'], center_lonlat['lat_rad'])
    return projection_matrix(centerxyz_m, ellipsoid_normal)


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
        >>> xyz_from_geographic(-0.1382, 0.9266, 99.115)
        array([3802111.62491437, -528822.82583168, 5076662.15079859])
        >>> coords = array([[-0.1382, 0.9266,  99.115],\
                            [ 0.2979, 0.9123, 114.708]])
        >>> xyz_from_geographic(coords[:,0], coords[:,1], coords[:,2]).T
        array([[3802111.62491437, -528822.82583168, 5076662.15079859],
               [3738960.12012956, 1147998.32536741, 5021398.44437063]])
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
