"""LOFAR antenna database

Module for manipulating LOFAR antenna databases. Typical usage is to create
an instance of a LofarAntennaDatabase:

>>> import lofarantpos, numpy
>>> db = lofarantpos.db.LofarAntennaDatabase()
>>> db.phase_centres['CS001LBA']
array([3826923.942,  460915.117, 5064643.229])
>>> numpy.set_printoptions(suppress=True)
>>> db.antenna_pqr('RS210LBA')[:5]
array([[ 0.        ,  0.        ,  0.        ],
       [-0.00006...,  2.55059...,  0.00185...],
       [ 2.24997...,  1.3499502 ,  0.00130...],
       [ 2.24982..., -1.35031..., -0.0004149 ],
       [ 0.00006..., -2.55059..., -0.00185...]])
>>> db.container_locations['CS002']
array([3826609.602,  460990.583, 5064879.514])
"""
import csv
import os
import pathlib

import numpy

from lofarantpos import geo


def install_prefix():
    path_elements = pathlib.PurePath(__file__).parts
    path_to_module = path_elements[:-2]
    if path_to_module[-1] == 'site-packages':
        if 'python' in path_to_module[-2] and 'lib' in path_to_module[-3]:
            return os.path.join(*(path_to_module[:-3]))
    return os.path.join(*path_to_module)


def parse_csv(file_name, data_type):
    """Read a CSV file and convert the elements to a given data type

    Args:
        file_name (str): name of file to be read
        data_type (type): type to convert each line of the CSV to (e.g. `PhaseCentre`)

    Returns:
        list: List of objects of the given type
    """
    return [data_type(row)
            for row_id, row in enumerate(csv.reader(open(file_name)))
            if row_id > 0]


def getcol(rows, col_name):
    return [row.__dict__[col_name]
            for row in rows]


def parse_hba_rotations(file_name):
    hba_rotations = {}
    for row in parse_csv(file_name, list):
        if row[2].strip() == '':
            hba_rotations[row[0] + 'HBA'] = float(row[1]) * numpy.pi / 180.0
        else:
            hba_rotations[row[0] + 'HBA0'] = float(row[1]) * numpy.pi / 180.0
            hba_rotations[row[0] + 'HBA1'] = float(row[2]) * numpy.pi / 180.0
    return hba_rotations


class Antenna(object):
    def __init__(self, csv_row):
        self.station = csv_row[0]
        self.antenna_type = csv_row[1]
        self.antenna_id = int(csv_row[2])
        self.etrs = numpy.array([float(csv_row[3]),
                                 float(csv_row[4]),
                                 float(csv_row[5])])
        self.rcu_x = int(csv_row[6])
        self.rcu_y = int(csv_row[7])

    def __repr__(self):
        return repr(self.__dict__)


class PhaseCentre(object):
    def __init__(self, csv_row):
        self.station = csv_row[0]
        self.field = csv_row[1]
        self.etrs = numpy.array([float(csv_row[2]),
                                 float(csv_row[3]),
                                 float(csv_row[4])])

    def __repr__(self):
        return repr(self.__dict__)


class ContainerLocation(object):
    def __init__(self, csv_row):
        self.station = csv_row[0]
        self.etrs = numpy.array([float(csv_row[14]),
                                 float(csv_row[15]),
                                 float(csv_row[16])])

    def __repr__(self):
        return repr(self.__dict__)


class RotationMatrix(object):
    def __init__(self, csv_row):
        self.station = csv_row[0]
        self.field = csv_row[1]
        self.matrix = numpy.array([float(x) for x in csv_row[2:]]).reshape((3, 3))

    def __repr__(self):
        return repr(self.__dict__)


class LofarAntennaDatabase(object):
    """Database with LOFAR antenna positions

    This database contains the LOFAR antenna positions in both ETRS and the
    station/field specific pqr coordinate system. The upstream source is
    the LOFAR svn repository at https://svn.astron.nl/LOFAR.

    Attributes:
        antennas (list): all antenna information
        phase_centres (dict): ETRS phase centres for each antenna field
        hba_rotations (dict): HBA rotations (in radians) for each antenna field
        pqr_to_etrs (dict): Rotation matrix from PQR to ETRS for each antenna field
    """

    def __init__(self, path_to_files=None):
        if path_to_files is None:
            # Install_prefix can end up to be some_path/lib/site_packages,
            # append to the search path the install_prefix minus last two directories
            search_path = [install_prefix(),
                           os.sep.join(install_prefix().split(os.sep)[:-2]),
                           '/usr/local', '/usr']
            for attempt in search_path:
                share = os.path.join(attempt, os.path.join('share', 'lofarantpos'))
                if os.path.exists(os.path.join(share, 'etrs-phase-centres.csv')):
                    break
        else:
            share = path_to_files
        self.phase_centres = {
            c.station + c.field: c.etrs
            for c in parse_csv(os.path.join(share, 'etrs-phase-centres.csv'),
                               PhaseCentre)}
        self.container_locations = {
            c.station: c.etrs
            for c in parse_csv(os.path.join(share, 'stationinfo.csv'), ContainerLocation)
        }
        self.antennas = parse_csv(os.path.join(share, 'etrs-antenna-positions.csv'),
                                  Antenna)
        pqr_to_etrs_rows = parse_csv(os.path.join(share, 'rotation_matrices.dat'),
                                     RotationMatrix)
        self.pqr_to_etrs = {m.station + m.field: m.matrix for m in pqr_to_etrs_rows}
        self.hba_rotations = parse_hba_rotations(os.path.join(share, 'hba-rotations.csv'))
        core_stations = numpy.unique([name[0:5] for name in self.phase_centres.keys()
                                      if 'CS' in name])
        for core_station in core_stations:
            self.pqr_to_etrs[core_station + 'HBA'] = self.pqr_to_etrs[core_station + 'HBA0']

    #            self.antennas[core_station+'HBA'] = numpy.concatenate([self.antennas[core_station+'HBA0'],
    #                                                                   self.antennas[core_station+'HBA0']],
    #                                                                  axis=0)

    def __repr__(self):
        return repr(self.__dict__)

    def antenna_etrs(self, field_name):
        """Return a list of all ETRS antenna coordinates for a given antenna field

        Args:
            field_name (str): Field name (e.g. 'CS001HBA0')

        Returns:
            array: array of ETRS coordinates
        """
        station = field_name[0:5].upper()
        subfield = field_name[5:].upper()
        antenna_ids = {'LBA': numpy.arange(2048),
                       'HBA': numpy.arange(2048),
                       'HBA0': numpy.arange(0, 24),
                       'HBA1': numpy.arange(24, 48)}
        return numpy.array(
            getcol(
                sorted([ant for ant in self.antennas
                        if ant.station == station
                        and ant.antenna_type == subfield[0:3]
                        and ant.antenna_id in antenna_ids[subfield]],
                       key=lambda x: x.antenna_id),
                'etrs'))

    def antenna_pqr(self, field_name):
        """Return a list of all PQR antenna coordinates for a given antenna field

        Args:
            field_name (str): Field name (e.g. 'CS001HBA0')

        Returns:
            array: array of PQR coordinates
        """
        return geo.transform(
            self.antenna_etrs(field_name),
            self.phase_centres[field_name],
            self.pqr_to_etrs[field_name].T)

    def hba_dipole_pqr(self, field_name):
        """Return a list of all PQR dipole coordinates for a given HBA antenna field

        Args:
            field_name (str): Field name (e.g. "CS001HBA0")

        Returns:
            array: array of PQR coordinates

        Example:
            >>> import lofarantpos.db
            >>> import numpy
            >>> db = lofarantpos.db.LofarAntennaDatabase()
            >>> db.hba_dipole_pqr("CS001HBA0")[:5]
            array([[ 1.9336444 , 15.284...  ,  0.00008769],
                   [ 3.075576  , 14.776116  ,  0.00008769],
                   [ 4.217508  , 14.267695  ,  0.00008769],
                   [ 5.3594... , 13.7592745 ,  0.00008769],
                   [ 1.4252236 , 14.142605  ,  0.00008769]], dtype=float32)
        """
        base_tile = numpy.array([[[-1.5, 1.5], [-0.5, 1.5], [+0.5, 1.5], [+1.5, +1.5]],
                                 [[-1.5, 0.5], [-0.5, 0.5], [+0.5, 0.5], [+1.5, +0.5]],
                                 [[-1.5, -0.5], [-0.5, -0.5], [+0.5, -0.5], [+1.5, -0.5]],
                                 [[-1.5, -1.5], [-0.5, -1.5], [+0.5, -1.5], [+1.5, -1.5]]],
                                dtype=numpy.float32)
        base_tile *= 1.25
        base_tile_delta_pqr = base_tile.reshape((-1, 2))
        rotation = self.hba_rotations[field_name]
        matrix = numpy.array([[numpy.cos(rotation), numpy.sin(rotation)],
                              [-numpy.sin(rotation), numpy.cos(rotation)]],
                             dtype=numpy.float32)
        rotated_tile_pqr = numpy.dot(matrix, base_tile_delta_pqr.T).T
        antenna_pqr = self.antenna_pqr(field_name)
        return numpy.array([[element[0] + ant[0], element[1] + ant[1], ant[2]]
                            for ant in antenna_pqr
                            for element in rotated_tile_pqr],
                           dtype=numpy.float32).reshape((-1, 3))

    def hba_dipole_etrs(self, field_name):
        """Return a list of all ETRS dipole coordinates for a given HBA antenna field

        Args:
            field_name (str): Field name (e.g. 'CS001HBA0')

        Returns:
            array: array of ETRS coordinates

        Example:
            >>> import lofarantpos.db
            >>> import numpy
            >>> db = lofarantpos.db.LofarAntennaDatabase()
            >>> db.hba_dipole_etrs("IE613HBA")[:5]
            array([[3801679.57332033, -528959.80788382, 5076969.80405122],
                   [3801680.56726901, -528959.55814198, 5076969.08837304],
                   [3801681.56121763, -528959.30839824, 5076968.37269509],
                   [3801682.55516625, -528959.0586545 , 5076967.65701715],
                   [3801679.7113895 , -528961.02799576, 5076969.57003303]])
        """
        return geo.transform(
            self.hba_dipole_pqr(field_name),
            numpy.zeros(3),
            self.pqr_to_etrs[field_name]) + \
               self.phase_centres[field_name][numpy.newaxis, :]
