import os, csv, pathlib
import numpy
from lofarantpos import geo

def install_prefix():
    path_elements = pathlib.PurePath(__file__).parts
    path_to_module = path_elements[:-2]
    if path_to_module[-1] == 'site-packages':
        if 'python' in path_to_module[-2]:
            if 'lib' in path_to_module[-3]:
                return os.path.join(*(path_to_module[:-3]))
    return os.path.join(*path_to_module)


def parse_csv(file_name, data_type):
    '''
    *Example*
    
    >>> parse_csv('data/etrs-phase-centres.csv', PhaseCentre)
    '''
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
            hba_rotations[row[0]+'HBA'] = float(row[1])*numpy.pi/180.0
        else:
            hba_rotations[row[0]+'HBA0'] = float(row[1])*numpy.pi/180.0
            hba_rotations[row[0]+'HBA1'] = float(row[2])*numpy.pi/180.0
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


class RotationMatrix(object):
    def __init__(self, csv_row):
        self.station = csv_row[0]
        self.field = csv_row[1]
        self.matrix = numpy.array([float(x) for x in csv_row[2:]]).reshape((3, 3))

    def __repr__(self):
        return repr(self.__dict__)


class LofarAntennaDatabase(object):
    def __init__(self, path_to_files=None):
        if path_to_files is None:
            share = os.path.join(install_prefix(), 'share/lofarantpos/')
        else:
            share = path_to_files
        self.phase_centres = {
            c.station+c.field: c.etrs
            for c in parse_csv(os.path.join(share, 'etrs-phase-centres.csv'),
                               PhaseCentre)}
        self.antennas = parse_csv(os.path.join(share, 'etrs-antenna-positions.csv'),
                                  Antenna)
        pqr_to_etrs_rows = parse_csv(os.path.join(share, 'rotation_matrices.dat'),
                            RotationMatrix)
        self.pqr_to_etrs = {m.station+m.field: m.matrix for m in pqr_to_etrs_rows}
        self.hba_rotations = parse_hba_rotations(os.path.join(share, 'hba-rotations.csv'))

    def __repr__(self):
        return repr(self.__dict__)


    def antenna_etrs(self, field_name):
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
        return geo.transform(
            self.antenna_etrs(field_name),
            self.phase_centres[field_name],
            self.pqr_to_etrs[field_name].T)
        
    
    def hba_dipole_pqr(self, field_name):
        base_tile= numpy.array([[[-1.5, 1.5], [-0.5, 1.5], [+0.5, 1.5], [+1.5, +1.5]],
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
        return numpy.array([[element[0]+ant[0], element[1]+ant[1], ant[2]]
                            for ant in antenna_pqr
                            for element in rotated_tile_pqr],
                           dtype=numpy.float32).reshape((-1,3))
    
    def hba_dipole_etrs(self, field_name):
        return  geo.transform(
            self.hba_dipole_pqr(field_name),
            numpy.zeros(3),
            self.pqr_to_etrs[field_name]) + \
            self.phase_centres[field_name][numpy.newaxis,:]
