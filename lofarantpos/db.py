import os, csv, pathlib
import numpy

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
