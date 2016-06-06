import csv, numpy

def install_prefix():
    return __file__

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
    for row in parse_csv('data/hba-rotations.csv', list):
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


