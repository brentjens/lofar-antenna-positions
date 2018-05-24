#!/usr/bin/env python

from distutils.core import setup
from lofarantpos import __version__

setup(name='lofar-antenna-positions',
      version      = __version__,
      description  = 'Access, query, and manipulate LOFAR antenna positions from ',
      author       = 'M.A. Brentjens',
      author_email = 'brentjens@astron.nl',
      packages     = ['lofarantpos'],
      requires     = ['numpy', 'pathlib'],
      scripts      = [],
      data_files   = [('share/lofarantpos/',
                       ['share/lofarantpos/etrs-antenna-positions.csv',
                        'share/lofarantpos/etrs-phase-centres.csv',
                        'share/lofarantpos/hba-rotations.csv',
                        'share/lofarantpos/normal_vectors.dat',
                        'share/lofarantpos/rotation_matrices.dat'])]
     )

