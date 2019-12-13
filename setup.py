#!/usr/bin/env python

from setuptools import setup
from lofarantpos import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

if __name__ == "__main__":
    setup(name='lofarantpos',
          version      = __version__,
          description  = 'Access, query, and manipulate LOFAR antenna positions',
          long_description = long_description,
          long_description_content_type="text/markdown",
          author       = 'M.A. Brentjens',
          author_email = 'brentjens@astron.nl',
          packages     = ['lofarantpos'],
          url          = "https://github.com/lofar-astron/lofar-antenna-positions",
          requires     = ['numpy', 'pathlib'],
          scripts      = [],
          classifiers  = [
              "Programming Language :: Python :: 3",
              "Programming Language :: Python :: 2",
              "License :: OSI Approved :: Apache Software License",
              "Operating System :: OS Independent",
          ],
          data_files   = [('share/lofarantpos',
                           ['share/lofarantpos/etrs-antenna-positions.csv',
                            'share/lofarantpos/etrs-phase-centres.csv',
                            'share/lofarantpos/hba-rotations.csv',
                            'share/lofarantpos/normal_vectors.dat',
                            'share/lofarantpos/rotation_matrices.dat'])]
         )

