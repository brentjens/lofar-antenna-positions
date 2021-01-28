# LOFAR antenna database

Module for manipulating LOFAR antenna databases. Typical usage is to create
an instance of a LofarAntennaDatabase:

```python
>>> import lofarantpos.db

>>> db = lofarantpos.db.LofarAntennaDatabase()

>>> cs001lba_etrs = db.phase_centres['CS001LBA']
array([3826923.942,  460915.117, 5064643.229])

>>> db.antenna_pqr('RS210LBA')[:5]
array([[ 0.        ,  0.        ,  0.        ],
       [-0.00006...,  2.55059...,  0.00185...],
       [ 2.24997...,  1.3499502 ,  0.00130...],
       [ 2.24982..., -1.35031..., -0.0004149 ],
       [ 0.00006..., -2.55059..., -0.00185...]])
```

Some functions are included to convert xyz coordinates to latitude, longitude,
height w.r.t. WGS84 ellipsoid.

```python
>>> from lofarantpos.geo import geographic_from_xyz

>>> geographic_from_xyz(cs001lba_etrs)
{'lon_rad': 0.11986275972340964,
 'lat_rad': 0.9234780446647385,
 'height_m': 50.162683041766286}
```

## Installation

This module can be pip-installed:
```
pip install lofarantpos
```

Alternatively, install it from this source with
```
python setup.py install
```

**Note** This package used to be called `lofar-antenna-positions`. It may be
necessary to uninstall `lofar-antenna-positions` before installing.
