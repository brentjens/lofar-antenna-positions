import unittest
from lofarantpos import geo
import numpy as np

class TestLofarGeo(unittest.TestCase):
    def test_xyz_from_geographic_and_back(self):
        geographic = [-0.1382, 0.9266, 99.115]
        xyz = geo.xyz_from_geographic(*geographic)
        g2 = geo.geographic_from_xyz(xyz)
        self.assertAlmostEqual(g2['lon_rad'], geographic[0])
        self.assertAlmostEqual(g2['lat_rad'], geographic[1])
        self.assertAlmostEqual(g2['height_m'], geographic[2])

