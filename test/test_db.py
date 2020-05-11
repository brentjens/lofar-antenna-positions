import unittest
from lofarantpos import geo, db
import numpy as np

class TestLofarGeo(unittest.TestCase):
    def setUp(self):
        self.db = db.LofarAntennaDatabase()

    def test_db_not_empty(self):
        self.assertTrue(len(self.db.antennas) > 1000)
        self.assertTrue(len(self.db.phase_centres) > 50)
        self.assertTrue(len(self.db.hba_rotations) > 50)
        self.assertTrue(len(self.db.pqr_to_etrs) > 50)
