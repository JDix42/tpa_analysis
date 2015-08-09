import numpy as np
import MDAnalysis as MDA
import sys

sys.path.append("/home/mbdxkjd7/code/tpa_analysis/lib")
from numpy.testing import *
from dihedral import dihedral_calc

class testDihedralCalc(TestCase):
    def setUp(self):
        
        self.d = dihedral_calc

        self.r1 = np.array(([[1,0,0]]))
        self.r2 = np.array(([[0,0,0]]))
        self.r3 = np.array(([[0,1,0]]))

    def tearDown(self):
        del self.d
        
    def _test_dihedral_angle(self, r1, r2, r3, r4, angle):
        a = self.d(r1, r2, r3, r4)
        
        a_deg = np.rad2deg(a)

        assert_almost_equal(a_deg[0], angle,   7, err_msg="Your dihedral angles don't quite add up (to what I don't know)")

    def test_zero_dihedral_angle(self):
        r4 = np.array(([[1,1,0]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, 0)

    def test_pl_ninety_dihedral_angle(self):
        r4 = np.array(([[0,1,1]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, 90)

    def test_ng_ninety_dihedral_angle(self):
        r4 = np.array(([[0,1,-1]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, -90)

    def test_oneeighty_dihedral_angle(self):
        r4 = np.array(([[-1,1,0]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, 180)

    def test_pl_fortyfive_dihedral_angle(self):
        r4 = np.array(([[1,1,1]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, 45)

    def test_pl_one_thirtyfive_dihedral_angle(self):
        r4 = np.array(([[-1,1,1]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, 135)

    def test_ng_fortyfive_dihedral_angle(self):
        r4 = np.array(([[1,1,-1]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, -45)

    def test_ng_one_thrityfive_dihedral_angle(self):
        r4 = np.array(([[-1,1,-1]]))
        self._test_dihedral_angle(self.r1, self.r2, self.r3, r4, -135)
