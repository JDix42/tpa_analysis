import numpy as np
import MDAnalysis as MDA
from numpy.testing import *
from tpa_dihedral_calc import dot_product, norm_vect, dihedral_calc

class TestDotProduct(TestCase):
    def setUp(self):
        
        self.dp = dot_product

    def tearDown(self):
        del self.dp

    def _dim_vect(self,r):
        d = dot_product(r,r)

        assert_equal( d, 4, err_msg="dot_product is in 1 dimension broken")

    def test_x_dim_vect(self):
        r = np.array(([[2,0,0]]))
        
        self._dim_vect(r)

    def test_y_dim_vect(self): 
        r = np.array(([[0,2,0]]))
        
        self._dim_vect(r)

    def test_z_dim_vect(self):
        r = np.array(([[0,0,2]]))
        
        self._dim_vect(r)

    def test_sum_vect(self):
        r = np.array(([[1,2,3]]))

        d = dot_product(r,r)

        assert_equal(d, 14, err_msg="sum in dot_product is broken")

    def test_block_vect(self):
        r = np.random.rand(100, 3)

        d = dot_product(r, r)

        r_sqar = np.power(r[:,:], 2)
        r_sum = np.sum(r_sqar[:,:], axis=1)

        assert_equal(np.shape(d), np.shape(r_sum), err_msg="The shape of the final dot product is not right")

        assert_almost_equal(d[:], r_sum[:], 8, err_msg="The dot_product sum does not add up (funny, huh?)")

class testNormVect(TestCase):
    def setUp(self):
        
        self.nv = norm_vect

    def tearDown(self):
        del self.nv

    def test_normalisation(self):
        r = np.random.rand(100,3)
        
        n = norm_vect(r)
        d = dot_product(n,n)
        
        ones = np.ones((100))
        
        assert_almost_equal(d, ones, 8, err_msg="Your vectors don't normalise too well")

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
