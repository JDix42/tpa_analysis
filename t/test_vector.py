import numpy as np
import MDAnalysis as MDA
import sys

sys.path.append("/home/mbdxkjd7/code/tpa_analysis/lib/")
from numpy.testing import *
from vector import dot_product, norm_vect

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


