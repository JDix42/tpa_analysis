import numpy as np
import MDAnalysis as MDA

def vect_calc(r1, r2):
    """
    This function calculates the vector from r1 to r2, ie. it
    does r2-r1

    r1 and r2 are assumed to be array of the same shape
    """
    r12x = r2[:,0] - r1[:,0]
    r12y = r2[:,1] - r1[:,1]
    r12z = r2[:,2] - r1[:,2]

    r12 = np.zeros((np.shape(r1)))
    r12[:,0] = r12x[:]
    r12[:,1] = r12y[:]
    r12[:,2] = r12z[:]

    return r12

def norm_vect(r):
    """
    This function normalizes the numpy array r using the standard
    identity matrix of 1
    """
    r_squar = np.power(r[:,:], 2)
    r_sum = np.sum(r_squar, axis=1)
    r_sqrt = np.sqrt(r_sum)

    r_sqrt_3coord = np.zeros((np.shape(r)))
    r_sqrt_3coord[:,0] = r_sqrt[:]
    r_sqrt_3coord[:,1] = r_sqrt[:]
    r_sqrt_3coord[:,2] = r_sqrt[:]

    r_norm = np.divide(r[:,:], r_sqrt_3coord[:,:])
    
    return r_norm

def dot_product(r1, r2):
    """ 
    This function calculates the dot product between vector r1 and
    vector r2
    """
    d1 = np.multiply(r1[:,0], r2[:,0])
    d2 = np.multiply(r1[:,1], r2[:,1])
    d3 = np.multiply(r1[:,2], r2[:,2])
    
    d = d1 + d2 + d3
    
    return d

