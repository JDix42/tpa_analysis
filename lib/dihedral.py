import numpy as np
import MDAnalysis as MDA
from vector import vect_calc, norm_vect, dot_product

def dihedral_calc(r1, r2, r3, r4):
    """
    This function calculates the  angle made by the vectors normal
    to the planes made by r1, r2, r3 and r2, r3, r4. 

    Each of the positions is assumed to be a 3 by n array for 
    n lists of dihedral sets

    r1
     \ 
      \
       r2---r3
              \
               \
               r4

    The angle for dihedral will go between +180 and -180 degrees
    """
    r12 = vect_calc(r1, r2)
    r23 = vect_calc(r2, r3)
    r34 = vect_calc(r3, r4)

    r123 = np.cross(r12[:], r23[:])
    r234 = np.cross(r23[:], r34[:])
    
    n123 = norm_vect(r123)
    n234 = norm_vect(r234)
    
    axis_2 = np.cross(n123[:], r23[:])
    axis_2_norm = norm_vect(axis_2)    

    x = dot_product(n123[:], n234[:])
    y = dot_product(axis_2_norm[:], n234[:])
    
    dih_angle = np.arctan2(y,x)

    return dih_angle
