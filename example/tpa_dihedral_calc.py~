import numpy as np
import MDAnalysis as MDA

def dihedral_calc(r1, r2, r3, r4):
    """
    This function calculates the  angle made by the vectors normal
    to the planes made by r1, r2, r3 and r2, r3, r4. 

    Each of the positions is assumed to be a 3 by n array for 
    n lists of dihedral sets
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

def time_test(u):
    """
    This function tests if an MDA universe, u, has time or not
    """
    is_there_time = 1
    try:
        u.trajectory.totaltime
    except KeyError:
        is_there_time = 0
        
    return is_there_time

def traj_dihedral_calc(ag, dih_list, max_sep, min_sep):
    """
    This function loops over a trajectory to calculate the dihedral of a set
    of atoms supplied in dih_list
    """
    sep_range = max_sep - min_sep
    
    ag_temp = ag
    time_range = ag_temp.trajectory.totaltime - ag_temp.trajectory.time
    
    sep_rate = sep_range/time_range
    
    for ts in ag.trajectory:
        TPA = ag.selectAtoms("resname TPA")
        #Select the 1st, 2nd, 3rd and 4th atoms in the dihedral groups
        r1_TPA = TPA[dih_list[:,0]].positions
        r2_TPA = TPA[dih_list[:,1]].positions
        r3_TPA = TPA[dih_list[:,2]].positions
        r4_TPA = TPA[dih_list[:,3]].positions

        dih_angle = dihedral_calc(r1_TPA, r2_TPA, r3_TPA, r4_TPA)
        dih_angle_deg = np.rad2deg(dih_angle)

        #Calculate the separation at the corresponding ts.time step
        sep = max_sep - sep_rate * ts.time
    
        dih_file.write(str(ts.time) + " " + str(sep) + " ")
    
        for angle in range(np.shape(dih_angle_deg)[0]):
                    dih_file.write(str(dih_angle_deg[angle]) + "  ")
    
        dih_file.write("\n") 

    return dih_angle_deg

def single_frame_dihedral_calc(ag, dih_list):
    """
    This function calculates the dihedral angles for the list in dih_list for
    a single frame of a simulation
    """
    TPA = u.selectAtoms("resname TPA")
    #Select the 1st, 2nd, 3rd and 4th atoms in the dihedral groups
    r1_TPA = TPA[dih_list[:,0]].positions
    r2_TPA = TPA[dih_list[:,1]].positions
    r3_TPA = TPA[dih_list[:,2]].positions
    r4_TPA = TPA[dih_list[:,3]].positions


    dih_angle = dihedral_calc(r1_TPA, r2_TPA, r3_TPA, r4_TPA)
    dih_angle_deg = np.rad2deg(dih_angle) 
        
    for angle in range(np.shape(dih_angle_deg)[0]):
        dih_file.write(str(dih_angle_deg[angle]) + "  ")
    
    dih_file.write("\n")    

    return dih_angle_deg

###
#Here marks the start of the run section of hte code
###
u=MDA.Universe("conf-pull.gro","conf-pull.trr")

#List of dihedral sets based on the atom index number
dih_list = np.array(([[0, 1, 3, 20],
                      [0, 2, 3, 18],
                      [3, 1, 2, 30],
                      [3, 1, 2, 28],
                      [3, 1, 0, 4],
                      [3, 1, 0, 10]]))

max_sep = 1.4
min_sep = 0.7


time = time_test(u)
dih_file=open("dihedral_angle.dat","w")

if time == 1:
    dih_angle_deg = traj_dihedral_calc(u, dih_list, max_sep, min_sep)  
else:   
    dih_angle_deg = single_frame_dihedral_calc(u, dih_list)     

print dih_angle_deg
