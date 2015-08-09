import numpy as np
import MDAnalysis as MDA
import sys

sys.path.append("/home/mbdxkjd7/code/tpa_analysis/lib")
from dihedral import dihedral_calc

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
#Here marks the start of the run section of the code
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
