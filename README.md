# tpa_analysis
Dihedral angle analysis codes for triphenylamine molecules - all in python
Author: James Dix - 2015

Pre-requisites: numpy, MDAnalysis, nosetests

The analysis code in this repo. is to analyse the internal diheral angle of phenyl rings in triphenylamine molecules
confined between two graphene sheets. It is also to analyse the angle between the plane of the phenyl rings and the 
plane made by the graphene sheet.

The simulations are all run in GROMACS 4.5.4, which is a molecular dynamics software package. To save reading in the
results of these simulations myself, MDAnalysis (see here for more details: https://github.com/MDAnalysis/mdanalysis)
has been used to read the results into a python useable format.
