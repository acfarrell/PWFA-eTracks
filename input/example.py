# Here is where the initial conditions of the simulation are defined
# This filename is the input parameter of the beamSim.py file

#Output File Name
output_fname = "FlatTop" 

#Field Map File Names
Er_fname = "Efield_r.h5"
Ez_fname = "Efield_z.h5"
Bphi_fname = "Bfield_phi.h5"
t0 = 858.95 #time in simulation where field maps were recorded

#Beam Parameters
sig_r = .47 #transverse beam spread in c/wp
sig_z = .47 #longitudinal beam spread in c/wp

#Impurity Profile
#Array the length of the field map in xi with values in units of np
import numpy as np
from include.getOsirisFields import axes
r,xi = axes(Er_fname)
n_imp =  np.zeros(len(xi))

#Total Charge to inject in C
injCharge = 1e-9 #1 nC of charge
