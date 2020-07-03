# Here is where the initial conditions of the simulation are defined
# This filename is the input parameter of the beamSim.py file

#Output File Name
output_fname = "FlatTop10K" 

#Field Map File Names
Er_fname = "EField_r.h5"
Ez_fname = "EField_z.h5"
Bphi_fname = "BField_phi.h5"
t0 = 858.95 #time in simulation where field maps were recorded

#Beam Parameters
sigma_r = .47 #transverse beam spread in c/wp
sigma_z = .47 #longitudinal beam spread in c/wp

#Impurity Profile
#Array the length of the field map in xi with values in units of np
import numpy as np
from include.getOsirisFields import axes
r,xi = axes(Er_fname,t0)
n_imp =  np.zeros(len(xi))
for i in range(len(xi)):
  if xi[i] < 8 and xi[i] > 6.75:
    n_imp[i] = 0.01
#Total Charge to inject in C
injCharge = 1e-9 #1 nC of charge
