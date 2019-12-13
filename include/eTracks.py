# Electron Tracking Script
# Author: Audrey Claire Farrell - audrey.farrell@stonybrook.edu
#   This script is designed to track the path of electrons injected in
#   electron-driven plasma wakefield accelerators.

# Python imports
import sys
import math
import numpy as np
import h5py as h5
import importlib
import matplotlib.pyplot as plt
import matplotlib.colors as col

# include file imports
from .getOsirisFields import axes, longE, transE, phiB 
from .getBounds import getBounds
from .plotTracks import plot

#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
N = 1e23                          #electron number density in 1/m^3
W = np.sqrt(N*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

# Retrieve simulated fields from OSIRIS simulations
try: r_sim
except NameError: r_sim = None
if r_sim is None:
  r_sim, z_sim, t0 = axes()
  Er_sim = transE()
  Ez_sim = longE()
  Bphi_sim = phiB()
try: bounds
except NameError: bounds = None
if bounds is None:
  bounds = getBounds()

def EField(r,z,axis,SHMmodel):
  # SHMmodel = true returns the electric field at position r (from Wei Lu's paper)
  # SHMmodel = false returns the simulated data from OSIRIS
  # axis = 1 refers to xi-axis (longitudinal) field
  # axis = 2 refers to r-axis (transverse) field
  if axis == 2:
    if SHMmodel:
      E = -1.0/4.0 * r
      return E
    else:
      zDex = find_nearest_index(z_sim, z)
      rDex = find_nearest_index(r_sim, r)
      return -Er_sim[rDex,zDex]
  elif axis == 1:
    if SHMmodel:
      return 0.0
    else:
      zDex = find_nearest_index(z_sim, z)
      rDex = find_nearest_index(r_sim, r)
      return -Ez_sim[rDex, zDex]

def BForce(r,z,v1,v2,axis,model):
  if model:# or z - t0 > 5:
    return 0.0

  zDex = find_nearest_index(z_sim, z)
  rDex = find_nearest_index(r_sim, r)
  BField =  Bphi_sim[rDex, zDex]
  if axis == 1:
    return -1.0 * v2 * BField
  else:
    return 1.0 * (v1 + 1) * BField

def Velocity(r, z, dt, v1, v2, axis, model):
  #returns the velocity from the momentum, in units of c
  F = (EField(r, z, axis, model) + BForce(r,z,v1,v2,axis,model))
  if axis == 1:
    v = v1 + 1
  else:
    v = v2
  if abs(v) > 1:
    print("Error: v exceeds light speed")
    return 0
  dv = (F * dt ) / (Gamma(v) + v**2 * Gamma(v)**3 )
  if axis == 1:
    return v - 1 + dv
  return v + dv

def Gamma(v):
  return  1 / math.sqrt(1 - v**2)

def outOfBounds(r,z):
  zDex = find_nearest_index(z_sim, z)
  rDex = find_nearest_index(r_sim, r)
  
  if bounds[rDex, zDex] == 1:
    print(' electron is out of bounds')
    return True
  return False

def GetTrajectory(r_0,pr_0,vr_0,z_0,pz_0,vz_0,SHM):
  #returns array of r v. t

  r_dat, z_dat, t_dat, xi_dat, E_dat = [],[],[],[],[]

  rn = r_0 # position in c/w_p
  pr0 = pr_0 # momentum in m_e c
  vrn = vr_0 # velocity in c
  t = t0 # start time in 1/w_p
  dt = .001 # time step in 1/w_p
  
  z0 = GetInitialZ(z_0,r_0)
  zn = z0
  pz0 = pz_0
  vzn = vz_0 - 1.0 
  
  old_r = r_0 - 1.0
  turnRad = r_0
  xin = zn - t0
  
  #Iterate through position and time using a linear approximation 
  #until the radial position begins decreasing
  i = 0 #iteration counter
  while rn > 0:

  #Determine Momentum and velocity at this time and position
    vrn = Velocity(rn, zn, dt, vzn, vrn,2, SHM)
        
    vzn = Velocity(rn, zn, dt, vzn, vrn, 1, SHM)

    #Add former data points to the data lists
    r_dat.append(rn)
    t_dat.append(t)
    z_dat.append(zn)
    xi_dat.append(xin)
    E_dat.append( EField(rn, zn, 2, SHM) )
    #print("z = ", zn)
    if rn > turnRad:
      turnRad = rn

    #Add the distance traveled in dt to r, increase t by dt
    zn += vzn * dt
    rn += vrn * dt
    t += dt
    xin = zn - t0
    i += 1
#    print("r = ",rn,  ", xi = ",xin, ", vz = ", vzn)

    if outOfBounds(rn,zn) or rn > 6 or xin < 0:
      return -1        
#  print("\n Turn Radius = ",turnRad)
  #return np.array(r_dat),np.array(z_dat),np.array(t_dat), np.array(xi_dat), np.array(E_dat)        
  return 1        

def GetInitialZ(z_0,r_0):
  if z_0 == -1:
    nhalf = int(len(E_sim[0])/2)
    r0dex = find_nearest_index(r_sim, r_0)
    halfE = E_sim[r0dex,nhalf:]
    mindex = np.argwhere(halfE == np.min(halfE))[0] + nhalf
    return z_sim[mindex][0]
  else:
    return z_0

def find_nearest_index(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

