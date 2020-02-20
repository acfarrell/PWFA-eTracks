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
import include.plotTracks as plotTracks 
import include.getOsirisFields as osiris

#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
N = 1e23                          #electron number density in 1/m^3
W = np.sqrt(N*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

# Retrieve simulated fields from OSIRIS simulations
r_sim, z_sim, t0 = osiris.axes()
Er_sim = osiris.transE()
Ez_sim = osiris.longE()
Bphi_sim = osiris.phiB()


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

def Velocity(r, z, dt, vi, vj, dvi, dvj, start, model):
  #returns the velocity from the momentum, in units of c in axis direction
  Fi = (EField(r, z, 1, model) + BForce(r,z,vi,vj,1,model))
  Fj = (EField(r, z, 2, model) + BForce(r,z,vi,vj,2,model))
  #Correct for frame moving in +z at speed of light
  vi = vi + 1
  v = math.sqrt(vi**2 + vj**2)
  print("velocity = ", v)
  if abs(v) > 1:
    print("Error: v exceeds light speed")
#    return 0.0,0.0
  vdv = (vi*dvi + vj*dvj)
  dvi = (Fi*dt / Gamma(v) - Gamma(v)**2 * vi *vj* dvj)/(1+Gamma(v)**2 * vi**2)
  dvj = (Fj*dt / Gamma(v) - Gamma(v)**2 * vj *vi* dvi)/(1+Gamma(v)**2 * vj**2)
  #dvi = dt/Gamma(v) * (Fj*Gamma(v)**2 * vi * vj - Fi*(1+Gamma(v)**2*vj**2))*((Gamma(v)**2*vi*vj)**2 - (1+Gamma(v)**2*vi**2)*(1+Gamma(v)**2*vj**2))
  #dvj = dt/Gamma(v) * (Fi*Gamma(v)**2 * vi * vj - Fj*(1+Gamma(v)**2*vi**2))*((Gamma(v)**2*vi*vj)**2 - (1+Gamma(v)**2*vi**2)*(1+Gamma(v)**2*vj**2))
  
  print("dvi = ",dvi,", dvj = ",dvj)
  return dvi, dvj

def Gamma(v):
  return  1 / math.sqrt(1 - v**2)

def GetTrajectory(r_0,pr_0,vr_0,z_0,pz_0,vz_0,SHM):
  #returns array of r v. t

  r_dat, z_dat, t_dat, xi_dat, E_dat = [],[],[],[],[]

  rn = r_0 # position in c/w_p
  pr0 = pr_0 # momentum in m_e c
  if pr0 > 1:
    vrn = vr_0
  else:
    vrn = pr0/Gamma(pr0) # velocity in c
  t = t0 # start time in 1/w_p
  dt = .0005 # time step in 1/w_p
  
  z0 = GetInitialZ(z_0,r_0)
  zn = z0
  pz0 = pz_0
  if pz_0 > 1:
    vzn = vz_0
  else:
    vzn = pz0/Gamma(pz0) - 1.0 
  print("\n Initial z = ",zn)
  
  dvz = 0.0
  dvr = 0.0

  old_r = r_0 - 1.0
  turnRad = r_0
  xin = zn - t0
    
  #Iterate through position and time using a linear approximation 
  #until the radial position begins decreasing
  i = 0 #iteration counter
  while rn > 0:

  #Determine Momentum and velocity at this time and position
    dvz, dvr = Velocity(rn, zn, dt, vzn, vrn, dvz, dvr, i % 2, SHM)
    vzn = vzn + dvz
    vrn = vrn + dvr

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
    #print("r = ",rn,  ", xi = ",xin, ", vz = ", vzn)

    if xin < 0 or rn > 6:
      print("Tracking quit due to xi or r out of range")
      return np.array(r_dat),np.array(z_dat),np.array(t_dat), np.array(xi_dat), np.array(E_dat)        
    if i > 10000000:
      print("Tracking quit due to more than 10K iterations")
      return np.array(r_dat),np.array(z_dat),np.array(t_dat), np.array(xi_dat), np.array(E_dat)        
  print("\n Turn Radius = ",turnRad)
  return np.array(r_dat),np.array(z_dat),np.array(t_dat), np.array(xi_dat), np.array(E_dat)        

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


def main():
  if len(sys.argv) == 2:
    input_fname = str(sys.argv[1])
    print("Using initial conditions from ",input_fname)
    init = importlib.import_module(input_fname)
    r_0 = init.r_0
    pr_0 = init.pr_0
    pz_0 = init.pz_0
    xi_0 = init.xi_0
    Model = init.SHModel
    vz_0 = init.vz_0
    vr_0 = init.vr_0
    track = init.track
    z_0 = xi_0 + t0
  elif len(sys.argv) == 1:
  #Get initial position and momentum from user input:
    r_0 = float(input("Initial radius (c/w_p): "))
    p_0 = float(input("Initial transverse momentum (m_e c): "))        
    z_0 = float(input("Initial z-position (c/w_p) (Enter -1 for position of injection from OSIRIS): "))
    Model = bool(input("Use SHM Model (True/False): "))
    track = 'med'
  else:
    print("Improper number of arguments. Expected 'python3 eTracks.py' or 'python3 eTracks.py <fname>'")
    return

  if Model:
    print("Using SHM Model")
  #Determine trajectory, creates n-length lists of data points
  r_dat, z_dat, t_dat, xi_dat, E_dat = GetTrajectory(r_0,pr_0,vr_0,z_0,pz_0,vz_0,Model)
  plotTracks.plot(r_dat,z_dat, t_dat,xi_dat, Er_sim, r_sim,z_sim,Model,track)

main()
