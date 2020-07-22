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
r_sim, xi_sim, t0 = osiris.axes()
Er_sim = osiris.transE()
Ez_sim = osiris.longE()
Bphi_sim = osiris.phiB()


def EField(r,xi,axis):
  # axis = 1 refers to xi-axis (longitudinal) field
  # axis = 2 refers to r-axis (transverse) field
  if axis == 2:
    xiDex = find_nearest_index(xi_sim, xi)
    rDex = find_nearest_index(r_sim, r)
    return -1*Er_sim[rDex,xiDex]
  elif axis == 1:
    xiDex = find_nearest_index(xi_sim, xi)
    rDex = find_nearest_index(r_sim, r)
    return -1*Ez_sim[rDex, xiDex]

def BForce(r,xi,vz,vr,axis):
  xiDex = find_nearest_index(xi_sim, xi)
  rDex = find_nearest_index(r_sim, r)
  BField =  Bphi_sim[rDex, xiDex]
  if axis == 1:
    return -1.0 * vr * BField
  else:
    return  1.0 * vz * BField

def Momentum(r, xi, dt, pr, pz):
  #returns the velocity from the momentum, in units of c in axis direction
  p = math.sqrt(pr**2 + pz**2)
  vr = Velocity(pr,p) 
  #Correct for frame moving in +z at speed of light
  vz = Velocity(pz,p) 
  
  Fz = (EField(r, xi, 1) + BForce(r,xi,vz,vr,1))
  Fr = (EField(r, xi, 2) + BForce(r,xi,vz,vr,2))
  #print("Fz = ",Fz,", Fr = ",Fr)
  pz = pz + Fz * dt
  pr = pr + Fr * dt
  p = math.sqrt(pr**2 + pz**2)
  #print("pz = ",pz,", pr = ",pr)
  return pz, pr, p

def Velocity(pi,p):
  v = pi / Gamma(p)
  return v

def Gamma(p):
  return  math.sqrt(1 + p**2)

def GetTrajectory(r_0,pr_0,vr_0,z_0,pz_0,vz_0,SHM):
  #returns array of r v. t

  r_dat, z_dat, t_dat, xi_dat, E_dat = [],[],[],[],[]
  p = math.sqrt(pr_0**2 + pz_0**2)
  rn = r_0 # position in c/w_p
  pr = pr_0 # momentum in m_e c
  vrn = pr_0/Gamma(p) # velocity in c
  t = t0 # start time in 1/w_p
  dt = .005 # time step in 1/w_p
  
  z0 = GetInitialZ(z_0,r_0)
  zn = z0
  pz = pz_0 
  vzn = pz/Gamma(p) 
  print("\n Initial z = ",zn)
  
  dvz = 0.0
  dvr = 0.0

  old_p = p + 1
  old_r = r_0 #- 1.0
  turnRad = r_0
  xin = zn - t0
  old_xi = xin + 1  
  
  #Iterate through position and time using a linear approximation 
  #until the radial position begins decreasing
  i = 0 #iteration counter
  while Gamma(p) < 100/.511 :
    old_p = p
  #Determine Momentum and velocity at this time and position
    pz, pr, p = Momentum(rn, xin, dt, pr, pz)
    vzn = Velocity(pz,p)  
    vrn = Velocity(pr,p)

    #Add former data points to the data lists
    r_dat.append(rn)
    t_dat.append(t)
    z_dat.append(zn)
    xi_dat.append(xin)
    E_dat.append( EField(rn, zn, 2) )
    
    old_xi = xin
    if rn > turnRad:
      turnRad = rn

    #Add the distance traveled in dt to r, increase t by dt
    zn += vzn * dt
    rn += vrn * dt
    t += dt
    xin = zn - t
    i += 1
    #print("r = ",rn,  ", xi = ",xin, ", vz = ", vzn)
    # Allow for crossing the beam axis
    if rn < 0:
      rn = -rn
      pr = -pr
    if xin < 0 or rn > 6:
      print("Tracking quit due to xi or r out of range")
      break
      #return np.array(r_dat),np.array(z_dat),np.array(t_dat), np.array(xi_dat), np.array(E_dat)        
    if i > 10000000:
      print("Tracking quit due to more than 10K iterations")
      break
      #return np.array(r_dat),np.array(z_dat),np.array(t_dat), np.array(xi_dat), np.array(E_dat)        
  print("\n Turn Radius = ",turnRad)
  r_back, xi_back = GetBackwardsTrajectory(rn,pr,vrn,zn,pz,vzn,t, False)
  return np.array(r_dat), np.array(xi_dat), r_back, xi_back        

def GetBackwardsTrajectory(r_0,pr_0,vr_0,z_0,pz_0,vz_0,t,SHM):
  #returns array of r v. t

  r_dat, z_dat, t_dat, xi_dat, E_dat = [],[],[],[],[]
  p = math.sqrt(pr_0**2 + pz_0**2)
  rn = r_0 # position in c/w_p
  pr = pr_0 # momentum in m_e c
  vrn = pr_0/Gamma(p) # velocity in c
 # start time in 1/w_p
  dt = -.005 # time step in 1/w_p
  
  z0 = GetInitialZ(z_0,r_0)
  zn = z0
  pz = pz_0 
  vzn = pz/Gamma(p) 
  print("\n Initial z = ",zn)
  
  dvz = 0.0
  dvr = 0.0

  old_p = p + 1
  old_r = r_0 #- 1.0
  turnRad = r_0
  xin = zn - t
  old_xi = xin + 1  
  
  #Iterate through position and time using a linear approximation 
  #until the radial position begins decreasing
  i = 0 #iteration counter
  while True :
    if p < 0.1 and p > old_p and xin > xi_sim[-1]/2.0:
      break
    old_p = p
  #Determine Momentum and velocity at this time and position
    pz, pr, p = Momentum(rn, xin, dt, pr, pz)
    vzn = Velocity(pz,p)  
    vrn = Velocity(pr,p)

    #Add former data points to the data lists
    r_dat.append(rn)
    t_dat.append(t)
    z_dat.append(zn)
    xi_dat.append(xin)
    E_dat.append( EField(rn, zn, 2) )
    
    old_xi = xin
    if rn > turnRad:
      turnRad = rn

    #Add the distance traveled in dt to r, increase t by dt
    zn += vzn * dt
    rn += vrn * dt
    t += dt
    xin = zn - t
    i += 1
    #print("r = ",rn,  ", xi = ",xin, ", vz = ", vzn)
    # Allow for crossing the beam axis
    if rn < 0:
      rn = -rn
      pr = -pr
    if xin < 0 or rn > 6:
      print("Tracking quit due to xi or r out of range")
      return np.array(r_dat), np.array(xi_dat)        
    if i > 10000000:
      print("Tracking quit due to more than 10K iterations")
      return np.array(r_dat), np.array(xi_dat)        
  print("\n Turn Radius = ",turnRad)
  return np.array(r_dat), np.array(xi_dat)        

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
  r_dat, xi_dat, r_back, xi_back = GetTrajectory(r_0,pr_0,vr_0,z_0,pz_0,vz_0,Model)
  plotTracks.plotBackTrack(r_dat,xi_dat, r_back, xi_back, Er_sim, r_sim,xi_sim,Model,track)

main()
