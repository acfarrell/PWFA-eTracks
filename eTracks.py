# Electron Tracking Script
# Author: Audrey Claire Farrell - audrey.farrell@stonybrook.edu
#   This script is designed to track the path of electrons injected in
#   electron-driven plasma wakefield accelerators.

import sys
import numpy as np
import plotTracks 
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.colors as col

#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
N = 1e23                          #electron number density in 1/m^3
W = np.sqrt(N*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

def EField(r):
  #returns the electric field at position r (from Wei Lu's paper)
  return -1.0/4.0 * r

def Momentum(r, dt, p_0):
  #Returns the momentum at t + dt, in units of m_e 
  return p_0 +  EField(r) * dt

def Velocity(p):
  #returns the velocity from the momentum, in units of c
  return p

def GetTrajectory(r_0,p_0):
  #returns array of r v. t
  r_dat = []
  t_dat = []

  rn = r_0 # position in c/w_p
  pn = p_0 # momentum in m_e c
  vn = Velocity(pn) # velocity in c
  t = 0.0 # start time in 1/w_p
  dt = .01 # time step in 1/w_p

  old_r = r_0 - 1.0
        
  #Iterate through position and time using a linear approximation 
  #until the radial position begins decreasing
  while rn > old_r:
#               print("r = " + str(rn) + " c/w_p, v = ", vn)

  #Determine Momentum and velocity at this time and position
    pn = Momentum(rn, dt, pn)
    vn = Velocity(pn)
          
    #Add former data points to the data lists
    r_dat.append(rn)
    t_dat.append(t)

    #Add the distance traveled in dt to r, increase t by dt
    old_r = rn
    rn += vn*dt
    t += dt
  print("\n Turn Radius = ",rn)
  return r_dat,t_dat        

def getEfromSim():
  f = h5.File("simulated_data/fields/e2-000066.h5","r")
  datasetNames = [n for n in f.keys()] #Two Datasets: AXIS and e2
  E_dat = f['e2'][:]
  E_dat = E_dat.astype(float)
  a1_bounds = f['AXIS']['AXIS1']#.astype(float)
  a2_bounds = f['AXIS']['AXIS2']#.astype(float)
  
  a1_dat = np.linspace(a1_bounds[0],a1_bounds[1],len(E_dat[0]))
  a2_dat = np.linspace(a2_bounds[0],a2_bounds[1],len(E_dat))

  fig, ax = plt.subplots()
  colors = ax.pcolormesh(a1_dat,a2_dat,E_dat,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=E_dat.min(),vmax=E_dat.max()),cmap="RdBu_r")
  cbar = fig.colorbar(colors,ax=ax)
  
  cbar.set_label('Electric Field ($m_e c\omega_p / e$)')
  ax.set_xlabel('z ($c/\omega_p$)')
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Transverse Electric Field from Simulation')
  plt.show()
def main():
  #Get initial position and momentum from user input:
  r_0 = float(input("Initial radius (c/w_p): "))
  p_0 = float(input("Initial transverse momentum (m_e c): "))
        
  #Determine trajectory, creates n-length lists of data points
  r_dat, t_dat = GetTrajectory(r_0,p_0)
  #Get number of data points
  n = len(r_dat)
  #Create list of E values for each transverse position
  E_dat = [[-1.0*EField(r) for t in t_dat] for r in r_dat] 
  #Plot the trajectory
  plotTracks.plot(r_dat, t_dat, E_dat)
#main()
getEfromSim()
