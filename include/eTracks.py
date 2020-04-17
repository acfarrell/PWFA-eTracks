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
import matplotlib.ticker as ticker

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
  r_sim, xi_sim, t0 = axes()
  Er_sim = transE()
  Ez_sim = longE()
  Bphi_sim = phiB()
try: bounds
except NameError: bounds = None
if bounds is None:
  bounds = getBounds()

def EField(r,xi,axis):
  # axis = 1 refers to xi-axis (longitudinal) field
  # axis = 2 refers to r-axis (transverse) field
  if axis == 2:
    xiDex = find_nearest_index(xi_sim, xi)
    rDex = find_nearest_index(r_sim, r)
    return -Er_sim[rDex,xiDex]
  elif axis == 1:
    xiDex = find_nearest_index(xi_sim, xi)
    rDex = find_nearest_index(r_sim, r)
    return -Ez_sim[rDex, xiDex]

def BForce(r,xi,v1,v2,axis):
  xiDex = find_nearest_index(xi_sim, xi)
  rDex = find_nearest_index(r_sim, r)
  BField =  Bphi_sim[rDex, xiDex]
  if axis == 1:
    return -1.0 * v2 * BField
  else:
    return 1.0 * v1 * BField

def Momentum(r, xi, dt, pr, pz):
  p = math.sqrt(pr**2 + pz**2)
  vr = Velocity(pr,p)
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

def outOfBounds(r,xi):
  xiDex = find_nearest_index(xi_sim, xi)
  rDex = find_nearest_index(r_sim, r)
  
  if bounds[rDex, xiDex] == 1:
    #    print(' electron is out of bounds')
    return True
  return False

def GetTrajectory(r_0,pr_0,vr_0,xi_0,pz_0,vz_0,plot,num):
  #returns array of r v. t
  r_dat, z_dat, t_dat, xi_dat = np.array([]),np.array([]),np.array([]),np.array([])
  p = math.sqrt(pr_0**2 + pz_0**2)
  rn = r_0 # position in c/w_p
  pr = pr_0 # momentum in m_e c
  vrn = pr_0/Gamma(p) # velocity in c
  t = t0 # start time in 1/w_p
  dt = .005 # time step in 1/w_p
  
  z0 = xi_0 + t0
  zn = xi_0 + t0
  pz = pz_0 
  vzn = pz/Gamma(p) 
  
  dvz = 0.0
  dvr = 0.0

  old_r = r_0 #- 1.0
  turnRad = r_0
  xin = xi_0
  
  esc = 1

  #Iterate through position and time using a linear approximation 
  #until the radial position begins decreasing
  i = 0 #iteration counter
  while rn > 0:
  
    #Determine Momentum and velocity at this time and position
    pz, pr, p = Momentum(rn, xin, dt, pr, pz)
    vzn = Velocity(pz,p)  
    vrn = Velocity(pr,p)

    #Add former data points to the data lists
    r_dat = np.append(r_dat, rn)
    t_dat = np.append(t_dat, t)
    z_dat = np.append(z_dat, zn)
    xi_dat = np.append(xi_dat, xin)
    #print("z = ", zn)
    if rn > turnRad:
      turnRad = rn

    #Add the distance traveled in dt to r, increase t by dt
    zn += vzn * dt
    rn += vrn * dt
    t += dt
    xin = zn - t
    i += 1
    
    if i > 10000 or rn > 6 or xin < 0 or xin > 9:
      break
    if outOfBounds(rn, xin): 
      esc, xiPos = -1, xin
  if plot:
    plotTest(r_dat, xi_dat, esc, num)
  if esc != -1:
    xiPos = xin
  del r_dat, xi_dat, z_dat, t_dat
  return esc, xiPos        

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

def plotTest(r,xi, esc, num):

  fig, ax = plt.subplots()

  #Make color axis of electric field
  colors = ax.pcolormesh(xi_sim,r_sim,Er_sim,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-Er_sim.max(),vmax=Er_sim.max()),cmap="RdBu_r")
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  cbar = fig.colorbar(colors,ax=ax,ticks=tick_locations, format=ticker.LogFormatterMathtext())
  cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')    
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')
    
  ax.plot(xi,r,'k',label = "Simulated Trajectory")

  plt.xlim(xi_sim[0], xi_sim[-1])
  ax.legend()
  if esc == 1:
    status = "Captured"
  else:
    status = "Escaped"
  ax.set_title('Electron Trajectory Marked '+status)
  fn = "plots/"+ str(num) + "_" + status + ".png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.close()
