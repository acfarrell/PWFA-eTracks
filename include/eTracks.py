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
import matplotlib.pylab as pl

# include file imports
from .getQuickPICFields import axes, longE, transE, phiB 
from .getBounds import getBounds
from .plotTracks import plot
import include.plotPhaseSpace as phaseSpace

trajectories = []

#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817e-12                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
N = 1e23                          #electron number density in 1/m^3
W = np.sqrt(N*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

#t0,r_sim,xi_sim,Er_sim,Ez_sim,Bphi_sim,bounds = 0,[],[],[],[],[],[]

# Retrieve fields from OSIRIS simulations
def InitFields(Er_dat,Ez_dat,Bphi_dat,r,xi,t):
  global t0 
  global r_sim, xi_sim
  global Er_sim
  global Ez_sim
  global Bphi_sim
  global trajectories
  trajectories = []
  t0 = t
  r_sim, xi_sim = r,xi
  Er_sim = Er_dat
  Ez_sim = Ez_dat
  Bphi_sim = Bphi_dat
  global bounds
  bounds = getBounds(Er_sim,r_sim,xi_sim,t0)

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

def Momentum_NoB(r, xi, dt, pr, pz):
  p = math.sqrt(pr**2 + pz**2)
  vr = Velocity(pr,p)
  vz = Velocity(pz,p)

  Fz = (EField(r, xi, 1)) #+ BForce(r,xi,vz,vr,1))
  Fr = (EField(r, xi, 2) )#+ BForce(r,xi,vz,vr,2))
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

def GetTrajectory(r_0,xi_0):
  #returns array of r v. t
  r_dat, z_dat, t_dat, xi_dat, vz_dat,pr_dat,Er_dat = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
  p = 0
  rn = r_0 # position in c/w_p
  pr = 0 # momentum in m_e c
  vrn = pr/Gamma(p) # velocity in c
  t = t0 # start time in 1/w_p
  dt = .005 # time step in 1/w_p
  
  z0 = xi_0 + t0
  zn = xi_0 + t0
  pz = 0 
  vzn = pz/Gamma(p) 
  dvz = 0.0
  dvr = 0.0

  old_r = r_0 #- 1.0
  turnRad = r_0
  xin = xi_0
  
  esc = 2

  #Iterate through position and time using a linear approximation 
  #until the radial position begins decreasing
  i = 0 #iteration counter
  # Iterate while electron energy is under 100 MeV
  while Gamma(p) < 100/.511:
  
    #Determine Momentum and velocity at this time and position
    pz, pr, p = Momentum(rn, xin, dt, pr, pz)
    vzn = Velocity(pz,p)  
    vrn = Velocity(pr,p)

    #Add former data points to the data lists
    r_dat = np.append(r_dat, rn)
    t_dat = np.append(t_dat, t)
    z_dat = np.append(z_dat, zn)
    vz_dat = np.append(vz_dat, vzn)
    xi_dat = np.append(xi_dat, xin)
    pr_dat = np.append(pr_dat, pr)
    Er_dat = np.append(Er_dat,-1* EField(rn,xin,2))
    #print("z = ", zn)
    if rn > turnRad:
      turnRad = rn

    #print("vz=",vzn) 
    #Add the distance traveled in dt to r, increase t by dt
    zn += vzn * dt
    rn += vrn * dt
    t += dt
    xin = zn - t
    i += 1
    
    # Allow for crossing the beam axis
    if rn < 0:
      rn = -rn
      pr = -pr
    if rn > 6 or xin < 0 or xin > 10:
      esc = -1
      break
    if outOfBounds(rn, xin): 
      esc = 1
  xiPos = xin
  global trajectories
  data = [r_dat,xi_dat,pr_dat,Er_dat]
  trajectories.append(data)  #print(esc)

  trackName = "r0="+str(round(r_0,3))
  phaseSpace.plot(r_dat,xi_dat,pr_dat, Er_sim,r_sim,xi_sim,trackName)
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
  
def plotVzTest(fname,t0):
      
  fig, axs = plt.subplots(2,sharex=True)

  #Make color axis of electric field
  colors = axs[0].pcolormesh(xi_sim,r_sim,Er_sim,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-Er_sim.max(),vmax=Er_sim.max()),cmap="RdBu_r")
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  #cbar = fig.colorbar(colors,ax=axs[0],ticks=tick_locations, format=ticker.LogFormatterMathtext())
  #cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')    
  axs[1].set_xlabel("$\\xi$ ($c/\omega_p$)")
  axs[0].set_ylabel('r ($c/\omega_p$)')
  axs[1].set_ylabel('$p_z$ ($c$)')
  
  for i in range(len(trajectories)):
    track = trajectories[i]
    r = track[0][:]
    xi = track[1][:]
    vz = track[2][:]
    axs[0].plot(xi,r,'k',alpha=0.5,linewidth=1)
    axs[1].plot(xi,vz,'k',alpha=0.5,linewidth=1)
  axs[0].set_xlim(xi_sim[0], xi_sim[-1])
  axs[0].set_ylim(0, r_sim[-1])
  axs[0].set_title('Ionized Electron Trajectories, t = '+str(t0)+'$\omega_p^{-1}$')
  fn = "plots/"+fname+"_pz.png"
  plt.savefig(fn,dpi=300)
  #plt.show()
def plotNoBTest(fname,t0):
      
  fig, axs = plt.subplots()

  #Make color axis of electric field
  #colors = axs.pcolormesh(xi_sim,r_sim,Er_sim,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-Er_sim.max(),vmax=Er_sim.max()),cmap="RdBu_r")
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  #cbar = fig.colorbar(colors,ax=axs[0],ticks=tick_locations, format=ticker.LogFormatterMathtext())
  #cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')    
  axs.set_xlabel("$\\xi$ ($c/\omega_p$)")
  axs.set_ylabel('$p_r$ ($m_e c$)')
  
  Btrack = trajectories[0]
  noBtrack = trajectories[1]
  Br = Btrack[0][:]
  Bxi = Btrack[1][:]
  Bpr = Btrack[2][:]
  noBr = noBtrack[0][:]
  noBxi = noBtrack[1][:]
  noBpr = noBtrack[2][:]
  axs.plot(Bxi,Bpr,'k',label="With B")
  axs.plot(noBxi,noBpr,'c--',label="Without B")
  axs.legend()
  #axs.set_xlim(xi_sim[0], xi_sim[-1])
  #axs.set_ylim(0, r_sim[-1])
  axs.set_title('Ionized Electron Trajectories, t = '+str(t0)+'$\omega_p^{-1}$')
  fn = "plots/"+fname+".png"
  plt.savefig(fn,dpi=300)
  #plt.show()
def plotTest(fname,t0,percentCaptured):
      
  fig, axs = plt.subplots()
  rainbow = pl.cm.jet(np.linspace(0,1,len(trajectories)))
  
  #Make color axis of electric field
  colors = axs.pcolormesh(xi_sim,r_sim,Er_sim,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-Er_sim.max(),vmax=Er_sim.max()),cmap="RdBu_r")
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  #cbar = fig.colorbar(colors,ax=axs[0],ticks=tick_locations, format=ticker.LogFormatterMathtext())
  #cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')    
  axs.set_xlabel("$\\xi$ ($c/\omega_p$)")
  axs.set_ylabel('r ($c/\omega_p$)')
  
  for i in range(len(trajectories)):
    track = trajectories[i]
    initXi = track[1][0]
    r = track[0][:]
    xi = track[1][:]
    p = track[2][:]
    axs.plot(xi,r,color=rainbow[-i-1])#,linewidth=1)
  axs.set_xlim(initXi - 1.5, initXi)#xi_sim[0], xi_sim[-1])
  axs.set_ylim(0, 1)#r_sim[-1])
  fig.suptitle("Ionized Electron Trajectories")
  axs.set_title('t = '+str(round(t0,2))+'$\omega_p^{-1}$, '+str(round(percentCaptured,2))+'% electrons captured')
  fn = "plots/"+fname+"_rainbow.png"
  plt.savefig(fn,dpi=300)
  #plt.show()
  phaseSpace.plotAll(trajectories)
def plotFieldTest(fname,t0,percentCaptured):
      
  fig, axs = plt.subplots(1,2,sharey=True)
  axs[0].set_ylabel("$E_r$ ($m_ec\omega_p/e$)")
  axs[1].set_ylabel("$E_r$ ($m_ec\omega_p/e$)")
  axs[0].set_xlabel('r ($c/\omega_p$)')
  axs[1].set_xlabel('$\\xi$ ($c/\omega_p$)')
  
  for i in range(len(trajectories)):
    track = trajectories[i]
    r = track[0][:]
    xi = track[1][:]
    p = track[2][:]
    E = track[3][:]
    axs[0].plot(r,E,'k')#,alpha=0.25,linewidth=1)
    axs[1].plot(xi,E,'k')#,alpha=0.25,linewidth=1)
  axs[0].plot(r_sim,Er_sim[:,int(len(xi_sim)/4)], 'r--',label="Transverse Field within the Wake")
  #axs[1].plot(xi_sim,Er_sim[:,int(len(xi_sim)/4)], 'r--',label="Transverse Field within the Wake")
  axs[0].legend()
  #axs[1].legend()
  axs[0].set_xlim(0,2)
  #axs.set_ylim()
  fig.suptitle("Transverse Electric Field Experienced by Electrons")
  #axs.set_title('t = '+str(round(t0,2))+'$\omega_p^{-1}$, '+str(round(percentCaptured,2))+'% electrons captured')
  fn = "plots/"+fname+"_transverseFields.png"
  fig.set_size_inches(12,8)
  plt.savefig(fn,dpi=300)
  #plt.show()
