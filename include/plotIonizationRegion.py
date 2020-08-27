import sys
import csv
import math
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tempfile import TemporaryFile as tmp
from scipy.special import gamma

#import getOsirisFields as osiris
import include.getQuickPICFields as quickPIC


#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817e-12                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
N = 5e22                          #electron number density in 1/m^3
WP = np.sqrt(N*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

#Er = quickPIC.transE("quickPIC/exslicexz_00000100.h5")
#Ez = quickPIC.longE("quickPIC/ezslicexz_00000100.h5")
#print(xi[1] - xi[0] C/WP)
r,xi, Er = quickPIC.spliceLowRes("quickPIC/exslicexz_00000117.h5")
r,xi, Ez = quickPIC.spliceLowRes("quickPIC/ezslicexz_00000117.h5")


print("Pixel width = ",(xi[1] - xi[0])/WP," seconds")
#print((M_E * C * i

def W(E0):
  if E0 == 0.0:
    return 0
  # Takes an electric field in normalized units and returns the ionization rate in s^-1
  # Define constants for calculating ionization rate for He
  gsE = 24.5 # unperturbed ground state energy in eV
  Z = 1 # charge number after ionization
  n = 0.746 # effective principal quantum number
  
  E = E0 * (M_E * C * WP)/EC * 1e-9  #Convert normalized field to GV/m

  W0 = 1.52e15 * 4**n * gsE / (n * gamma(2*n)) * (20.5*gsE**(3/2))**(2*n-1)
  try:
    W = W0/((E)**(2*n-1)) * math.exp(-6.83*(gsE**(3/2))/E)
  except:
    return 0.0
  #print(W)
  return W 

def ionRatio(i,j):
        #dt = 14.7/2
  integral = 0
  for n in range(len(xi)-1,j,-1):
    En = math.sqrt((Er[i,n])**2 + (Ez[i,n])**2)
    integral = integral +  W(En) * (xi[n] - xi[n-1])/WP
  ratio = 1 - math.exp(-1*integral)#
  #print(ratio)
  if ratio > 1.0:
    return 1.0
  return ratio

def saveIonizationRegion(r0,xi0,Er0,Ez0):
        #print("Saving Ionization Region Data File")
  global r, xi, Er, Ez, maxRatio
  r = r0
  xi = xi0
  Er = Er0
  Ez = Ez0
  max_r = 0
  max_ratio = 0
  max_E = 0
  eRatio  = np.zeros((len(r),len(xi)))
  for i in range(int(len(r)/2.0), int(3*len(r)/4)):
    integral = 0
    for j in range(len(xi)-1,0,-1):
      En = math.sqrt((Er[i,j])**2 + (Ez[i,j])**2)
      integral = integral +  W(En) * (xi[j] - xi[j-1])/WP
      ratio = 1 - math.exp(-1*integral)#
  #print(ratio)
      #if ratio > 1.0:
        #ratio = 1.0
      #print('Row ',i, "/",len(r),", Column ", j,"/",len(xi), end="\r", flush=True)
      if ratio > .01:
        if r[i] > max_r:
          max_r = r[i]
        if ratio > max_ratio and xi[j] > 1:
          max_ratio = ratio
        if En > max_E and xi[j] > 1:
          max_E = En
        eRatio[i,j] = ratio
  np.savez('ionizationRegion.npz', eRatio=eRatio, max_E=max_E,max_r=max_r, max_ratio=max_ratio)
  #print("Ionization Data Saved")
  return

def clipIonizationRatio():
  global eRatio
  for i in range(len(eRatio)):
    for j in range(len(eRatio[0])):
      if eRatio[i,j] >1.0:
        eRatio[i,j] = 1.0
  return

def calculateTimeStep():
  global dt, eRatio
  dat = np.load('ionizationRegion.npz')
  max_ratio = dat['max_ratio']
  eRatio = dat['eRatio']
  max_E = dat['max_E']
  try:
    dt = max_ratio * WP / W(max_E)
  except:
    dt = 0
  #print("Maximum Radius of Ionization Region = ", max_r)
  #timeStep=.075
  print("Maximum Ionization Fraction = ", max_ratio)
  print("Time Step = ",dt)
  return #dt, eRatio

def init(r,xi,Er,Ez,t0,fname):
  saveIonizationRegion(r,xi,Er,Ez)
  calculateTimeStep()
  clipIonizationRatio()
  maxXi, maxR = plotIonizationRegion(fname,t0)
  return maxXi, maxR

#dt, eRatio = calculateTimeStep()
#dt = (xi[1] - xi[0]) 
#print(dt)#0.1
def Wdt(i,j):
  global dt
  En = math.sqrt((Er[i,j])**2 + (Ez[i,j])**2)
  ratio = W(En) * dt /WP  #1 - math.exp(-1*integral)
  #print(W(En)/WP)
  if ratio > 1.0:
    return 1.0
  return ratio

#init(r,xi,Er,Ez)
def plotIonizationRegion(fname,t0):
  global eRatio, maxRpos, maxZpos
  max_r = 0
  maxRatio= 0
  maxXi = 0
  WDT  = np.zeros((len(r),len(xi)))
  for j in range(len(xi)):
    for i in range(int(len(r)/2.0), int(3*len(r)/4)):
            #print('Row ',i, "/",len(r),", Column ", j,"/",len(xi), end="\r", flush=True)
      ratio =  Wdt(i,j)
      if ratio > maxRatio and xi[j] > 3:
        maxRatio = ratio
        maxXi = xi[j]
      if ratio > .01:
        if r[i] > max_r and xi[j] > 3:
          max_r = r[i]
        #print(ratio)
        WDT[i][j] = ratio

  fig, ax  = plt.subplots()#2,figsize=(9,6),sharex=True)
  E = Er 
  eRatio = np.ma.masked_where(eRatio == 0, eRatio)
  WDT = np.ma.masked_where(WDT == 0, WDT)
  cmap = plt.cm.OrRd
  cmap.set_bad(color = (1,1,1,0))
  

  #for ax in axs:
  colors = ax.pcolormesh(xi ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-4,vmax=4),cmap="RdBu_r")
  
  #colors2 = ax.pcolormesh(xi ,r,eRatio, cmap=cmap)
  colors3 = ax.pcolormesh(xi ,r,WDT,vmin=0.,vmax=1.0, cmap=cmap)
  #tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  #cbar = fig.colorbar(colors,ax=axs[0],ticks=tick_locations, format=ticker.LogFormatterMathtext())
  #cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')
  #cbar2 = fig.colorbar(colors2,ax=axs[0])
  #cbar2.set_label('$N_e/N_0$, Integrated Fraction of Ionized Atoms')
  cbar3 = fig.colorbar(colors3,ax=ax)
  cbar3.set_label('$W\Delta t$, Fraction of Ionized Atoms')

  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Ionization Region, t='+str(round(t0,1))+"$\omega_p^{-1}$")
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  plt.text(1,4,"Maximum Ionization Fraction = "+str(round(maxRatio,3))+" at $\\xi = $"+str(round(maxXi,2))) 
  #plt.xlim(xi[0], xi[-1])
  #plt.xlim(5,9)
  ax.set_ylim(0,r[-1])
  #axs[1].set_ylim(0,r[-1])
  #axs[0].set_xlim(3.75,4.75)#r[-1])
  #axs[1].set_xlim(3.75,4.75)#r[-1])
  fn = "plots/"+fname+"_ionizationRegion.png"
  plt.savefig(fn,dpi=200)
  #plt.show()
  return maxXi, max_r
#plotIonizationRegion()
def plotIonizationAtRadius():
        #dat = np.load('data.npz')
  #xi = dat['xi']
  #escaped = dat['esc']
  #trail = []#dat['beam']
  #drive = []
  eRatio  = np.zeros(len(r))
  #n = 0
  for j in range(len(xi)):
    if abs(xi[j] - 4) <0.01:
      print("Found xi=4")
      for i in range(int(len(r)/2.0), len(r)):
            #if escaped[i,j] == 1:
#        ri = r[i]
#        xii = xi[j]
#        xif = dat['beam'][j]
#        drive.append(xi[j])
#        trail.append(dat['beam'][j])
#        n += 1
        print('Row ',i, "/",len(r),", Column ", j,"/",len(xi), end="\r", flush=True)
      #if Er[i,j] < -0.5:
        ratio =  ionRatio(i,j)
        #if ratio > .1:
        #              print(ratio)
        eRatio[i] = ratio
  
  fig, ax  = plt.subplots(figsize=(9,6))
  E = Er 
  #eRatio = np.ma.masked_where(eRatio == 0, eRatio)
  #cmap = plt.cm.OrRd
  #cmap.set_bad(color = (1,1,1,0))

  #colors = ax.pcolormesh(xi ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  #colors2 = ax.pcolormesh(xi ,r,eRatio, cmap=cmap)
  #tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  #cbar = fig.colorbar(colors,ax=ax,ticks=tick_locations, format=ticker.LogFormatterMathtext())
  #cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')
  #cbar2 = fig.colorbar(colors2,ax=ax)
  #cbar2.set_label('$N_e/N_0$, Fraction of Ionized Atoms')
  ax.plot(r, eRatio)
  ax.set_xlabel('r ($c/\omega_p$)')
  ax.set_title('Ionization Fraction at $\\xi=4$')
  ax.set_ylabel("Ionization Fraction")
  
  #plt.xlim(xi[0], xi[-1])
  #plt.xlim(5,9)
  plt.xlim(0,.25)#r[-1])
  fn = "ionizationRatiosvRadius.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()
  return
#plotIonizationAtRadius()
