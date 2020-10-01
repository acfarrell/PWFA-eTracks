import sys
import math
import csv
import time
import importlib
import numpy as np
from tempfile import TemporaryFile as tmp
import random as rand
from scipy.special import gamma
import include.eTracks as eTracks 
import include.getOsirisFields as osiris
import include.getQuickPICFields as quickPIC
import include.plotIonizationRegion as ionization


# Initialize simulation information from user input as global variables
if len(sys.argv) == 2:
  input_fname = str(sys.argv[1])
else:
  input_fname = "input.example"
#print("Using initial conditions from ",input_fname)
init = importlib.import_module(input_fname)

#Output File Name
output_fname = init.output_fname

#Field Map File Names
Er_fname = init.Er_fname
Ez_fname = init.Ez_fname
Bphi_fname = init.Bphi_fname
t0 = init.t0

#Beam Parameters
sig_r = init.sigma_r #transverse beam spread in c/wp
sig_z = init.sigma_z #longitudinal beam spread in c/wp

#Impurity Profile
n_imp = init.n_imp #longitudinal helium density profile in np

#Total amount of charge to inject
injCharge = init.injCharge

#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817e-12                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
NP = 1e23                          #electron number density in 1/m^3
WP = np.sqrt(NP*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

#Er = quickPIC.transE(Er_fname)
#Ez = quickPIC.longE(Ez_fname)
#Bphi = quickPIC.longE(Bphi_fname)

r,xi, Er = quickPIC.spliceLowRes(Er_fname)
r,xi, Ez = quickPIC.spliceLowRes(Ez_fname)
r,xi, Bphi = quickPIC.spliceLowRes(Bphi_fname)
dy = quickPIC.getDepth()



def main():
  pinchXi, pinchR = ionization.init(r,xi,Er,Ez,t0,output_fname)
  #ionization.plotIonizationRegion()
  eTracks.InitFields(Er,Ez,Bphi,r,xi,t0)
  nrows = len(r)
  ncols = len(xi)
  escaped = np.zeros((nrows,ncols))

  driveBeamProf = []
  trailBeamProf = []
  
  eCount = 0
  num = 0
  ionizedCount = 0
  totIonFrac = 0
  capturedCount = 0
  print("Simulating t = ", t0)
  eTot = int(injCharge / EC) #number of electrons to inject

  #print("Simulating ",eTot," injected electrons.")
  for j in range(int(ncols/4),ncols -1):
    for i in range(int(nrows/2),nrows -1):
      ratio = ionization.Wdt(i,j)
      if ratio > 0.01:
        totIonFrac += ratio
        ionizedCount += 1
  print(ionizedCount," pixels in ionization region.")
  if ionizedCount > 0:
    meanIonFrac = totIonFrac / ionizedCount
    #Calculate volume of each pixel in SI units
    pixVolume = dy * (xi[1] - xi[0]) * (r[1] - r[0]) * C / WP
    pinchVolume = pixVolume * ionizedCount
    # Calculate impurity density 
    impDensity = (eTot / pinchVolume)/ (meanIonFrac * .6641)
    print("Impurity density = ", impDensity," m^-3")
    
  #while eCount < eTot: 
  else:
    return
  for j in range(int(ncols/4),ncols -1):
    for i in range(int(nrows/2),nrows -1):
            #i = rand.randint(int(nrows/2),nrows-1)
    #j = rand.randint(int(ncols/4),ncols-1)
    #r0 = 0.1#0.01 + (pinchR / eTot) * (eCount+1)
      ratio = ionization.Wdt(i,j) 
      if ratio > 0.01:
        ionCount = int(ratio * impDensity * pixVolume)
        escaped[i,j], xiPos = eTracks.GetTrajectory(r[i],xi[j],True)
        for electron in range(ionCount):
          print('Row ',i, "/",nrows,", Column ", int(j - ncols/2),"/",int(ncols/2)," : ",eCount," electrons", end="\r", flush=True)
          eCount += 1
        #if eCount > eTot:
          #break
          #escaped[i,j], xiPos = eTracks.GetTrajectory(r[i],xi[j],True)
  #    escaped[i,j], xiPos = eTracks.GetTrajectory(r0,pinchXi,False)#r[i],xi[j])
          if escaped[i,j] == 1 or escaped[i,j] == 2 :
            capturedCount += 1
            trailBeamProf.append(xiPos)
            driveBeamProf.append(xi[j])
        #print('Row ',i, "/",nrows,", Column ", int(j - ncols/2),"/",int(ncols/2)," : ",eCount," electrons", end="\r", flush=True)
  np.savez(output_fname, r=r,xi=xi, esc=escaped, drive=driveBeamProf,trail=trailBeamProf)
  percentCaptured = 100 * capturedCount / eCount
  eTracks.plotTest(output_fname,t0,percentCaptured)
  eTracks.plotMomentaTest(output_fname,t0,percentCaptured)
  #eTracks.plotNoBTest(output_fname,t0)
  #plot()
  print("\n Finished ")

  

start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))
