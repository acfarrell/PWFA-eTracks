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


# Initialize simulation information from user input as global variables
if len(sys.argv) == 2:
  input_fname = str(sys.argv[1])
else:
  input_fname = "input.example"
print("Using initial conditions from ",input_fname)
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

Er = quickPIC.transE(Er_fname)
Ez = quickPIC.longE(Ez_fname)
Bphi = quickPIC.longE(Bphi_fname)

r,xi = quickPIC.axes(Er_fname)

def W(E0):
  # Takes an electric field in normalized units and returns the ionization rate in s^-1
  # Define constants for calculating ionization rate for He
  gsE = 24.5 # unperturbed ground state energy in eV
  Z = 1 # charge number after ionization
  n = 0.746 # effective principal quantum number
  
  E = E0 * (M_E * C * WP)/EC * 1e-9  #Convert normalized field to GV/m

  W0 = 1.52e15 * 4**n * gsE / (n * gamma(2*n)) * (20.5*gsE**(3/2))**(2*n-1)
  
  W = W0/((E)**(2*n-1)) * math.exp(-6.83*gsE**(3/2)/E)
  return W 

def ionRatio(i,j):
  dt = 14.7/WP
  #integral = 0
  #for n in range(1,j):
  En = math.sqrt((Er[i,j])**2 + (Ez[i,j])**2)
  #  integral = integral +  W(En) * (xi[n] - xi[n-1])/WP
  ratio = W(En) * dt#1 - math.exp(-1*integral)
  if ratio > 1.0:
    return 1.0
  return ratio


def main():
  eTracks.InitFields(Er,Ez,Bphi,r,xi,t0)
  nrows = len(r)
  ncols = len(xi)
  escaped = np.zeros((nrows,ncols))

  driveBeamProf = []
  trailBeamProf = []
  
  eCount = 0
  num = 0
  ionizedCount = 0

  eTot = 50#int(injCharge / EC) #number of electrons to inject print("Simulating ",eTot," injected electrons.")
  for j in range(len(xi)):
    for i in range(int(len(r)/2)):
      ratio = ionRatio(i,j)
      if ratio > 0.1:
        ionizedCount += 1
    
  while eCount < eTot: 
    if ionizedCount ==0:
      break
    i = rand.randint(int(nrows/2),nrows-1)
    j = rand.randint(int(ncols/4),ncols-1)
    
    ratio = ionRatio(i,j) 
    if ratio > 0.1:
      prob = rand.uniform(0,1)
      if prob < ratio:
        eCount += 1
        if eCount > eTot:
          break
        escaped[i,j], xiPos = eTracks.GetTrajectory(r[i],xi[j])
        if escaped[i,j] == 1 or escaped[i,j] == 2 :
          trailBeamProf.append(xiPos)
          driveBeamProf.append(xi[j])
        print('Row ',i, "/",nrows,", Column ", int(j - ncols/2),"/",int(ncols/2)," : ",eCount," electrons", end="\r", flush=True)
  np.savez(output_fname, r=r,xi=xi, esc=escaped, drive=driveBeamProf,trail=trailBeamProf)
  eTracks.plotTest(output_fname,t0)
  eTracks.plotVzTest(output_fname,t0)
  #plot()
  print("\n Finished ")

  

start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))
