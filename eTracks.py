# Electron Tracking Script
# Author: Audrey Claire Farrell - audrey.farrell@stonybrook.edu
#   This script is designed to track the path of electrons injected in
#   electron-driven plasma wakefield accelerators.

import sys
import numpy as np
import plotTracks 

#Definition of Constants
M_E = 9.109e-31                   #electron rest mass in kg
EC = 1.60217662e-19               #electron charge in C
EP_0 = 8.854187817                #vacuum permittivity in C/(V m)
C = 299892458                     #speed of light in vacuum in m/s
N = 1e23                          #electron number density in 1/m^3
W = np.sqrt(N*EC**2/(M_E*EP_0))   #plasma frequency in 1/s

#Code Structure:
# -Function for determining transverse E Field
# -Function for determining radial position at time t
#     -Takes input of present time position (t, r) and time step dt
#     -Returns position at t+dt
# -Function for plotting the position as a function of xi, on a color
#   axis of E field
# -Function for creating animated gif from resultant images

def EField(r):
        #returns the electric field at position r (from Wei Lu's paper)
        return -1.0*W*M_E/(4*EC) * r

def Momentum(r, dt, p_0):
        #Returns the momentum at t + dt
        return p_0 + EC * EField(r) * dt

def Velocity(p):
        #returns the velocity from the momentum
        return p/M_E

def GetTrajectory(r_0,p_0):
        #returns array of r v. t
        


def main():
        #create nxn array of xi v. r, with values of the electric field
        n = 20
        r = np.linspace(0,1,n)
        xi = np.linspace(0,1,n)
        E_dat = 20 * [0]
        
        for i in range(0,20):
                E_dat[i] = 20 * [EField(r[i])]
        plotTracks.plot(E_dat)
main()
