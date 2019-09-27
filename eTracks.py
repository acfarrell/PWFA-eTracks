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

def EField(r):
        #returns the electric field at position r (from Wei Lu's paper)
        return -1.0/4.0 * r

def Momentum(r, dt, p_0):
        #Returns the momentum at t + dt, in units of m_e c
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


def main():
        #Get initial position and momentum from user input:
        r_0 = float(input("Initial radius (c/w_p): "))
        p_0 = float(input("Initial transverse momentum (m_e c): "))
        
        #Determine trajectory, creates n-length lists of data points
        r_dat, t_dat = GetTrajectory(r_0,p_0)
        #Get number of data points
        n = len(r_dat)
        #Create list of E values for each transverse position
        E_dat = [[EField(r) for r in r_dat] for t in t_dat] 
        #Plot the trajectory
        plotTracks.plot(r_dat, t_dat, E_dat)
main()
