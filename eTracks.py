# Electron Tracking Script
# Author: Audrey Claire Farrell - audrey.farrell@stonybrook.edu
#   This script is designed to track the path of electrons injected in
#   electron-driven plasma wakefield accelerators.

import sys
import numpy as np
import matplotlib as plot

#Definition of Constants
M_E = 9.109e-31             #electron rest mass in kg
E = 1.60217662e-19          #electron charge in C
EP_0 = 8.854187817          #vacuum permittivity in C/(V m)
C = 299892458               #speed of light in vacuum in m/s
N = 1e23                    #electron number density in 1/m^3
W = np.sqrt(N*E**2/(M_E*EP_0))  #plasma frequency in 1/s

#Code Structure:
# -Function for determining transverse E Field
# -Function for determining radial position at time t
#     -Takes input of present time position (t, r) and time step dt
#     -Returns position at t+dt
# -Function for plotting the position as a function of xi, on a color
#   axis of E field
# -Function for creating animated gif from resultant images

def E(r):
        return -1.0*W**2*M_E/(4*E) * r

