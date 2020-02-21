#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm
import include.plotSimTracks as plotSimTracks

def plot(t,dvi, dvj, dviRatio, dvjRatio):
  plt.style.use('seaborn-poster')

  fig, ax = plt.subplots()
    
  ax.plot(t,dviRatio,'r',label = "$dv_\\xi$ Error")
  ax.plot(t,dvjRatio,'b',label = "$dv_r$ Error")

  ax.set_xlabel("$\\xi$ $(c/\omega_p)$")
  ax.set_ylabel("$(dv_i^0 - dv_i)/dv_i$")
  
  ax.set_title("Error in Numerical Approximation for Velocity")

  ax.legend()
  fn = "plots/error.png"
  plt.savefig(fn,transparent=True)
  plt.show()
