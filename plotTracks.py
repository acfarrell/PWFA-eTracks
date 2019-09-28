#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm

def plot(r,t,E):
  plt.style.use('seaborn-poster')

  fig, ax = plt.subplots(figsize=(10,10))
  #Make color axis of electric field
  colors = ax.contourf(t,r,E,levels=100,cmap="YlGnBu",alpha=0.5)
  #Plot radial position v. time
  ax.plot(t,r,'k')
  
  #Labels
  cbar = fig.colorbar(colors,ax=ax)
  cbar.set_label('Electric Field ($m_e c\omega_p/e$)')
  ax.set_xlabel("Time ($\omega_p^{-1}$)")
  ax.set_ylabel("Radius ($c/\omega_p$)")
  ax.set_title("Electron Radial Trajectory")
  plt.show()
