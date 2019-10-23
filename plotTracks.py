#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm

def plot(r,z,t,E,r_sim,z_sim,SHM):
  plt.style.use('seaborn-poster')

  fig, ax = plt.subplots()
  #Make color axis of electric field
  if SHM:
    colors = ax.pcolormesh(z,r,E,cmap="YlGnBu")
    cbar = fig.colorbar(colors,ax=ax)
    cbar.set_label('Electric Field ($m_e c\omega_p/e$)')
    ax.set_xlabel("z ($c/\omega_p$)")
    ax.set_ylabel("r ($c/\omega_p$)")
    ax.plot(z,r,'k')
  else:
    colors = ax.pcolormesh(z_sim,r_sim,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
    
    cbar = fig.colorbar(colors,ax=ax)
    cbar.set_label('Electric Field ($m_e c\omega_p / e$)')
    ax.set_xlabel('z ($c/\omega_p$)')
    ax.set_ylabel('r ($c/\omega_p$)')
    ax.set_title('Transverse Electric Field from Simulation')
    
    ax.plot(z,r,'k')

  ax.set_title("Electron Radial Trajectory")
  plt.show()
