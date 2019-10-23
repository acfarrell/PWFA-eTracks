#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm

def plot(r,z,t,E,r_sim,z_sim,SHM):
  plt.style.use('seaborn-poster')

  fig, ax = plt.subplots(figsize=(10,10))
  #Make color axis of electric field
  if SHM:
    colors = ax.contourf(t,r,E,levels=100,cmap="YlGnBu",alpha=0.5)
    ax.plot(t,r,'k')
  else:
    #colors = ax.pcolormesh(z_sim,r_sim,np.transpose(E),norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=E.min(),vmax=E.max()),cmap="RdBu_r")
    #cbar = fig.colorbar(colors,ax=ax)
  
    #cbar.set_label('Electric Field ($m_e c\omega_p / e$)')
    ax.set_xlabel('z ($c/\omega_p$)')
    ax.set_ylabel('r ($c/\omega_p$)')
    ax.set_title('Transverse Electric Field from Simulation')
    ax.plot(z,r,'k')
  #Plot radial position v. time
  
  #Labels
  #cbar = fig.colorbar(colors,ax=ax)
  #cbar.set_label('Electric Field ($m_e c\omega_p/e$)')
  #ax.set_xlabel("Time ($\omega_p^{-1}$)")
  #ax.set_ylabel("Radius ($c/\omega_p$)")
  ax.set_title("Electron Radial Trajectory")
  plt.show()
