#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import include.plotSimTracks as plotSimTracks

def plot(r, z, t, xi, E, r_sim, z_sim, SHM, track):

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
    colors = ax.pcolormesh(z_sim - 858.95 ,r_sim,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
    tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
    cbar = fig.colorbar(colors,ax=ax,ticks=tick_locations, format=ticker.LogFormatterMathtext())
    cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')    
    ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
    ax.set_ylabel('r ($c/\omega_p$)')
    ax.set_title('Electron trajectory from simulation v. OSIRIS for ' + track + ' track')
    
    ax.plot(xi,r,'k',label = "Simulated")

  plt.xlim(z_sim[0]- 858.95, z_sim[-1]-858.95)
  xi_OSIRIS, r_OSIRIS = plotSimTracks.get_xir(track)
  ax.plot(xi_OSIRIS, r_OSIRIS, 'r:', label="OSIRIS")
  ax.legend()
  if SHM:
    model = "SHM"
  else:
    model = "simE"
  fn = "plots/"+model + "_"+track +".png"
  plt.savefig(fn,transparent=True)
  #plt.show()
