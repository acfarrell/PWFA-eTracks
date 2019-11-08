#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sea
import matplotlib.cm as cm
import include.plotSimTracks as plotSimTracks

def plot(r, z, t, xi, E, r_sim, z_sim, SHM, track):
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
    colors = ax.pcolormesh(z_sim - 858.95 ,r_sim,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
    
    cbar = fig.colorbar(colors,ax=ax)
    cbar.set_label('Transverse Electric Field ($m_e c\omega_p / e$)')
    ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
    ax.set_ylabel('r ($c/\omega_p$)')
    ax.set_title('Electron trajectory from simulation v. OSIRIS for ' + track + ' track')
    
    ax.plot(xi,r,'k',label = "Simulated")

  plt.xlim(z_sim[0]- 858.95, z_sim[-1]-858.95)
  xi_OSIRIS, r_OSIRIS = plotSimTracks.get_xir(track)
  ax.plot(xi_OSIRIS, r_OSIRIS, 'k--', label="OSIRIS")
  ax.legend()
  if SHM:
    model = "SHM"
  else:
    model = "simE"
  fn = "plots/"+model + "_"+track +".png"
  plt.savefig(fn,transparent=True)
  plt.show()
