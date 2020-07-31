 
#Script for generating plots of electron trajectories

import numpy as np
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import include.plotSimTracks as plotSimTracks

def plot(r, xi, pr, E, r_sim, xi_sim, trackName):

  fig, axs = plt.subplots(2,1)
  #Make color axis of electric field
  colors = axs[0].pcolormesh(xi_sim ,r_sim,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
    
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  #cbar = fig.colorbar(colors,ax=axs[0],ticks=tick_locations, format=ticker.LogFormatterMathtext())
  #cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)') 
  axs[0].set_xlabel("$\\xi$ ($c/\omega_p$)")
  axs[0].set_ylabel('r ($c/\omega_p$)')
  axs[0].set_title('Transverse Phase Space for ' + trackName + ' Track')
    
  axs[0].plot(xi,r,'r',label = "$r(\\xi)$")
  axs[0].set_xlim(xi_sim[0], xi_sim[-1])
  axs[0].set_ylim(0,5)
  ax2 = axs[0].twinx()
  ax2.plot(xi,pr,'b',label = "$p_r(\\xi)$")
  ax2.set_ylabel(
  "$p_r$, Transverse Momentum ($m_e c$)")
  axs[0].legend()
  ax2.legend(loc='upper left')

# Second plot of transverse phase space
  axs[1].plot(r,pr,'k')
  axs[1].set_xlabel('r ($c/\omega_p$)')
  axs[1].set_ylabel('$p_r$ ($m_e c$)')

  fig.set_size_inches(10,8)
  fn = "plots/phaseSpace_"+trackName +".png"
  plt.savefig(fn,dpi=200,transparent=True)
  #plt.show()

def plotAll(data):
  fig, ax = plt.subplots()
  colors = pl.cm.jet(np.linspace(0,1,len(data)))

# Second plot of transverse phase space
  for i in range(len(data)):
    track = data[i]
    r = track[0]
    pr = track[-1]
    ax.plot(r,pr,color=colors[-i - 1])
    
  ax.set_xlabel('r ($c/\omega_p$)')
  ax.set_ylabel('$p_r$ ($m_e c$)')
  ax.set_title('Transverse Phase Space for Electrons at Various Initial Radii')
  fig.set_size_inches(10,8)
  fn = "plots/phaseSpace_allTracks.png"
  plt.savefig(fn,dpi=200,transparent=True)
  plt.show()
