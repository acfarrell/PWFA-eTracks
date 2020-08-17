 
#Script for generating plots of electron trajectories

import numpy as np
import math
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import include.plotSimTracks as plotSimTracks

from scipy.signal import argrelextrema

def plot(r, xi, pr, E, r_sim, xi_sim, trackName):
  print("Plotting Phase Space and Trajectory")
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
  axs[1].plot(pr,r,'k')
  axs[1].set_ylabel('r ($c/\omega_p$)')
  axs[1].set_xlabel('$p_r$ ($m_e c$)')

  fig.set_size_inches(10,8)
  fn = "plots/phaseSpace_"+trackName +".png"
  plt.savefig(fn,dpi=200,transparent=True)
  plotPhaseSpace(r,pr,trackName)
  plotInteractionRegion(r,xi,pr,E,r_sim,xi_sim,trackName)
  #plt.show()
  plt.close()
def plotPhaseSpace(r,pr, trackName):

  fig, ax = plt.subplots()
  ax.set_title("Transverse Phase Space for "+trackName+", Initial Kick")
# plot of transverse phase space for initial kick from drive beam
  #find first maximum of transverse momentum
  idx = argrelextrema(pr,np.greater)[0][0]
  print("Index of Kick = ",idx)

  ax.plot(pr[:idx],r[:idx],'k')
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_xlabel('$p_r$ ($m_e c$)')
  
  maxp = pr[idx]
  idx9 = find_nearest_index(pr[:idx],0.9*maxp)

  ax.plot(pr[idx],r[idx],'ro',label="$p_{max}$")
  ax.plot(pr[idx9],r[idx9],'bo',label="$0.9p_{max}$")
  ax.axhline(y=r[idx],color='r',linestyle='--')
  ax.axhline(y=r[idx9],color='b',linestyle='--')
  ax.legend()

  fig.set_size_inches(10,8)
  fn = "plots/phaseSpaceKick_"+trackName +".png"
  plt.savefig(fn,dpi=200,transparent=True)
  #plt.show()
  plt.close()
def plotInteractionRegion(r, xi, pr, E, r_sim, xi_sim, trackName):

  idx = argrelextrema(pr,np.greater)[0][0]
  maxp = pr[idx]
  idx9 = find_nearest_index(pr[:idx],0.9*maxp)
  
  fig, ax = plt.subplots()
  #Make color axis of electric field
  colors = ax.pcolormesh(xi_sim ,r_sim,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
    
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  #cbar = fig.colorbar(colors,ax=axs[0],ticks=tick_locations, format=ticker.LogFormatterMathtext())
  #cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)') 
  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Interaction Region for ' + trackName + ' Track')
    
  ax.plot(xi,r,'k',label = "$r(\\xi)$")
  ax.set_xlim(0,8)
  ax.set_ylim(0,3)
  ax.plot(xi[idx],r[idx],'ro',label="($\\xi,r$) where $p=p_{max}$")
  ax.plot(xi[idx9],r[idx9],'bo',label="($\\xi,r$) where $p = 0.9p_{max}$")
  ax.axvline(x=xi[idx],color='r',linestyle='--')
  ax.axvline(x=xi[idx9],color='b',linestyle='--')
  ax.axhline(y=r[idx],color='r',linestyle='--')
  ax.axhline(y=r[idx9],color='b',linestyle='--')
  ax.legend()
  
  fig.set_size_inches(10,8)
  fn = "plots/interactionRegion_"+trackName +".png"
  plt.savefig(fn,dpi=200,transparent=True)
  plt.close()
  #plt.show()
def find_nearest_index(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
      return idx
def plotAll(data):
  fig, ax = plt.subplots()
  colors = pl.cm.jet(np.linspace(0,1,len(data)))

# Second plot of transverse phase space
  for i in range(len(data)):
    track = data[i]
    r = track[0]
    xi = track[1]
    pr = track[2]
    idx = argrelextrema(pr,np.greater)[0][0]
    mindiff = 1
    for ix in range(len(xi)):
      diff = abs(xi[0] - xi[ix] - 1.5)
      if diff < mindiff:
        mindiff = diff
        idx15 = ix
    #idx15 = find_nearest_index(xi,3.5)
    print(idx15)
    ax.plot(r[:idx15],pr[:idx15],color=colors[-i - 1])
    
  ax.set_xlabel('$r$ ($c/\omega_p$)')
  ax.set_ylabel('$p_r$ ($m_e c$)')
  ax.set_title('Transverse Phase Space for Electrons at Various Initial Radii')
  fig.set_size_inches(10,8)
  fn = "plots/phaseSpace_allTracks.png"
  plt.savefig(fn,dpi=200,transparent=True)
  plt.show()
  plt.close()
