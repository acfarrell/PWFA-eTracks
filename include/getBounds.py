# Iterative test of boundary conditions for electrons inside the wake
# of plasma wakefield accelerators.

import numpy as np
import math
#import eTracks
from .getOsirisFields import axes, transE
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import matplotlib.ticker as ticker

def getBounds(Er,r,xi,t0):
  thresh = 0.#00001
  bounds = np.zeros((len(Er), len(Er[0])),dtype=int)
  for j in range(len(Er[0]) - 1):
    col = Er[:,j]
    #print(col) 
    maxE = np.amax(col)
    #print(maxE)
    maxE_dex = (np.where(col == maxE))[0]
    #print(maxE_dex[0])
    if maxE < 0.1 or r[maxE_dex[0]] < 1.0:
      continue
    i = maxE_dex[0]
    while i < len(col) - 1:
      bounds[i,j] = 1
      i+= 1

  return bounds

def plotBounds():
  fig, ax  = plt.subplots()
  r, z, t0 = axes()
  E = transE()
  
  bounds = getBounds()

  #define binary color map
  binary_cmaplist = [(0.,0.,0.,0.),(0.,0.,0.,.3)]
  binary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', binary_cmaplist, 2)
  colors = ax.pcolormesh(z - t0 ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  tick_locations=[x*0.01 for x in range(2,10)]+ [x*0.01 for x in range(-10,-1)] + [x*0.1 for x in range(-10,10)] +[ x for x in range(-10,10)]
  cbar = fig.colorbar(colors,ax=ax,ticks=tick_locations, format=ticker.LogFormatterMathtext())
  cbar.set_label('$E_r$, Transverse Electric Field ($m_e c\omega_p / e$)')


  colors2 = ax.pcolormesh(z - t0 ,r,bounds,cmap=binary_cmap)
  
  cbar2 = fig.colorbar(colors2,ax=ax)
  cbar2.set_ticks([])
  cbar2.set_label('Allowed/Not Allowed')

  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Allowed Region for Boundary Conditions')
    
  
  plt.xlim(z[0]- t0, z[-1]-t0)
  fn = "bounds.png"
  plt.savefig(fn,dpi=300,transparent=True)
  plt.show()

#plotBounds()
