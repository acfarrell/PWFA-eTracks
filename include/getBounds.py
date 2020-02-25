# Iterative test of boundary conditions for electrons inside the wake
# of plasma wakefield accelerators.

import numpy as np
import math
#import eTracks
from .getOsirisFields import axes, transE
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import seaborn as sea
import matplotlib.cm as cm

def getBounds():
  Er = transE()
  r,xi,t0 = axes()
  thresh = 0.#00001
  bounds = np.zeros((len(Er), len(Er[0])),dtype=int)
  for j in range(len(Er[0]) - 1):
    col = Er[:,j]
    
    maxE = np.amax(col)
    maxE_dex = np.where(col == maxE)[0]
    if maxE < 0.1 or r[maxE_dex] < 1.0:
      continue
    i = maxE_dex
    while i < len(col) - 1:
      bounds[i,j] = 1
      i+= 1

  return bounds

def plotBounds():
  plt.style.use('seaborn-poster')
  fig, ax  = plt.subplots()
  r, z, t0 = axes()
  E = transE()
  
  bounds = getBounds()

  #define binary color map
  binary_cmaplist = [(0.,0.,0.,0.),(0.,0.,0.,.3)]
  binary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', binary_cmaplist, 2)
  colors = ax.pcolormesh(z - t0 ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  colors2 = ax.pcolormesh(z - t0 ,r,bounds,cmap=binary_cmap)
  cbar = fig.colorbar(colors,ax=ax)
  cbar.set_label('Transverse Electric Field ($m_e c\omega_p / e$)')
  
  cbar2 = fig.colorbar(colors2,ax=ax)
  cbar2.set_ticks([])
  cbar2.set_label('Allowed/Not Allowed')

  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Allowed Region for Boundary Conditions')
    
  
  plt.xlim(z[0]- t0, z[-1]-t0)
  fn = "bounds.png"
  plt.savefig(fn,transparent=True)
  plt.show()


