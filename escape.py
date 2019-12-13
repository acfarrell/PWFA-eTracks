import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col

import include.eTracks as eTracks 
import include.getOsirisFields as osiris

Er = osiris.transE()
r,z,t0 = osiris.axes()

def plot(r,z,t0,escaped):
  plt.style.use('seaborn-poster')
  fig, ax  = plt.subplots()
  E = osiris.transE()

  #define binary color map
  ternary_cmaplist = [(1.0,0.,0.,1.0),(0.,0.,0.,0.0),(0.0,1.0,0.0,1.0)]
  ternary_cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', ternary_cmaplist, 3)
  
  colors = ax.pcolormesh(z - t0 ,r,E,norm=col.SymLogNorm(linthresh=0.03,linscale=0.03,vmin=-E.max(),vmax=E.max()),cmap="RdBu_r")
  
  colors2 = ax.pcolormesh(z - t0 ,r,escaped,cmap=ternary_cmap)
  
  cbar = fig.colorbar(colors,ax=ax)
  cbar.set_label('Transverse Electric Field ($m_e c\omega_p / e$)')
  
  cbar2 = fig.colorbar(colors2,ax=ax)
  cbar2.set_ticks([])
  cbar2.set_label('Bound, Untested, Escaped')

  ax.set_xlabel("$\\xi$ ($c/\omega_p$)")
  ax.set_ylabel('r ($c/\omega_p$)')
  ax.set_title('Captured and Escaped Electron Test')
    
  
  plt.xlim(z[0]- t0, z[-1]-t0)
  fn = "escaped.png"
  plt.savefig(fn,transparent=True)
  plt.show()

  return


def main():
  nrows = len(Er)
  ncols = len(Er[0])
  escaped = np.zeros((nrows,ncols))
  for i in range(len(Er) -  1):
    j = int(len(Er[0])/2)
    while j < ncols:
      if Er[i,j] < 0:
        escaped[i,j] = eTracks.GetTrajectory(r[i],0,0,z[j],0,0,False)
      j += 1
  plot(r,z,t0,escaped)
main()
