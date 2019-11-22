import matplotlib.pyplot as plt
import numpy as np
import include.plotSimTracks as trackData
import include.getOsirisFields as fieldData
import math

def find_nearest_index(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def plot_Residuals(track):
  trackData.get_xir(track)
  
  rField, zField, t0 = fieldData.axes()
  ErField = fieldData.transE()
  EzField = fieldData.longE()
  BphiField = fieldData.phiB()

  ErTrack = trackData.E_r
  EzTrack = trackData.E_z
  BphiTrack = trackData.B

  rTrack = trackData.r
  zTrack = trackData.z
  xiTrack = trackData.xi
  
  fig, ax = plt.subplots(3, sharex=True)
  fig.suptitle('Residual Fields Along OSIRIS Track for '+ track + ' Track')
  ErDiff = []
  EzDiff = []
  BphiDiff = []
  for i in range(len(xiTrack)):
    zDex = find_nearest_index(zField, zTrack[i])
    rDex = find_nearest_index(rField, rTrack[i])
    
    ErDiff.append( ErField[rDex,zDex]- ErTrack[i])
    EzDiff.append( EzField[rDex,zDex] - EzTrack[i])
    BphiDiff.append(BphiField[rDex,zDex] - BphiTrack[i])

  ax[0].plot(xiTrack, ErDiff,'r')
  ax[1].plot(xiTrack, EzDiff,'r')
  ax[2].plot(xiTrack, BphiDiff,'r')
  
  x = np.linspace(0,8,50)
  plt.xlim(0,8)
  ax[0].set_ylim(-10,10)
  ax[1].set_ylim(-10,10)
  ax[2].set_ylim(-10,10)
  plt.xlabel('$\\xi\, (c/\omega_p)$') 
  #Add Dashed line at 0
  ax[0].plot(x,[0.0 for n in range(50)],':k')
  ax[1].plot(x,[0.0 for n in range(50)],':k')
  ax[2].plot(x,[0.0 for n in range(50)],':k')

  ax[0].set_ylabel('$E_r^{data} - E_r^{track}$')
  ax[1].set_ylabel('$E_z^{data} - E_z^{track}$')
  ax[2].set_ylabel('$B_\phi^{data} - B_\phi^{track}$')
  plt.savefig("plots/fieldResiduals_"+track,transparent=True)
  plt.show()

def main():
  track = str(input("Which track? "))
  plot_Residuals(track)
  return
main()
