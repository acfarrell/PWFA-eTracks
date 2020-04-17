import sys
import csv
import time
import numpy as np
from tempfile import TemporaryFile as tmp

import include.eTracks as eTracks 
import include.getOsirisFields as osiris

Er = osiris.transE()
r,xi,t0 = osiris.axes()

def main():
  nrows = len(Er)
  ncols = len(Er[0])
  escaped = np.zeros((nrows,ncols))

  driveBeamProf = []
  trailBeamProf = []
  
  fname = 'data.npz'
  eCount = 0
  num = 0
  for i in range(nrows):
    #if i == 1:
      #continue
    for j in range(int(ncols/2) , ncols):
      if Er[i,j] < -0.5:
        eCount += 1
        plotTrack = False
        if eCount % 1000 == 0:
          plotTrack = True
          num +=1
        
        escaped[i,j], xiPos = eTracks.GetTrajectory(r[i],0,0,xi[j],0,0,plotTrack, num)
        if escaped[i,j] == 1:
          trailBeamProf.append(xiPos)
      print('Row ',i, "/",nrows,", Column ", int(j - ncols/2),"/",int(ncols/2)," : ",eCount," electrons", end="\r", flush=True)
  np.savez(fname, r=r,xi=xi, esc=escaped, beam=trailBeamProf)
  #plot()
  print("\n Finished ")
start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))
#plot()
