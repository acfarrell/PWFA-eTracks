import sys
import math
import csv
import time
import numpy as np
from tempfile import TemporaryFile as tmp

import include.eTracks as eTracks 
import include.getOsirisFields as osiris

Er = osiris.transE()
r,xi,t0 = osiris.axes()

def find_nearest_index(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def main():
  nrows = len(Er)
  ncols = len(Er[0])
  startIdx_r = find_nearest_index(r,0.1)
  endIdx_r = find_nearest_index(r,0.2)
  startIdx_xi = find_nearest_index(xi,1.2)
  endIdx_xi = find_nearest_index(xi,1.4)
  escaped = np.zeros((nrows,ncols))

  driveBeamProf = []
  trailBeamProf = []
 
  gamma = 50/.511
  pz0 = math.sqrt(gamma**2 - 1)
  
  fname = 'backtrackUniform.npz'
  eCount = 0
  num = 0
  for i in range(startIdx_r, endIdx_r):
    #if i == 1:
      #continue
    for j in range(startIdx_xi, endIdx_xi):
      eCount += 1
      
      escaped[i,j], xiPos = eTracks.GetTrajectory(r[i],0.0001,xi[j],pz0)
      if escaped[i,j] == 1:
        trailBeamProf.append(xiPos)
      print('Row ',i, "/",nrows,", Column ", int(j - ncols/2),"/",int(ncols/2)," : ",eCount," electrons", end="\r", flush=True)
  np.savez(fname, r=r,xi=xi, esc=escaped, beam=trailBeamProf)
  eTracks.plotTest()
  print("\n Finished ")
start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))
#plot()
